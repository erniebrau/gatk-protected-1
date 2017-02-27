#
# A workflow for creating a Panel of Normals given a list of normal samples. Supports both WGS and WES samples. This was tested on a3c7368 commit of gatk-protected.
#
# Notes:
#
# - Input file (normal_bam_list) must contain file paths to bam and bam index files separated by tabs in the following format:
#   normal_bam_1	bam_idx_1
#   normal_bam_2	bam_idx_2
#   ...
# #
# - target file (which must be in tsv format) is only used with WES workflow - WGS workflow generates its own targets (so user can pass any string as an argument)
#
# - THIS SCRIPT SHOULD BE CONSIDERED OF "BETA" QUALITY
#
# - Example invocation:
#    java -jar cromwell.jar run pon_gatk_workflow.wdl myParameters.json
# - See pon_gatk_workflow_template.json for a template json file to modify with your own parameters (please save
#    your modified version with a different filename and do not commit to gatk-protected repo).
#
# - Some call inputs might seem out of place - consult with the comments in task definitions for details
#
#############

workflow pon_gatk_workflow {
    # Workflow input files
    File normal_bam_list
    Array[Array[String]] bam_file_names = read_tsv(normal_bam_list)
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    File gatk_jar

    # Workflow options
    Int wgsBinSize
    String combined_entity_id = "combined_read_counts"
    Boolean noQC

    # PoN name
    String pon_entity_id

  scatter (row in bam_file_names) {

    call GetBamFileName {
      input:
          input_bam=row[0]
    }


    call WholeGenomeCoverage {
      input:
          entity_id=GetBamFileName.name,
          input_bam=row[0],
          bam_idx=row[1],
          ref_fasta=ref_fasta,
          ref_fasta_fai=ref_fasta_fai,
          ref_fasta_dict=ref_fasta_dict,
          gatk_jar=gatk_jar,
          wgsBinSize=wgsBinSize,
          mem=10
    }
  }

  call AggregateTargetFiles {
    input:
        target_file_list=WholeGenomeCoverage.gatk_target_file
  }

  call CombineReadCounts {
    input:
        combined_entity_id=combined_entity_id,
        coverage_file_list=WholeGenomeCoverage.gatk_coverage_file,
        gatk_jar=gatk_jar,
        mem=4
  }

  call AnnotateTargets {
    input:
        entity_id=combined_entity_id,
        gatk_jar=gatk_jar,
        target_file=AggregateTargetFiles.target_file,
        ref_fasta=ref_fasta,
        ref_fasta_fai=ref_fasta_fai,
        ref_fasta_dict=ref_fasta_dict,
        mem=4
  }

  call CorrectGCBias {
    input:
        entity_id=combined_entity_id,
        gatk_jar=gatk_jar,
        coverage_file=CombineReadCounts.combined_read_counts,
        annotated_targets=AnnotateTargets.annotated_targets,
        mem=4
  }

  call CreatePanelOfNormals {
    input:
        pon_entity_id=pon_entity_id,
        read_counts_file=CorrectGCBias.coverage_file_gcbias_corrected,
        gatk_jar=gatk_jar,
        noQC=noQC,
        mem=4
  }
}


# Helper task to get the name of the given bam file
task GetBamFileName {
    File input_bam

    command <<<
        echo $(basename "${input_bam}" .bam)
     >>>

    output {
        String name=read_string(stdout())
    }
}


# Calculate coverage on Whole Genome Sequence using Spark. This task automatically creates an
# output file, containing targets
task WholeGenomeCoverage {
    String entity_id
    File coverage_file 
    File input_bam
    File bam_idx
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    File gatk_jar
    Int wgsBinSize
    Int mem

    command {
        java -Xmx${mem}g -jar ${gatk_jar} SparkGenomeReadCounts --outputFile ${entity_id}.coverage.tsv \
                --reference ${ref_fasta} -keepxy --input ${input_bam} --binsize ${wgsBinSize}
    }

    output {
        File gatk_coverage_file = "${entity_id}.coverage.tsv"
        File gatk_target_file = "${entity_id}.coverage.tsv.targets.tsv"
    }
}

# Helper task to aggregate the array of targets that is output by the scatter block into a single file
task AggregateTargetFiles {
    Array[File]+ target_file_list

    command <<<
        ln -s ${target_file_list[0]} targets.tsv
    >>>

    output {
        File target_file = "targets.tsv"
    }
}

# Combine a set of read-count input files into a single multicolumn output file
task CombineReadCounts {
    String combined_entity_id
    Array[File] coverage_file_list
    File gatk_jar
    Int mem

    command {
        java -Xmx${mem}g -jar ${gatk_jar} CombineReadCounts --input ${sep=" --input " coverage_file_list} \
         --maxOpenFiles 1000 --output ${combined_entity_id}.tsv
    }

    output {
        File combined_read_counts = "${combined_entity_id}.tsv"
    }
}

# Add new columns to an existing target table with various targets
# Note that this task is optional
task AnnotateTargets {
    String entity_id
    File target_file
    File gatk_jar
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    Int mem

    # If GC correction is disabled, then an empty file gets passed downstream
    command {
          java -Xmx${mem}g -jar ${gatk_jar} AnnotateTargets --targets ${target_file} --reference ${ref_fasta} --output ${entity_id}.annotated.tsv
    }

    output {
        File annotated_targets = "${entity_id}.annotated.tsv"
    }
}

# Correct coverage for sample-specific GC bias effects
# Note that this task is optional
task CorrectGCBias {
    String entity_id
    File coverage_file
    File annotated_targets
    File gatk_jar
    Int mem

    # If GC correction is disabled, then the coverage file gets passed downstream unchanged
    command {
          then java -Xmx${mem}g -jar ${gatk_jar} CorrectGCBias --input ${coverage_file} \
           --output ${entity_id}.gc_corrected_coverage.tsv --targets ${annotated_targets}
    }

    output {
        File coverage_file_gcbias_corrected = "${entity_id}.gc_corrected_coverage.tsv"
    }
}

# Create panel of normals given a collection of read counts for control samples
task CreatePanelOfNormals {
    String pon_entity_id
    File read_counts_file
    File gatk_jar
    Boolean noQC
    Int mem

    command {
        # If there are no removed samples the output file still needs to be created
        touch "${pon_entity_id}.pon.removed_samples.txt"
        java -Xmx${mem}g -jar ${gatk_jar} CreatePanelOfNormals --extremeColumnMedianCountPercentileThreshold 2.5 \
         --truncatePercentileThreshold 0.1 --input ${read_counts_file} --output ${pon_entity_id}.pon \
         --noQC ${noQC}
    }

    output {
        File output_pon = "${pon_entity_id}.pon"
        File removed_samples = "${pon_entity_id}.pon.removed_samples.txt"
        File target_weights = "${pon_entity_id}.pon.target_weights.txt"
    }
}
