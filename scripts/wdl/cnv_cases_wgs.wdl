#
# Case sample workflow for a list of pairs of case-control samples. Includes GATK CNV and ACNV. Supports both WGS and WES samples.
#
# Notes:
#
# - the input file(input_bam_list) must contain a list of tab separated values in the following format(one or more lines must be supplied):
# tumor_entity  tumor_bam  tumor_bam_index  normal_entity  normal_bam  normal_bam_index  <--first input
# tumor_entity  tumor_bam  tumor_bam_index  normal_entity  normal_bam  normal_bam_index  <--second input
# etc...
#
#
# - file names will use the entity ID specified, but inside the file, the bam SM tag will typically be used.
#
# - target file (which must be in tsv format) is only used with WES workflow, WGS workflow generates its own targets (so user can pass any string as an argument)
#
# - the WGS PoN must be generated with WGS samples
# 
# - THIS SCRIPT SHOULD BE CONSIDERED OF "BETA" QUALITY
#
# - Example invocation:
#    java -jar cromwell.jar run cnv_cases.wdl myParameters.json
# - See cnv_cases.json for a template json file to modify with your own parameters (please save
#    your modified version with a different filename and do not commit to gatk-protected repo).
#
# - Some call inputs might seem out of place - consult with the comments in task definitions for details
#
#############

workflow wgs_cnv_workflow {
    # Workflow input files
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    File input_bam_list
    Array[Array[String]] bam_list_array = read_tsv(input_bam_list)
    File PoN
    String gatk_jar

    # Workflow output directories and other options
    Int wgsBinSize

  scatter (row in bam_list_array) {

    call WholeGenomeCoverage {
      input:
          entity_id=row[0],
          input_bam=row[1],
          bam_idx=row[2],
          ref_fasta=ref_fasta,
          ref_fasta_fai=ref_fasta_fai,
          ref_fasta_dict=ref_fasta_dict,
          gatk_jar=gatk_jar,
          wgsBinSize=wgsBinSize,
          mem=8
    }

    call AnnotateTargets {
      input:
          entity_id=row[0],
          gatk_jar=gatk_jar,
          target_file=WholeGenomeCoverage.gatk_target_file,
          ref_fasta=ref_fasta,
          ref_fasta_fai=ref_fasta_fai,
          ref_fasta_dict=ref_fasta_dict,
          mem=4
    }

    call CorrectGCBias {
      input:
          entity_id=row[0],
          gatk_jar=gatk_jar,
          coverage_file=WholeGenomeCoverage.gatk_coverage_file,
          annotated_targets=AnnotateTargets.annotated_targets,
          mem=4
    }

    call NormalizeSomaticReadCounts {
      input:
          entity_id=row[0],
          coverage_file=CorrectGCBias.gatk_cnv_coverage_file_gcbias,
          padded_target_file=WholeGenomeCoverage.gatk_target_file,
          pon=PoN,
          gatk_jar=gatk_jar,
          mem=8
    }

    call PerformSegmentation {
      input:
          entity_id=row[0],
          gatk_jar=gatk_jar,
          tn_file=NormalizeSomaticReadCounts.tn_file,
          mem=2
    }

    call PlotSegmentedCopyRatio {
      input:
          entity_id=row[0],
          gatk_jar=gatk_jar,
          tn_file=NormalizeSomaticReadCounts.tn_file,
          pre_tn_file=NormalizeSomaticReadCounts.pre_tn_file,
          segments_file=PerformSegmentation.seg_file,
          ref_fasta_dict=ref_fasta_dict,
          output_dir="plots/${row[0]}/",
          mem=4
    }
  }
}

# Calculate coverage on Whole Genome Sequence using Spark.
# This task automatically creates a target output file.
task WholeGenomeCoverage {
    String entity_id
    File input_bam
    File bam_idx
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    String gatk_jar
    Int wgsBinSize
    Int mem

    command {
        java -Xmx${mem}g -jar ${gatk_jar} SparkGenomeReadCounts --outputFile ${entity_id}.coverage.tsv \
                --reference ${ref_fasta} --input ${input_bam} --sparkMaster local[1] --binsize ${wgsBinSize}
    }

    output {
        File gatk_coverage_file = "${entity_id}.coverage.tsv"
        File gatk_target_file = "${entity_id}.coverage.tsv.targets.tsv"
    }
}

# Add new columns to an existing target table with various targets
# Note that this task is optional 
task AnnotateTargets {
    String entity_id
    File target_file
    String gatk_jar
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
    String gatk_jar
    Int mem

    # If GC correction is disabled, then the coverage file gets passed downstream unchanged
    command {
          java -Xmx${mem}g -jar ${gatk_jar} CorrectGCBias --input ${coverage_file} \
            --output ${entity_id}.gc_corrected_coverage.tsv --targets ${annotated_targets}
    }

    output {
        File gatk_cnv_coverage_file_gcbias = "${entity_id}.gc_corrected_coverage.tsv"
    }
}

# Perform tangent normalization (noise reduction) on the proportional coverage file.
task NormalizeSomaticReadCounts {
    String entity_id
    File coverage_file
    File padded_target_file
    File pon
    String gatk_jar
    Int mem

    command {
        java -Xmx${mem}g -jar ${gatk_jar} NormalizeSomaticReadCounts --input ${coverage_file} \
         --targets ${padded_target_file} --panelOfNormals ${pon} --factorNormalizedOutput ${entity_id}.fnt.tsv --tangentNormalized ${entity_id}.tn.tsv \
         --betaHatsOutput ${entity_id}.betaHats.tsv --preTangentNormalized ${entity_id}.preTN.tsv  --help false --version false --verbosity INFO --QUIET false
    }

    output {
        File tn_file = "${entity_id}.tn.tsv"
        File pre_tn_file = "${entity_id}.preTN.tsv"
        File betahats_file = "${entity_id}.betaHats.tsv"
    }
}

# Segment the tangent normalized coverage profile.
task PerformSegmentation {
    String entity_id
    String gatk_jar
    File tn_file
    Int mem

    command {
        java -Xmx${mem}g -jar ${gatk_jar} PerformCopyRatioSegmentation --tangentNormalized ${tn_file} \
         -S ${entity_id}.seg  -initialNumStates 7
    }

    output {
        File seg_file = "${entity_id}.seg"
    }
}


# Create plots of coverage data and copy-ratio estimates
task PlotSegmentedCopyRatio {
    String entity_id
    String gatk_jar
    File tn_file
    File pre_tn_file
    File segments_file
    File ref_fasta_dict
    String output_dir
    Int mem

    command {
        mkdir -p ${output_dir} && \
        java -Xmx${mem}g -jar ${gatk_jar} PlotSegmentedCopyRatio --tangentNormalized ${tn_file} \
         --preTangentNormalized ${pre_tn_file} --segments ${segments_file} \
         -SD ${ref_fasta_dict} \
         --output ${output_dir} --outputPrefix ${entity_id}
    }

    output {
        File segments_plot="${output_dir}/${entity_id}_FullGenome.png"
        File before_after_normalization_plot="${output_dir}/${entity_id}_Before_After.png"
        File before_after_cr_lim_4="${output_dir}/${entity_id}_Before_After_CR_Lim_4.png"
    }
}
