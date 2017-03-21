# Given tumor-normal pairs compare Mutect calls in tumor-only mode to those of tumor-normal mode
import "mutect2.wdl" as m2


task Concordance {
  File gatk4_jar
  File? gatk4_jar_override
  File? intervals
  File truth_vcf
  File truth_vcf_idx
  File eval_vcf
  File eval_vcf_idx

  command {

        GATK_JAR=${gatk4_jar}
        if [[ "${gatk4_jar_override}" == *.jar ]]; then
            GATK_JAR=${gatk4_jar_override}
        fi

      java -jar $GATK_JAR Concordance ${"-L " + intervals} \
        -truth ${truth_vcf} -eval ${eval_vcf} -tpfn "true_positives_and_false_negatives.vcf" \
        -tpfp "true_positives_and_false_positives.vcf" \
        -summary summary.tsv
  }

    runtime {
        memory: "5 GB"
    }

  output {
        File tpfn = "true_positives_and_false_negatives.vcf"
        File tpfn_idx = "true_positives_and_false_negatives.vcf.idx"
        File tpfp = "true_positives_and_false_positives.vcf"
        File tpfp_idx = "true_positives_and_false_positives.vcf.idx"
        File summary = "summary.tsv"
  }
}

workflow Mutect2EvaluateTumorOnly {
	File gatk4_jar
	Int scatter_count
	# pair_list file is a tsv file with the following six columns in this order.
	# tumor_bam, tumor_bam_index, tumor_sample_name, normal_bam, normal_bam_index, normal_sample_name
	File pair_list
	Array[Array[String]] pairs = read_tsv(pair_list)
	File intervals
	File ref_fasta
	File ref_fasta_index
	File ref_dict
	File? pon
	File? pon_index
	File? dbsnp
	File? dbsnp_index
	File? cosmic
	File? cosmic_index
	Boolean is_run_orientation_bias_filter
	Boolean is_run_oncotator
	String oncotator_docker
	String m2_docker
	File? gatk4_jar_override
	Int preemptible_attempts
	Array[String] artifact_modes

	scatter(pair in pairs) {
		call m2.Mutect2 as TumorNormal {
			input:
				gatk4_jar=gatk4_jar,
				intervals=intervals,
				ref_fasta=ref_fasta,
				ref_fasta_index=ref_fasta_index,
				ref_dict=ref_dict,
				tumor_bam=pair[0],
				tumor_bam_index=pair[1],
				tumor_sample_name=pair[2],
				normal_bam=pair[3],
				normal_bam_index=pair[4],
				normal_sample_name=pair[5],
				pon=pon,
				pon_index=pon_index,
				scatter_count=scatter_count,
				dbsnp=dbsnp,
				dbsnp_index=dbsnp_index,
				cosmic=cosmic,
				cosmic_index=cosmic_index,
				picard_jar = picard_jar,
                is_run_orientation_bias_filter = is_run_orientation_bias_filter,
                is_run_oncotator=is_run_oncotator,
                oncotator_docker=oncotator_docker,
                m2_docker = m2_docker,
                gatk4_jar_override = gatk4_jar_override,
                preemptible_attempts = preemptible_attempts,
                artifact_modes = artifact_modes
		}

		call m2.Mutect2 as TumorOnly {
        			input:
        				gatk4_jar=gatk4_jar,
        				intervals=intervals,
        				ref_fasta=ref_fasta,
        				ref_fasta_index=ref_fasta_index,
        				ref_dict=ref_dict,
        				tumor_bam=pair[0],
        				tumor_bam_index=pair[1],
        				tumor_sample_name=pair[2],
        				pon=pon,
        				pon_index=pon_index,
        				scatter_count=scatter_count,
        				dbsnp=dbsnp,
        				dbsnp_index=dbsnp_index,
        				cosmic=cosmic,
        				cosmic_index=cosmic_index,
        				picard_jar = picard_jar,
                        is_run_orientation_bias_filter = is_run_orientation_bias_filter,
                        is_run_oncotator=is_run_oncotator,
                        oncotator_docker=oncotator_docker,
                        m2_docker = m2_docker,
                        gatk4_jar_override = gatk4_jar_override,
                        preemptible_attempts = preemptible_attempts,
                        artifact_modes = artifact_modes
        		}

		call Concordance {
		    input:
		        gatk4_jar = gatk4_jar,
                gatk4_jar_override = gatk4_jar_override,
                intervals = intervals,
                truth_vcf = TumorNormal.filtered_vcf, #note, no orientation bias since it's optional output
                truth_vcf_idx = TumorNormal.filtered_vcf_index,
                eval_vcf = TumorNormal.filtered_vcf,
                eval_vcf_idx = TumorNormal.filtered_vcf_index,
		}
	}

	output {
		Array[File] tpfn = Concordance.tpfn
        Array[File] tpfn_idx = Concordance.tpfn_idx
        Array[File] tpfp = Concordance.tpfp
        Array[File] tpf_idx = Concordance.tpfp_idx
        Array[File] summary = Concordance.summary
	}
}
