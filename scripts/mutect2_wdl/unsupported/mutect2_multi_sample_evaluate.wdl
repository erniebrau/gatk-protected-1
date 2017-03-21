# Note: input json is identical to mutect2_multi_sample.wdl except for the additioanl columns for the truth vcfs
# in the pairs list file

#  Run Mutect 2 on a list of tumors or tumor-normal pairs and evaluate calls against truth vcfs
#
#  Description of inputs
#  gatk4_jar: java jar file containing gatk 4 (protected)
#  intervals: genomic intervals
#  ref_fasta, ref_fasta_index, ref_dict: reference genome, index, and dictionary
#  pon, pon_index: optional panel of normals and index in vcf format containing known false positves
#  scatter_count: number of parallel jobs when scattering over intervals
#  dbsnp, dbsnp_index: optional database of known germline variants
#  cosmic, cosmic_index: optional database of known somatic variants
#  variants_for_contamination, variants_for_contamination_index: vcf of common variants with allele frequencies fo calculating contamination
#  is_run_orientation_bias_filter: if true, run the orientation bias filter post-processing step
#  pair_list: a tab-separated table with no header in the following format:
#   TUMOR_1_BAM</TAB>TUMOR_1_BAM_INDEX</TAB>TUMOR_1_SAMPLE</TAB>NORMAL_1_BAM</TAB>NORMAL_1_BAM_INDEX</TAB>NORMAL_1_SAMPLE</TAB>TRUTH_VCF_1</TAB>TRUTH_VCF_1_IDX
#   TUMOR_2_BAM</TAB>TUMOR_2_BAM_INDEX</TAB>TUMOR_2_SAMPLE</TAB>NORMAL_2_BAM</TAB>NORMAL_2_BAM_INDEX</TAB>NORMAL_2_SAMPLE</TAB>TRUTH_VCF_2</TAB>TRUTH_VCF_2_IDX
#   . . .
#  Temporarily, while waiting for a Cromwell bug to be resolved, tumor-only input looks like
#  TUMOR_1_BAM</TAB>TUMOR_1_BAM_INDEX</TAB>TUMOR_1_SAMPLE</TAB>NO_NORMAL</TAB>NO_NORMAL</TAB>NO_NORMAL</TAB>TRUTH_VCF_1</TAB>TRUTH_VCF_1_IDX
#  TUMOR_2_BAM</TAB>TUMOR_2_BAM_INDEX</TAB>TUMOR_2_SAMPLE</TAB>NO_NORMAL</TAB>NO_NORMAL</TAB>NO_NORMAL</TAB>TRUTH_VCF_2</TAB>TRUTH_VCF_2_IDX
#   . . .
#  That is, you actually write out "NO NORMAL" thrice per row

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


workflow Mutect2_Multi {
    # gatk4_jar needs to be a String input to the workflow in order to work in a Docker image
	String gatk4_jar
	Int scatter_count
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
	File? variants_for_contamination
    File? variants_for_contamination_index
	Boolean is_run_orientation_bias_filter
	Boolean is_run_oncotator
    String m2_docker
    String oncotator_docker
    File? gatk4_jar_override
    Int preemptible_attempts
    File? onco_ds_tar_gz
    String? onco_ds_local_db_dir
    Array[String] artifact_modes

	scatter( row in pairs ) {
	    #The non-hack way, but there's a bug
	    #      In WDL, variables inside the block can be used outside the block.
	    #      If the conditional block is run, they retain their values.
	    #      Otherwise, they evaluate to null, which in WDL is equivalent to an empty optional
	    #      If we simply tried to use eg row[3] below it could cause an out-of-bounds exception
        #if(length(pairs[n]) == 6) {
        #    File normal_bam = row[3]
        #    File normal_bam_index = row[4]
        #    String normal_sample_name = row[5]
        #}

            call m2.Mutect2 {
                input:
                    gatk4_jar = gatk4_jar,
                    intervals = intervals,
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    ref_dict = ref_dict,
                    tumor_bam = row[0],
                    tumor_bam_index = row[1],
                    tumor_sample_name = row[2],
                    normal_bam = row[3],
                    normal_bam_index = row[4],
                    normal_sample_name = row[5],
                    pon = pon,
                    pon_index = pon_index,
                    scatter_count = scatter_count,
                    dbsnp = dbsnp,
                    dbsnp_index = dbsnp_index,
                    cosmic = cosmic,
                    cosmic_index = cosmic_index,
                    variants_for_contamination = variants_for_contamination,
                    variants_for_contamination_index = variants_for_contamination_index,
                    is_run_orientation_bias_filter = is_run_orientation_bias_filter,
                    is_run_oncotator=is_run_oncotator,
                    oncotator_docker=oncotator_docker,
                    m2_docker = m2_docker,
                    gatk4_jar_override = gatk4_jar_override,
                    preemptible_attempts = preemptible_attempts,
                    onco_ds_tar_gz = onco_ds_tar_gz,
                    onco_ds_local_db_dir = onco_ds_local_db_dir,
                    artifact_modes = artifact_modes
            }

            call Concordance {
            		    input:
            		        gatk4_jar = gatk4_jar,
                            gatk4_jar_override = gatk4_jar_override,
                            intervals = intervals,
                            truth_vcf = row[6]
                            truth_vcf_idx = row[7]
                            eval_vcf = Mutect2.filtered_vcf,
                            eval_vcf_idx = Mutect2.filtered_vcf_index,
            		}
    }

    output {
        Array[File] summary = Concordance.summary
        Array[File] tpfp = Concordance.tpfp
        Array[File] tpfn = Concordance.tpfn
    }
}