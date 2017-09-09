# Workflow for creating a GATK CNV Panel of Normals given a list of normal samples. Supports both WGS and WES.
#
# Notes:
#
# - Input file (normal_bams_list) must contain file paths to bam and bam index files separated by tabs in the following format:
#    normal_bam_1    bam_idx_1
#    normal_bam_2    bam_idx_2
#    ...
#
# - The target file (targets) is required for the WES workflow and should be a TSV file with the column headers:
#    contig    start    stop    name
#   These targets will be padded on both sides by the amount specified by PadTargets.padding (default 250).
#
# - If a target file is not provided, then the WGS workflow will be run instead and the specified value of
#   wgs_bin_size (default 10000) will be used.
#
# - Example invocation:
#    java -jar cromwell.jar run cnv_somatic_panel_workflow.wdl myParameters.json
#   See cnv_somatic_panel_workflow_template.json for a template json file to modify with your own parameters (please save
#   your modified version with a different filename and do not commit to the gatk repository).
#
#############

import "cnv_common_tasks.wdl" as CNVTasks

workflow CNVSomaticPanelWorkflow {
    # Workflow input files
    File? targets
    File normal_bams_list
    Array[Array[String]]+ normal_bams = read_tsv(normal_bams_list)
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    File gatk_jar
    String gatk_docker

    # PoN name
    String pon_entity_id

    # If true, AnnotateTargets will be run to create GC annotations and
    # explicit GC correction will be performed by CreateReadCountPanelOfNormals before PCA is performed
    Boolean do_explicit_gc_correction = false

    # If no target file is input, then do WGS workflow
    Boolean is_wgs = !defined(targets)

    if (!is_wgs) {
        call CNVTasks.PadTargets {
            input:
                # The task will fail if targets is not defined when it gets here, but that should not be allowed to happen.
                targets = select_first([targets, ""]),
                gatk_jar = gatk_jar,
                gatk_docker = gatk_docker
        }
    }

    scatter (normal_bam in normal_bams) {
        call CNVTasks.CollectReadCounts {
            input:
                padded_targets = PadTargets.padded_targets,
                bam = normal_bam[0],
                bam_idx = normal_bam[1],
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta_dict = ref_fasta_dict,
                gatk_jar = gatk_jar,
                gatk_docker = gatk_docker
        }
    }

	if (do_explicit_gc_correction) {
		call CNVTasks.AnnotateTargets {
            input:
                entity_id = pon_entity_id,
                targets = CollectReadCounts.read_counts[0],
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta_dict = ref_fasta_dict,
                gatk_jar = gatk_jar,
                gatk_docker = gatk_docker
        }
	}

    call CreateReadCountPanelOfNormals {
        input:
            pon_entity_id = pon_entity_id,
            read_count_files = if is_wgs then CollectReadCounts.read_counts_hdf5 else CollectReadCounts.read_counts,
            annotated_intervals = AnnotateTargets.annotated_targets,
            gatk_jar = gatk_jar,
            gatk_docker = gatk_docker
    }

    output {
        File read_count_pon = CreateReadCountPanelOfNormals.read_count_pon
    }
}

# Create read-count panel of normals
task CreateReadCountPanelOfNormals {
    String pon_entity_id
    Array[File] read_count_files
    Float? minimum_interval_median_percentile = 25.
    Float? maximum_zeros_in_sample_percentage = 2.
    Float? maximum_zeros_in_interval_percentage = 5.
    Float? extreme_sample_median_percentile = 2.5
    Float? extreme_outlier_truncation_percentile = 0.1
    Int? number_of_eigensamples = 20
    File? annotated_intervals   #null = do not perform explicit GC correction by default
    String gatk_jar

    # Runtime parameters
    Int? mem = 4
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    command {
        java -Xmx${mem}g -jar ${gatk_jar} CreateReadCountPanelOfNormals \
            --input ${sep=" --input " read_count_files} \
            --minimumIntervalMedianPercentile ${minimum_interval_median_percentile} \
            --maximumZerosInSamplePercentage ${maximum_zeros_in_sample_percentage} \
            --maximumZerosInIntervalPercentage ${maximum_zeros_in_interval_percentage} \
            --extremeSampleMedianPercentile ${extreme_sample_median_percentile} \
            --extremeOutlierTruncationPercentile ${extreme_outlier_truncation_percentile} \
            --numberOfEigensamples ${number_of_eigensamples} \
            --annotatedIntervals ${annotated_intervals} \
            --output ${pon_entity_id}.pon.hdf5
    }

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File read_count_pon = "${pon_entity_id}.pon.hdf5"
    }
}