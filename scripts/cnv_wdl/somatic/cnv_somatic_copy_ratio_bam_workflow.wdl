# Subworkflow for running GATK CNV on a single BAM. Supports both WGS and WES.
#
# Notes:
#
# - The padded target file (padded_targets) is required for the WES workflow and should be a TSV file with the column headers:
#    contig    start    stop    name
#
# - If a target file is not provided, then the WGS workflow will be run instead and the specified value of
#   wgs_bin_size (default 10000) will be used.
#
#############

import "cnv_common_tasks.wdl" as CNVTasks

workflow CNVSomaticCopyRatioBAMWorkflow {
    # Workflow input files
    File? padded_targets
    File bam
    File bam_idx
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    File read_count_pon
    Boolean do_gc_correction
    String gatk_jar
    String gatk_docker

    call CNVTasks.CollectReadCounts {
        input:
            padded_targets = padded_targets,
            bam = bam,
            bam_idx = bam_idx,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            gatk_jar = gatk_jar,
            gatk_docker = gatk_docker
    }

    call DenoiseReadCounts {
        input:
            entity_id = CollectReadCounts.entity_id,
            read_counts = CollectReadCounts.read_counts,
            read_count_panel_of_normals = read_count_panel_of_normals,
            gatk_jar = gatk_jar,
            gatk_docker = gatk_docker
    }

    call ModelSegments {
        input:
            entity_id = CollectReadCounts.entity_id,
            denoised_copy_ratios = DenoiseReadCounts.denoised_copy_ratios,
            gatk_jar = gatk_jar,
            gatk_docker = gatk_docker
    }

    call CallSegments {
        input:
            entity_id = CollectReadCounts.entity_id,
            denoised_copy_ratios = DenoiseReadCounts.denoised_copy_ratios,
            copy_ratio_segments = ModelSegments.copy_ratio_segments,
            gatk_jar = gatk_jar,
            gatk_docker = gatk_docker
    }

    call PlotSegmentedCopyRatio  {
        input:
            entity_id = CollectReadCounts.entity_id,
            standardized_copy_ratios = DenoiseReadCounts.standardized_copy_ratios,
            denoised_copy_ratios = DenoiseReadCounts.denoised_copy_ratios,
            called_copy_ratio_segments = CallSegments.called_copy_ratio_segments,
            ref_fasta_dict = ref_fasta_dict,
            gatk_jar = gatk_jar,
            gatk_docker = gatk_docker
    }

    output {
        String entity_id = CollectReadCounts.entity_id
        File read_counts = CollectReadCounts.read_counts
        File standardized_copy_ratio = DenoiseReadCounts.standardized_copy_ratio
        File denoised_copy_ratio = DenoiseReadCounts.denoised_copy_ratio
        File called_copy_ratio_segments = CallSegments.called_copy_ratio_segments
        File segmented_copy_ratio_plot = PlotSegmentedCopyRatio.segmented_copy_ratio_plot
        File copy_ratio_before_after_denoising_plot = PlotSegmentedCopyRatio.copy_ratio_before_after_denoising_plot
        File copy_ratio_before_after_denoising_lim_4_plot = PlotSegmentedCopyRatio.copy_ratio_before_after_denoising_lim_4_plot
    }
}

# Denoise the coverage
task DenoiseReadCounts {
    String entity_id
    File read_counts
    File read_count_panel_of_normals
    Int? number_of_eigensamples
    String gatk_jar

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    command {
        java -Xmx${default=4 mem}g -jar ${gatk_jar} NormalizeSomaticReadCounts \
            --input ${read_counts} \
            --panelOfNormals ${read_count_panel_of_normals} \
            --numberOfEigensamples ${default="null" number_of_eigensamples} \
            --standardizedCR ${entity_id}.standardizedCR.tsv \
            --denoisedCR ${entity_id}.denoisedCR.tsv
    }

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(read_count_pon, "GB")) + 50]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File standardized_copy_ratios = "${entity_id}.standardizedCR.tsv"
        File denoised_copy_ratios = "${entity_id}.denoisedCR.tsv"
    }
}

# Segment the denoised copy-ratio profile
task ModelSegments {
    String entity_id
    File denoised_copy_ratios
    Int? max_num_segments_per_chromosome
    Float? kernel_variance
    Int? kernel_approximation_dimension
    Array[Int]? window_sizes
    Float? num_changepoints_penalty_linear_factor
    Float? num_changepoints_penalty_log_linear_factor
    String gatk_jar

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    command {
        java -Xmx${default=4 mem}g -jar ${gatk_jar} ModelSegments \
            --input ${denoised_copy_ratios} \
            --maxNumSegmentsPerChromosome ${default="50" max_num_segments_per_chromosome} \
            --kernelVariance ${default="0." kernel_variance} \
            --kernelApproximationDimension ${default="100" kernel_approximation_dimension} \
            --windowSizes ${default="[8, 16, 32, 64, 128, 256]" window_sizes} \
            --numChangepointsPenaltyLinearFactor ${default="1." num_changepoints_penalty_linear_factor} \
            --numChangepointsPenaltyLogLinearFactor ${default="1." num_changepoints_penalty_log_linear_factor} \
            --output ${entity_id}.seg
    }

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File copy_ratio_segments = "${entity_id}.seg"
    }
}

# Make calls (amplified, neutral, or deleted) on each segment
task CallSegments {
    String entity_id
    File denoised_copy_ratios
    File copy_ratio_segments
    String gatk_jar

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    command {
        java -Xmx${default=4 mem}g -jar ${gatk_jar} CallSegments \
            --tangentNormalized ${denoised_copy_ratios} \
            --segments ${copy_ratio_segments} \
            --legacy false \
            --output ${entity_id}.called
    }

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File called_copy_ratio_segments = "${entity_id}.called"
    }
}

# Create plots of coverage data and copy-ratio estimates
task PlotSegmentedCopyRatio {
    String entity_id
    File standardized_copy_ratios
    File denoised_copy_ratios
    File called_copy_ratio_segments
    File ref_fasta_dict
    String? output_dir
    String gatk_jar

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    # If optional output_dir not specified, use "."
    String output_dir_ = select_first([output_dir, "."])

    command {
        mkdir -p ${output_dir_}; \
        java -Xmx${default=4 mem}g -jar ${gatk_jar} PlotSegmentedCopyRatio \
            --preTangentNormalized ${standardized_copy_ratios} \
            --tangentNormalized ${denoised_copy_ratios} \
            --segments ${called_copy_ratio_segments} \
            -SD ${ref_fasta_dict} \
            --output ${output_dir_} \
            --outputPrefix ${entity_id}
    }

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File segmented_copy_ratio_plot = "${output_dir_}/${entity_id}_FullGenome.png"
        File copy_ratio_before_after_denoising_plot = "${output_dir_}/${entity_id}_Before_After.png"
        File copy_ratio_before_after_denoising_lim_4_plot = "${output_dir_}/${entity_id}_Before_After_CR_Lim_4.png"
    }
}