# Tasks common to both the CNV somatic panel and case workflows.
#
#############

# Pad targets in the target file by the specified amount (this was found to improve sensitivity and specificity)
task PadTargets {
    File targets
    Int? padding = 250
    String gatk_jar

    # Runtime parameters
    Int? mem = 1
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    # Determine output filename
    String filename = select_first([targets, ""])
    String base_filename = basename(filename, ".tsv")

    command <<<
        echo ${filename}; \
        java -Xmx${mem}g -jar ${gatk_jar} PadTargets \
            --targets ${targets} \
            --padding ${padding} \
            --output ${base_filename}.padded.tsv
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 2]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File padded_targets = "${base_filename}.padded.tsv"
    }
}

# Collect read counts
task CollectReadCounts {
    File? padded_targets
    File bam
    File bam_idx
    Boolean? keep_non_autosomes = false
    Boolean? disable_all_read_filters = false
    Boolean? disable_sequence_dictionary_validation = true
    Boolean? keep_duplicate_reads = true
    Int? wgs_bin_size = 10000
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    String gatk_jar

    # Runtime parameters
    Int? mem = 4
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    # If no padded target file is input, then do WGS workflow
    Boolean is_wgs = !defined(padded_targets)

    # Sample name is derived from the bam filename
    String base_filename = basename(bam, ".bam")

    String read_counts_tsv_filename = "${base_filename}.readCounts.tsv"
    String read_counts_hdf5_filename = if is_wgs then "${base_filename}.readCounts.hdf5" else ""

    command <<<
        if [ ${is_wgs} = true ]
            then
                java -Xmx${mem}g -jar ${gatk_jar} SparkGenomeReadCounts \
                    --input ${bam} \
                    --reference ${ref_fasta} \
                    --binsize ${wgs_bin_size} \
                    --keepXYMT ${keep_non_autosomes} \
                    --disableToolDefaultReadFilters ${disable_all_read_filters} \
                    --disableSequenceDictionaryValidation ${disable_sequence_dictionary_validation} \
                    $(if [ ${keep_duplicate_reads} = true ]; then echo " --disableReadFilter NotDuplicateReadFilter "; else echo ""; fi) \
                    --output ${read_counts_tsv_filename} \
                    --writeHdf5
            else
                java -Xmx${mem}g -jar ${gatk_jar} CalculateTargetCoverage \
                    --input ${bam} \
                    --reference ${ref_fasta} \
                    --targets ${padded_targets} \
                    --groupBy SAMPLE \
                    --transform RAW \
                    --targetInformationColumns FULL \
                    --interval_set_rule UNION \
                    --interval_merging_rule OVERLAPPING_ONLY \
                    --interval_padding 0 \
                    --secondsBetweenProgressUpdates 10.0 \
                    --disableToolDefaultReadFilters ${disable_all_read_filters} \
                    --disableSequenceDictionaryValidation ${disable_sequence_dictionary_validation} \
                    $(if [ ${keep_duplicate_reads} = true ]; then echo " --disableReadFilter NotDuplicateReadFilter "; else echo ""; fi) \
                    --output ${read_counts_tsv_filename}
        fi
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(bam, "GB")) + 50]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        String entity_id = base_filename
        File read_counts = read_counts_tsv_filename
        File read_counts_hdf5 = read_counts_hdf5_filename   #"" if is_wgs = false
    }
}

# Create a target file with GC annotations
task AnnotateTargets {
    String entity_id
    File targets
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    String gatk_jar

    # Runtime parameters
    Int? mem = 4
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    command {
        java -Xmx${mem}g -jar ${gatk_jar} AnnotateTargets \
            --targets ${targets} \
            --reference ${ref_fasta} \
            --interval_merging_rule OVERLAPPING_ONLY \
            --output ${entity_id}.annotated.tsv
    }

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(ref_fasta, "GB")) + 50]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File annotated_targets = "${entity_id}.annotated.tsv"
    }
}