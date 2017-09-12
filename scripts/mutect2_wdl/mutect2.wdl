#  Run Mutect 2 on a single tumor-normal pair or on a single tumor sample.
#
#  Description of inputs
#  gatk4_jar: java jar file containing gatk 4 (protected)
#  intervals: genomic intervals
#  ref_fasta, ref_fasta_index, ref_dict: reference genome, index, and dictionary
#  tumor_bam, tumor_bam_index: self-explanatory
#  normal_bam, normal_bam_index: self-explanatory
#  pon, pon_index: optional panel of normals and index in vcf format containing known false positves
#  scatter_count: number of parallel jobs when scattering over intervals
#  gnomad, gnomad_index: optional database of known germline variants, obtainable from http://gnomad.broadinstitute.org/downloads
#  variants_for_contamination, variants_for_contamination_index: vcf of common variants with allele frequencies fo calculating contamination
#  is_run_orientation_bias_filter: if true, run the orientation bias filter post-processing step
#  is_run_oncotator: if true, annotate the M2 VCFs using oncotator (to produce a TCGA MAF).  Important:  This requires a docker image and should
#   not be run in environments where docker is unavailable (e.g. SGE cluster on a Broad on-prem VM).  Access to docker
#   hub is also required, since the task will download a public docker image.
#
#
# This WDL needs to decide whether to use the ``gatk_jar`` or ``gatk_jar_override`` for the jar location.  As of cromwell-0.24,
#   this logic *must* go into each task.  Therefore, there is a lot of duplicated code.  This allows users to specify a jar file
#   independent of what is in the docker file.  See the README.md for more info.
#


workflow Mutect2 {
  # gatk4_jar needs to be a String input to the workflow in order to work in a Docker image
  String gatk4_jar
  File? intervals
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File tumor_bam
  File tumor_bam_index
  File? normal_bam
  File? normal_bam_index
  File? pon
  File? pon_index
  Int scatter_count
  File? gnomad
  File? gnomad_index
  File? variants_for_contamination
  File? variants_for_contamination_index
  Boolean is_run_orientation_bias_filter
  Boolean is_run_oncotator

  File? gatk4_jar_override
  File? onco_ds_tar_gz
  String? onco_ds_local_db_dir
  Array[String] artifact_modes
  File picard_jar
  String? m2_extra_args
  String? m2_extra_filtering_args
  String? sequencing_center
  String? sequence_source
  File? default_config_file

  Int? preemptible_attempts
  String gatk_docker
  String basic_bash_docker = "ubuntu:16.04"
  String oncotator_docker

  call SplitIntervals {
    input:
      gatk4_jar = gatk4_jar,
      scatter_count = scatter_count,
      intervals = intervals,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      gatk4_jar_override = gatk4_jar_override,
      preemptible_attempts = preemptible_attempts,
      gatk_docker = gatk_docker
  }

  scatter (subintervals in SplitIntervals.interval_files ) {
    call M2 {
      input: 
        gatk4_jar = gatk4_jar,
        intervals = subintervals,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        tumor_bam = tumor_bam,
        tumor_bam_index = tumor_bam_index,
        normal_bam = normal_bam,
        normal_bam_index = normal_bam_index,
        pon = pon,
        pon_index = pon_index,
        gnomad = gnomad,
        gnomad_index = gnomad_index,
        gatk4_jar_override = gatk4_jar_override,
        preemptible_attempts = preemptible_attempts,
        gatk_docker = gatk_docker,
        m2_extra_args = m2_extra_args
    }
  }

  call MergeVCFs {
    input:
      gatk4_jar = gatk4_jar,
      input_vcfs = M2.output_vcf,
      gatk4_jar_override = gatk4_jar_override,
      preemptible_attempts = preemptible_attempts,
      gatk_docker = gatk_docker
  }

  if (is_run_orientation_bias_filter) {
      call CollectSequencingArtifactMetrics {
        input:
            preemptible_attempts = preemptible_attempts,
            gatk_docker = gatk_docker,
            tumor_bam = tumor_bam,
            tumor_bam_index = tumor_bam_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            picard_jar = picard_jar
      }
  }

  call Filter {
    input:
      gatk4_jar = gatk4_jar,
      gatk4_jar_override = gatk4_jar_override,
      unfiltered_vcf = MergeVCFs.output_vcf,
      intervals = intervals,
      gatk_docker = gatk_docker,
      preemptible_attempts = preemptible_attempts,
      pre_adapter_metrics = CollectSequencingArtifactMetrics.pre_adapter_metrics,
      tumor_bam = tumor_bam,
      tumor_bam_index = tumor_bam_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      artifact_modes = artifact_modes,
      variants_for_contamination = variants_for_contamination,
      variants_for_contamination_index = variants_for_contamination_index,
      m2_extra_filtering_args = m2_extra_filtering_args
  }


  if (is_run_oncotator) {
        call oncotate_m2 {
            input:
                m2_vcf = Filter.filtered_vcf,
                preemptible_attempts = preemptible_attempts,
                oncotator_docker = oncotator_docker,
                onco_ds_tar_gz = onco_ds_tar_gz,
                onco_ds_local_db_dir = onco_ds_local_db_dir,
                sequencing_center = sequencing_center,
                sequence_source = sequence_source,
                default_config_file = default_config_file
        }
  }

  output {
        File unfiltered_vcf = MergeVCFs.output_vcf
        File unfiltered_vcf_index = MergeVCFs.output_vcf_index
        File filtered_vcf = Filter.filtered_vcf
        File filtered_vcf_index = Filter.filtered_vcf_index
        File contamination_table = Filter.contamination_table

        # select_first() fails if nothing resolves to non-null, so putting in "null" for now.
        File? oncotated_m2_maf = select_first([oncotate_m2.oncotated_m2_maf, "null"])
        File? preadapter_detail_metrics = select_first([CollectSequencingArtifactMetrics.pre_adapter_metrics, "null"])
  }
}

task M2 {
  String gatk4_jar
  File? intervals
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File tumor_bam
  File tumor_bam_index
  File? normal_bam
  File? normal_bam_index
  File? pon
  File? pon_index
  File? gnomad
  File? gnomad_index
  File? gatk4_jar_override
  String? m2_extra_args

  # Runtime parameters
  Int? mem
  String gatk_docker
  Int? preemptible_attempts
  Int? disk_space_gb

  command <<<

  # Use GATK Jar override if specified
  GATK_JAR=${gatk4_jar}
  if [[ "${gatk4_jar_override}" == *.jar ]]; then
      GATK_JAR=${gatk4_jar_override}
  fi

  java -Xmx4g -jar $GATK_JAR GetSampleName -I ${tumor_bam} -O tumor_name.txt
  output_vcf=`cat tumor_name.txt`-tumor-only.vcf

  if [[ "_${normal_bam}" == *.bam ]]; then
      java -Xmx4g -jar $GATK_JAR GetSampleName -I ${normal_bam} -O normal_name.txt
      normal_command_line="-I ${normal_bam} -normal `cat normal_name.txt`"
      output_vcf=`cat tumor_name.txt`-vs-`cat normal_name.txt`.vcf
  fi



  java -Xmx4g -jar $GATK_JAR Mutect2 \
    -R ${ref_fasta} \
    -I ${tumor_bam} \
    -tumor `cat tumor_name.txt` \
    $normal_command_line \
    ${"--germline-resource " + gnomad} \
    ${"-pon " + pon} \
    ${"-L " + intervals} \
    -O $output_vcf \
    ${m2_extra_args}
  >>>

  runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
  }

  output {
    File output_vcf = select_first(glob("*.vcf"))
  }
}

task MergeVCFs {
  String gatk4_jar
  Array[File] input_vcfs

  #the scattered M2 task yields identical vcf names that have nothing to do with the scatter index, so we just pick the first
  String output_vcf_name = select_first(input_vcfs)
  File? gatk4_jar_override

  # Runtime parameters
  Int? mem
  String gatk_docker
  Int? preemptible_attempts
  Int? disk_space_gb


  # using MergeVcfs instead of GatherVcfs so we can create indices
  # WARNING 2015-10-28 15:01:48 GatherVcfs  Index creation not currently supported when gathering block compressed VCFs.
  command {
    # Use GATK Jar override if specified
    GATK_JAR=${gatk4_jar}
    if [[ "${gatk4_jar_override}" == *.jar ]]; then
        GATK_JAR=${gatk4_jar_override}
    fi

    java -Xmx2g -jar $GATK_JAR MergeVcfs -I ${sep=' -I ' input_vcfs} -O ${output_vcf_name}
  }

  runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
  }

  output {
    File output_vcf = "${output_vcf_name}"
    File output_vcf_index = "${output_vcf_name}.idx"
  }
}

task CollectSequencingArtifactMetrics {
  File tumor_bam
  File tumor_bam_index
  File ref_fasta
  File ref_fasta_index
  File picard_jar

  # Runtime parameters
  Int? mem
  String gatk_docker
  Int? preemptible_attempts
  Int? disk_space_gb

  command {
        set -e
        java -Xmx4G -jar ${picard_jar} CollectSequencingArtifactMetrics I=${tumor_bam} O="metrics" R=${ref_fasta}

        # Convert to GATK format
        sed -r "s/picard\.analysis\.artifacts\.SequencingArtifactMetrics\\\$PreAdapterDetailMetrics/org\.broadinstitute\.hellbender\.tools\.picard\.analysis\.artifacts\.SequencingArtifactMetrics\$PreAdapterDetailMetrics/g" \
            "metrics.pre_adapter_detail_metrics" > "gatk.pre_adapter_detail_metrics"
  }

  runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
  }

  output {
    File pre_adapter_metrics = "gatk.pre_adapter_detail_metrics"
  }
}

task Filter {
  String gatk4_jar
  File? gatk4_jar_override
  File unfiltered_vcf
  String filtered_vcf = sub(unfiltered_vcf, "\\.vcf$", "-filtered.vcf")
  File? intervals
  File? pre_adapter_metrics
  File? tumor_bam
  File? tumor_bam_index
  File? ref_fasta
  File? ref_fasta_index
  Array[String]? artifact_modes
  File? variants_for_contamination
  File? variants_for_contamination_index
  String? m2_extra_filtering_args

  # Runtime parameters
  Int? mem
  String gatk_docker
  Int? preemptible_attempts
  Int? disk_space_gb

  command {
    set -e

    # Use GATK Jar override if specified
    GATK_JAR=${gatk4_jar}
    if [[ "${gatk4_jar_override}" == *.jar ]]; then
        GATK_JAR=${gatk4_jar_override}
    fi

    touch contamination.table
    if [[ "${variants_for_contamination}" == *.vcf ]]; then
        java -Xmx4g -jar $GATK_JAR GetPileupSummaries -I ${tumor_bam} ${"-L " + intervals} -V ${variants_for_contamination} -O pileups.table
        java -Xmx4g -jar $GATK_JAR CalculateContamination -I pileups.table -O contamination.table
        contamination_cmd="--contamination-table contamination.table"
    fi

    java -Xmx4g -jar $GATK_JAR FilterMutectCalls -V ${unfiltered_vcf} \
      	    -O filtered.vcf $contamination_cmd \
      	    ${m2_extra_filtering_args}

    # FilterByOrientationBias must come after all of the other filtering.
    if [[ ! -z "${pre_adapter_metrics}" ]]; then
        java -Xmx4g -jar $GATK_JAR FilterByOrientationBias -A ${sep=" -A " artifact_modes} \
            -V filtered.vcf -P ${pre_adapter_metrics} --output ${filtered_vcf}
    else
        mv filtered.vcf ${filtered_vcf}
        mv filtered.vcf.idx "${filtered_vcf}.idx"
    fi
  }

  runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
  }

  output {
    File filtered_vcf = "${filtered_vcf}"
    File filtered_vcf_index = "${filtered_vcf}.idx"
    File contamination_table = "contamination.table"
  }
}

task SplitIntervals {
  String gatk4_jar
  Int scatter_count
  File? intervals
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File? gatk4_jar_override

  # Runtime parameters
  Int? mem
  String gatk_docker
  Int? preemptible_attempts
  Int? disk_space_gb

  command {
    # fail if *any* command below (not just the last) doesn't return 0, in particular if GATK SplitIntervals fails
    set -e

    # Use GATK Jar override if specified
    GATK_JAR=${gatk4_jar}
    if [[ "${gatk4_jar_override}" == *.jar ]]; then
        GATK_JAR=${gatk4_jar_override}
    fi

    mkdir interval-files
    java -jar $GATK_JAR SplitIntervals -R ${ref_fasta} ${"-L " + intervals} -scatter ${scatter_count} -O interval-files
    cp interval-files/*.intervals .
  }

  runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 3]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
  }

  output {
    Array[File] interval_files = glob("*.intervals")
  }
}



task oncotate_m2 {
    File m2_vcf
    String oncotator_docker
    File? onco_ds_tar_gz
    String? onco_ds_local_db_dir
    String? oncotator_exe
    String? sequencing_center
    String? sequence_source
    File? default_config_file

  # Runtime parameters
  Int? mem
  String oncotator_docker
  Int? preemptible_attempts
  Int? disk_space_gb

    command <<<

          # fail if *any* command below (not just the last) doesn't return 0, in particular if wget fails
          set -e

          #get tumor and normal sample names from the M2 vcf header
          #normal sample gives an empty string if not present
          case_id=grep 'tumor_sample=' $vcf | sed 's/\#\#tumor_sample\=//'
          control_id=grep 'normal_sample=' $vcf | sed 's/\#\#normal_sample\=//'

          # local db dir is a directory and has been specified
          if [[ -d "${onco_ds_local_db_dir}" ]]; then
              echo "Using local db-dir: ${onco_ds_local_db_dir}"
              echo "THIS ONLY WORKS WITHOUT DOCKER!"
              ln -s ${onco_ds_local_db_dir} onco_dbdir

          elif [[ "${onco_ds_tar_gz}" == *.tar.gz ]]; then
              echo "Using given tar file: ${onco_ds_tar_gz}"
              mkdir onco_dbdir
              tar zxvf ${onco_ds_tar_gz} -C onco_dbdir --strip-components 1

          else
              echo "Downloading and installing oncotator datasources from Broad FTP site..."
              # Download and untar the db-dir
              wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/oncotator/oncotator_v1_ds_April052016.tar.gz
              tar zxvf oncotator_v1_ds_April052016.tar.gz
              ln -s oncotator_v1_ds_April052016 onco_dbdir
          fi

        ${default="/root/oncotator_venv/bin/oncotator" oncotator_exe} --db-dir onco_dbdir/ -c $HOME/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt  \
            -v ${m2_vcf} $case_id.maf.annotated hg19 -i VCF -o TCGAMAF --skip-no-alt --infer-onps --collapse-number-annotations --log_name oncotator.log \
            -a Center:${default="Unknown" sequencing_center} \
            -a source:${default="Unknown" sequence_source} \
            -a normal_barcode:$control_id \
            -a tumor_barcode:$case_id \
            ${"--default_config " + default_config_file}
    >>>

    runtime {
        docker: "${oncotator_docker}"
        memory: select_first([mem, 3]) + " GB"
        bootDiskSizeGb: 12
        disks: "local-disk " + select_first([disk_space_gb, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File oncotated_m2_maf=select_first(glob("*..maf.annotated"))
    }
}
