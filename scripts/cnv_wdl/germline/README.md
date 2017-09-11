## Running the Germline CNV WDL

### Which WDL should you use?
- Building a panel of normals (PoN): ``cnv_germline_panel_creation_workflow.wdl``
- Calling events on a single normal sample: ``cnv_germline_single_sample_calling_workflow.wdl``
- Calling events on a cohort of normal samples: ``cnv_germline_cohort_calling_workflow.wdl``

#### Setting up parameter json file for a run

To get started, copy the relevant ``*_template.json`` for the workflow you wish to run and adjust parameters accordingly.
You can find all required resource inputs needed to run the workflows in the ``/resources`` directory. These inputs could be run out-of-the-box.

*Please note that there are task-level parameters that do not appear in the template files.  These are set to reasonable values by default, but can also be adjusted if desired.

#### Fields of germline CNV panel of normals creation workflow

  ``CNVGermlinePanelWorkflow.sex_genotypes`` -- path to table of per-sample sex genotypes
  ``CNVGermlinePanelWorkflow.contig_ploidy_annotations`` --  path to the germline contig ploidy annotations table; located in ``/resources`` directory
  ``CNVGermlinePanelWorkflow.transition_prior_table`` -- path to copy number transition priors table; located in ``/resources`` directory
  ``CNVGermlinePanelWorkflow.transition_matrix_XY_Y`` -- path to copy number transition prior for Y contig for XY-genotyped samples; located in ``/resources`` directory
  ``CNVGermlinePanelWorkflow.transition_matrix_XX_X`` -- path to copy number transition prior for X contig for XX-genotyped samples; located in ``/resources`` directory
  ``CNVGermlinePanelWorkflow.transition_matrix_XY_X`` -- path to copy number transition prior for X contig for XY-genotyped samples; located in ``/resources`` directory
  ``CNVGermlinePanelWorkflow.transition_matrix_XX_Y`` -- path to copy number transition prior for Y contig for XX-genotyped samples; located in ``/resources`` directory
  ``CNVGermlinePanelWorkflow.transition_matrix_autosomal`` -- path to transition prior on autosomal loci; located in ``/resources`` directory,
  ``CNVGermlinePanelWorkflow.normal_bams_list`` -- TSV file consisting of corresponding bam and corresponding index files as described in cnv_germline_panel_creation_workflow.wdl
  ``CNVGermlinePanelWorkflow.pon_output_path`` -- name of the final output directory
  ``CNVGermlinePanelWorkflow.num_latents`` -- (advanced) maximum number of principal components. Must be strictly less than the number of samples. The recommended value is 20 ~ 30 for large cohorts. For smaller cohorts, use 0.5 * number of samples. Unnecessary principal components are automatically pruned during PoN creation
  ``CNVGermlinePanelWorkflow.ref_fasta`` -- path to reference fasta file
  ``CNVGermlinePanelWorkflow.ref_fasta_dict`` -- path to reference dict file
  ``CNVGermlinePanelWorkflow.ref_fasta_fai`` -- path to reference fasta fai file
  ``CNVGermlinePanelWorkflow.gatk_jar`` -- absolute path to gatk-protected.jar
  ``CNVGermlinePanelWorkflow.targets`` -- (optional) Target file (NOT in BED format) corresponding to the genomic loci of enriched targets in WES sample (e.g. Agilent, Illumina, etc). Please run ConvertBedToTargetFile to convert a BED file to a target file. If provided, then WES workflow will be run; otherwise, WGS workflow will be run


#### Fields of germline CNV single sample calling workflow

The reference used must be the same between PoN and case samples.

  ``CNVGermlineSingleSampleWorkflow.sex_genotypes`` -- path to table of per-sample sex genotypes
  ``CNVGermlineSingleSampleWorkflow.contig_ploidy_annotations`` --  path to the germline contig ploidy annotations table; located in ``/resources`` directory
  ``CNVGermlineSingleSampleWorkflow.transition_prior_table`` -- path to copy number transition priors table; located in ``/resources`` directory
  ``CNVGermlineSingleSampleWorkflow.transition_matrix_XY_Y`` -- path to copy number transition prior for Y contig for XY-genotyped samples; located in ``/resources`` directory
  ``CNVGermlineSingleSampleWorkflow.transition_matrix_XX_X`` -- path to copy number transition prior for X contig for XX-genotyped samples; located in ``/resources`` directory
  ``CNVGermlineSingleSampleWorkflow.transition_matrix_XY_X`` -- path to copy number transition prior for X contig for XY-genotyped samples; located in ``/resources`` directory
  ``CNVGermlineSingleSampleWorkflow.transition_matrix_XX_Y`` -- path to copy number transition prior for Y contig for XX-genotyped samples; located in ``/resources`` directory
  ``CNVGermlineSingleSampleWorkflow.transition_matrix_autosomal`` -- path to transition prior on autosomal loci; located in ``/resources`` directory,
  ``CNVGermlineSingleSampleWorkflow.output_path`` -- name of the final output directory
  ``CNVGermlineSingleSampleWorkflow.num_latents`` -- (advanced) maximum number of principal components. Must be strictly less than the number of samples. The recommended value is 20 ~ 30 for large cohorts. For smaller cohorts, use 0.5 * number of samples. Unnecessary principal components are automatically pruned during PoN creation
  ``CNVGermlineSingleSampleWorkflow.model_path`` -- absolute path of the PoN model (posterior_finals directory of the panel creation output)
  ``CNVGermlineSingleSampleWorkflow.normal_bam`` -- path to the normal bam file
  ``CNVGermlineSingleSampleWorkflow.normal_bam_idx`` -- path to the corresponding bam index file
  ``CNVGermlineSingleSampleWorkflow.ref_fasta`` -- path to reference fasta file
  ``CNVGermlineSingleSampleWorkflow.ref_fasta_dict`` -- path to reference dict file
  ``CNVGermlineSingleSampleWorkflow.ref_fasta_fai`` -- path to reference fasta fai file
  ``CNVGermlineSingleSampleWorkflow.gatk_jar`` -- absolute path to gatk-protected.jar
  ``CNVGermlineSingleSampleWorkflow.targets`` -- (optional) Target file (NOT in BED format) corresponding to the genomic loci of enriched targets in WES sample (e.g. Agilent, Illumina, etc). Please run ConvertBedToTargetFile to convert a BED file to a target file. If provided, then WES workflow will be run; otherwise, WGS workflow will be run


#### Fields of germline CNV cohort calling workflow

The reference used must be the same between PoN and case samples.

  ``CNVGermlineCohortCallingWorkflow.sex_genotypes`` -- path to table of per-sample sex genotypes
  ``CNVGermlineCohortCallingWorkflow.contig_ploidy_annotations`` --  path to the germline contig ploidy annotations table; located in ``/resources`` directory
  ``CNVGermlineCohortCallingWorkflow.transition_prior_table`` -- path to copy number transition priors table; located in ``/resources`` directory
  ``CNVGermlineCohortCallingWorkflow.transition_matrix_XY_Y`` -- path to copy number transition prior for Y contig for XY-genotyped samples; located in ``/resources`` directory
  ``CNVGermlineCohortCallingWorkflow.transition_matrix_XX_X`` -- path to copy number transition prior for X contig for XX-genotyped samples; located in ``/resources`` directory
  ``CNVGermlineCohortCallingWorkflow.transition_matrix_XY_X`` -- path to copy number transition prior for X contig for XY-genotyped samples; located in ``/resources`` directory
  ``CNVGermlineCohortCallingWorkflow.transition_matrix_XX_Y`` -- path to copy number transition prior for Y contig for XX-genotyped samples; located in ``/resources`` directory
  ``CNVGermlineCohortCallingWorkflow.transition_matrix_autosomal`` -- path to transition prior on autosomal loci; located in ``/resources`` directory
  ``CNVGermlineCohortCallingWorkflow.output_path`` -- name of the final output directory
  ``CNVGermlineCohortCallingWorkflow.num_latents`` -- (advanced) maximum number of principal components. Must be strictly less than the number of samples. The recommended value is 20 ~ 30 for large cohorts. For smaller cohorts, use 0.5 * number of samples. Unnecessary principal components are automatically pruned during PoN creation
  ``CNVGermlineCohortCallingWorkflow.model_path`` -- absolute path of the PoN model (posterior_finals directory of the panel creation output)
  ``CNVGermlineCohortCallingWorkflow.normal_bams_list`` -- TSV file consisting of corresponding bam and corresponding index files as described in cnv_germline_cohort_calling_workflow.wdl
  ``CNVGermlineCohortCallingWorkflow.ref_fasta`` -- path to reference fasta file
  ``CNVGermlineCohortCallingWorkflow.ref_fasta_dict`` -- path to reference dict file
  ``CNVGermlineCohortCallingWorkflow.ref_fasta_fai`` -- path to reference fasta fai file
  ``CNVGermlineCohortCallingWorkflow.gatk_jar`` -- absolute path to gatk-protected.jar
  ``CNVGermlineCohortCallingWorkflow.targets`` -- (optional) Target file (NOT in BED format) corresponding to the genomic loci of enriched targets in WES sample (e.g. Agilent, Illumina, etc). Please run ConvertBedToTargetFile to convert a BED file to a target file. If provided, then WES workflow will be run; otherwise, WGS workflow will be run

In addition, there are several task-level parameters that may be set by advanced users; for example:

- ``CNVGermlineCohortCallingWorkflow.CollectReadCounts.wgs_bin_size`` -- Size of bins (in bp) for WGS coverage collection.  *This must be the same value used for all samples.*  Ignored if not running WGS.
- ``CNVGermlineCohortCallingWorkflow.PadTargets.padding`` -- Amount of padding (in bp) to add to both sides of targets for WES coverage collection.  *This must be the same value used for all samples.*  Ignored if not running WES.



Further explanation of these task-level parameters may be found by invoking the ``--help`` documentation available in the gatk-protected.jar for each tool.