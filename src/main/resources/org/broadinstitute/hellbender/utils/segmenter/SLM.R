# Useful for debugging
#options(error = quote({dump.frames(dumpto = "CBS_dump", to.file = TRUE); q()}))

# Library used for segmentation
library(SLMSeg)
library(naturalsort)

library(optparse)
option_list <- list(
    make_option(c("--sample_name", "-sample_name"), dest="sample_name", action="store"),
    make_option(c("--denoised_copy_ratio_file", "-denoised_copy_ratio_file"), dest="denoised_copy_ratio_file", action="store"),
    make_option(c("--output_file", "-output_file"), dest="output_file", action="store"),
    make_option(c("--log2_input", "-log2"), dest="log2_input", action="store"),
    make_option(c("--omega" , "-omega" ), dest="omega" , action="store", type="double"),
    make_option(c("--eta" , "-eta" ), dest="eta" , action="store", type="double"),
    make_option(c("--step_eta" , "-step_eta" ), dest="step_eta" , action="store", type="double"),
	make_option(c("--min_length", "-min_len"), dest="min_length", action="store", type="double")
)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

sample_name = opt[["sample_name"]]
denoised_copy_ratio_file = opt[["denoised_copy_ratio_file"]]
output_file = opt[["output_file"]]
log2_input = as.logical(opt[["log2_input"]])
omega = opt[["omega" ]]
eta = opt[["eta" ]]
step_eta = opt[["step_eta" ]]
min_length = opt[["min_length"]]

# Use a function for debugging purposes
segment_data = function(sample_name, denoised_copy_ratio_file, output_file, log2_input, omega, eta, step_eta, min_length) {
	# Read in file and extract needed data
	cr_file = read.table(denoised_copy_ratio_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
	contig = cr_file[, "contig"]
	pos = cr_file[, "stop"]
	# Ability to switch between copy-ratio and log2 copy-ratio
	if (log2_input) {
	    dat = cr_file[, sample_name]
	} else {
	    dat = log2(cr_file[, sample_name])
	}

	# Perform segmentation
	set.seed(25)

    segmented = segment(smooth.CNA(cna_dat), min.width=min_width, weights=weights,alpha=alpha, nperm=nperm, p.method=pmethod, kmax=kmax, nmin=nmin, eta=eta, trim=trim, undo.splits=undosplits, undo.prune=undoprune, undo.SD=undoSD)$output

	# Ensure that there are no too-small values which will be problematic for downstream tools.
	segmented[,"seg.mean"] = 2^segmented[,"seg.mean"]
	segmented[segmented[,"seg.mean"]<.Machine$double.eps,"seg.mean"] = .Machine$double.eps

	# Convention for column names
	colnames(segmented) = c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")

	# Undo conversion of sample_name to R format (retain dashes, spaces, etc.)
	segmented[,"Sample"] = sample_name

	# Change segment start to start of first target in segment (rather than end of first target in segment)
	segmented$rowID = 1:nrow(segmented)
	first_targets = merge(segmented[,c("rowID","Chromosome","Start")], tn[,c("contig","start","stop")], by.x=c("Chromosome","Start"), by.y=c("contig","stop"))
	# merge does not maintain order; restore order of original segmented dataframe
	first_targets = first_targets[order(first_targets$rowID), ]
	segmented[,"Start"] = first_targets[,"start"]
	segmented$rowID = NULL

    # Order based on contig (already ordered based on start)
    sorting = unique(naturalsort(segmented[,"Chromosome"]))
    segmented$Chromosome=factor(segmented$Chromosome, levels=sorting)
    segmented = segmented[order(segmented[,"Chromosome"]),]

	# Output seg file
	print(paste("Writing segment file: ", output_file, sep=""))
	write.table(segmented, file=output_file, sep="\t", quote=FALSE, row.names=FALSE)
}

segment_data(sample_name, tn_file, output_file, log_input, weights_file, min_width, alpha, nperm, pmethod, kmax, nmin, eta, trim, undosplits, undoprune, undoSD)
