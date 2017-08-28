package org.broadinstitute.hellbender.tools.copynumber.legacy;

import com.google.common.collect.ImmutableSet;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation.CopyRatioKernelSegmenter;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.stream.Collectors;

/**
 * TODO
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Model segmented copy ratio from denoised read counts.",
        oneLineSummary = "Model segmented copy ratio from denoised read counts",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class ModelSegments extends CommandLineProgram {
    private static final String MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_LONG_NAME = "maxNumSegmentsPerChromosome";
    private static final String MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_SHORT_NAME = "maxNumSegPerChr";

    private static final String KERNEL_VARIANCE_LONG_NAME = "kernelVariance";
    private static final String KERNEL_VARIANCE_SHORT_NAME = "kernVar";

    private static final String KERNEL_APPROXIMATION_DIMENSION_LONG_NAME = "kernelApproximationDimension";
    private static final String KERNEL_APPROXIMATION_DIMENSION_SHORT_NAME = "kernApproxDim";

    private static final String WINDOW_SIZES_LONG_NAME = "windowSizes";
    private static final String WINDOW_SIZES_SHORT_NAME = "winSizes";

    private static final String NUM_CHANGEPOINTS_PENALTY_LINEAR_FACTOR_LONG_NAME = "numChangepointsPenaltyLinearFactor";
    private static final String NUM_CHANGEPOINTS_PENALTY_LINEAR_FACTOR_SHORT_NAME = "numChangepointsPenLin";

    private static final String NUM_CHANGEPOINTS_PENALTY_LOG_LINEAR_FACTOR_LONG_NAME = "numChangepointsPenaltyLogLinearFactor";
    private static final String NUM_CHANGEPOINTS_PENALTY_LOG_LINEAR_FACTOR_SHORT_NAME = "numChangepointsPenLogLin";

    @Argument(
            doc = "Input file containing denoised copy-ratio profile (output of DenoiseReadCounts).",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME
    )
    private String inputDenoisedCopyRatioProfileFile;

    @Argument(
            doc = "Output file for segmented copy-ratio model result.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputSegmentsFile;

    @Argument(
            doc = "Maximum number of segments allowed per chromosome.",
            fullName = MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_LONG_NAME,
            shortName = MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_SHORT_NAME,
            optional = true,
            minValue = 1
    )
    private int maxNumSegmentsPerChromosome = 100;

    @Argument(
            doc = "Variance of Gaussian kernel.  If zero, a linear kernel will be used.",
            fullName = KERNEL_VARIANCE_LONG_NAME,
            shortName = KERNEL_VARIANCE_SHORT_NAME,
            optional = true,
            minValue = 0.
    )
    private double kernelVariance = 0.;

    @Argument(
            doc = "Dimension of kernel approximation.  A subsample containing this number of data points " +
                    "will be taken from the copy-ratio profile and used to construct the approximation for each chromosome.  " +
                    "If the total number of datapoints in a chromosome is greater " +
                    "than this number, then all datapoints in the chromosome will be used.",
            fullName = KERNEL_APPROXIMATION_DIMENSION_LONG_NAME,
            shortName = KERNEL_APPROXIMATION_DIMENSION_SHORT_NAME,
            optional = true,
            minValue = 1
    )
    private int kernelApproximationDimension = 100;

    @Argument(
            doc = "Window sizes to use for calculating local changepoint costs.  " +
                    "For each window size, the cost for each datapoint to be a changepoint will be calculated " +
                    "assuming that it demarcates two adjacent segments of that size.  " +
                    "Including small (large) window sizes will increase sensitivity to small (large) events.  " +
                    "Duplicate values will be ignored.",
            fullName = WINDOW_SIZES_LONG_NAME,
            shortName = WINDOW_SIZES_SHORT_NAME,
            optional = true,
            minValue = 1
    )
    private List<Integer> windowSizes;

    @Argument(
            doc = "Linear factor A for the penalty on the number of changepoints per chromosome.  " +
                    "Adds a penalty of the form A * C, where C is the number of changepoints in the chromosome, " +
                    "to the cost function for each chromosome.  " +
                    "Must be either zero or greater than or equal to 1.",
            fullName = NUM_CHANGEPOINTS_PENALTY_LINEAR_FACTOR_LONG_NAME,
            shortName = NUM_CHANGEPOINTS_PENALTY_LINEAR_FACTOR_SHORT_NAME,
            optional = true,
            minValue = 0.
    )
    private double numChangepointsPenaltyLinearFactor = 1.;

    @Argument(
            doc = "Log-linear factor B for the penalty on the number of changepoints per chromosome.  " +
                    "Adds a penalty of the form  B * C * log (N / C), where C is the number of changepoints in the chromosome and " +
                    "N is the number of datapoints in the chromosome, to the cost function for each chromosome.  " +
                    "Must be either zero or greater than or equal to 1.",
            fullName = NUM_CHANGEPOINTS_PENALTY_LOG_LINEAR_FACTOR_LONG_NAME,
            shortName = NUM_CHANGEPOINTS_PENALTY_LOG_LINEAR_FACTOR_SHORT_NAME,
            optional = true,
            minValue = 0.
    )
    private double numChangepointsPenaltyLogLinearFactor = 1.;

    public ModelSegments() {
        windowSizes = Arrays.asList(8, 16, 32, 64, 128, 256);
    }

    @Override
    public Object doWork() {
        //validate arguments
        Utils.validateArg(numChangepointsPenaltyLinearFactor == 0. || numChangepointsPenaltyLinearFactor >= 1.,
                "Linear factor for the penalty on the number of changepoints per chromosome must be either zero or greater than or equal to 1.");
        Utils.validateArg(numChangepointsPenaltyLogLinearFactor == 0. || numChangepointsPenaltyLogLinearFactor >= 1.,
                "Log-linear factor for the penalty on the number of changepoints per chromosome must be either zero or greater than or equal to 1.");

        //TODO clean this up once updated ReadCountCollection is available
        final String sampleName = ReadCountCollectionUtils.getSampleNameForCLIsFromReadCountsFile(new File(inputDenoisedCopyRatioProfileFile));
        final ReadCountCollection denoisedCopyRatioProfile;
        try {
            denoisedCopyRatioProfile = ReadCountCollectionUtils.parse(new File(inputDenoisedCopyRatioProfileFile));
        } catch (final IOException ex) {
            throw new UserException.BadInput("Could not read input file.");
        }


        //segment
        final List<SimpleInterval> segments = new CopyRatioKernelSegmenter(denoisedCopyRatioProfile)
                .findSegments(maxNumSegmentsPerChromosome, kernelVariance, kernelApproximationDimension, ImmutableSet.copyOf(windowSizes).asList(),
                        numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor);

        //TODO add copy-ratio modeller (+ smooth segments)

        //TODO output

        return "SUCCESS";
    }
}