package org.broadinstitute.hellbender.tools.copynumber.legacy;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hdf5.HDF5Library;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.rsvd.SVDDenoisingUtils;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;

/**
 * //TODO
 * Tool to create a panel of normals (PoN) given read counts for control samples.
 *
 * <p>
 * The input read counts are... 
 * </p>
 *
 * <h3>Examples</h3>
 *
 * <p>
 *     The following command is for either whole exome sequencing (WES) or whole genome sequencing (WGS) data.
 * </p>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" CreateReadCountPanelOfNormals \
 *   --input gc_corrected_coverages.tsv \
 *   --output panel_of_normals.pon
 * </pre>
 *
 * <p>
 * In addition to the resulting PoN, this command produces a .pon.removed_samples.txt file of samples removed for quality control (QC)
 * and a .pon.target_weights.txt file that gives the inverse variance per target that can optionally be passed to PerformSegmentation.
 * </p>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Create a panel of normals for copy-ratio denoising given the read counts for samples in the panel.",
        oneLineSummary = "Create a panel of normals for copy-ratio denoising",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
//public class CreateReadCountPanelOfNormals extends SparkCommandLineProgram {
public class CreateReadCountPanelOfNormals extends CommandLineProgram {
    private static final long serialVersionUID = 1L;

    //parameter names
    private static final String GC_ANNOTATED_INTERVAL_FILE_LONG_NAME = "gcAnnotatedIntervals";
    private static final String GC_ANNOTATED_INTERVAL_FILE_SHORT_NAME = "gcAnnot";
    private static final String MINIMUM_INTERVAL_MEDIAN_PERCENTILE_LONG_NAME = "minimumIntervalMedianPercentile";
    private static final String MINIMUM_INTERVAL_MEDIAN_PERCENTILE_SHORT_NAME = "minIntervalMedPct";
    private static final String MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE_LONG_NAME = "maximumZerosInSamplePercentage";
    private static final String MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE_SHORT_NAME = "maxZerosInSamplePct";
    private static final String MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE_LONG_NAME = "maximumZerosInIntervalPercentage";
    private static final String MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE_SHORT_NAME = "maxZerosInIntervalPct";
    private static final String EXTREME_SAMPLE_MEDIAN_PERCENTILE_LONG_NAME = "extremeSampleMedianPercentile";
    private static final String EXTREME_SAMPLE_MEDIAN_PERCENTILE_SHORT_NAME = "extSampleMedPct";
    private static final String EXTREME_OUTLIER_TRUNCATION_PERCENTILE_LONG_NAME = "extremeOutlierTruncationPercentile";
    private static final String EXTREME_OUTLIER_TRUNCATION_PERCENTILE_SHORT_NAME = "extOutTruncPct";
    private static final String NUMBER_OF_EIGENSAMPLES_LONG_NAME = "numberOfEigensamples";
    private static final String NUMBER_OF_EIGENSAMPLES_SHORT_NAME = "numEigen";
    private static final String INTERVAL_WEIGHTS_LONG_NAME = "intervalWeights";
    private static final String INTERVAL_WEIGHTS_SHORT_NAME = "weights";

    //default values for filtering (taken from ReCapSeg)
    private static final double DEFAULT_MINIMUM_INTERVAL_MEDIAN_PERCENTILE = 25.0;
    private static final double DEFAULT_MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE = 2.0;
    private static final double DEFAULT_MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE = 5.0;
    private static final double DEFAULT_EXTREME_SAMPLE_MEDIAN_PERCENTILE = 2.5;
    private static final double DEFAULT_EXTREME_OUTLIER_TRUNCATION_PERCENTILE = 0.1;
    private static final int DEFAULT_NUMBER_OF_EIGENSAMPLES = 10;

    private static final String INTERVAL_WEIGHTS_FILE_SUFFIX = ".interval_weights.txt";

    @Argument(
            doc = "Input read-count files containing integer read counts in genomic intervals for all samples in the panel of normals.  " +
                    "Intervals must be identical and in the same order for all samples.  " +
                    "Duplicate samples are not removed.",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            minElements = 1
    )
    private List<File> inputReadCountFiles = new ArrayList<>();

    @Argument(
            doc = "Input annotated-interval file containing annotations for GC content in genomic intervals (output of AnnotateTargets).  " +
                    "Intervals must be identical to and in the same order as those in the input read-count files.",
            fullName = GC_ANNOTATED_INTERVAL_FILE_LONG_NAME,
            shortName = GC_ANNOTATED_INTERVAL_FILE_SHORT_NAME,
            optional = true
    )
    private File inputGCAnnotatedIntervals = null;

    @Argument(
            doc = "Output file name for the panel of normals.  " +
                    "Output is given in the HDF5 file format.  Contents can be viewed with the hdfview program.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputPanelOfNormalsFile;

    @Argument(
            doc = "Genomic intervals with a median (across samples) of fractional coverage below this percentile are filtered out.  " +
                    "(This is the first filter applied.)",
            fullName  = MINIMUM_INTERVAL_MEDIAN_PERCENTILE_LONG_NAME,
            shortName = MINIMUM_INTERVAL_MEDIAN_PERCENTILE_SHORT_NAME,
            minValue = 0.,
            maxValue = 100.
    )
    private double minimumIntervalMedianPercentile = DEFAULT_MINIMUM_INTERVAL_MEDIAN_PERCENTILE;

    @Argument(
            doc = "Samples with an amount of zero-coverage genomic intervals above this percentage are filtered out.  " +
                    "(This is the second filter applied.)",
            fullName = MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE_LONG_NAME,
            shortName = MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE_SHORT_NAME,
            minValue = 0.,
            maxValue = 100.
    )
    private double maximumZerosInSamplePercentage = DEFAULT_MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE;

    @Argument(
            doc = "Genomic intervals with an amount of zero-coverage samples above this percentage are filtered out.  " +
                    "(This is the third filter applied.)",
            fullName = MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE_LONG_NAME,
            shortName = MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE_SHORT_NAME,
            minValue = 0.,
            maxValue = 100.
    )
    private double maximumZerosInIntervalPercentage = DEFAULT_MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE;

    @Argument(
            doc = "Samples with a median (across intervals) of fractional coverage normalized by genomic-interval medians  " +
                    "below this percentile or above the complementary percentile are filtered out.  " +
                    "(This is the fourth filter applied.)",
            fullName = EXTREME_SAMPLE_MEDIAN_PERCENTILE_LONG_NAME,
            shortName = EXTREME_SAMPLE_MEDIAN_PERCENTILE_SHORT_NAME,
            minValue = 0.,
            maxValue = 100.
    )
    private double extremeSampleMedianPercentile = DEFAULT_EXTREME_SAMPLE_MEDIAN_PERCENTILE;

    @Argument(
            doc = "Fractional coverages normalized by genomic-interval medians that are " +
                    "below this percentile or above the complementary percentile are set to the corresponding percentile value.  " +
                    "(This is applied after all filters.)",
            fullName = EXTREME_OUTLIER_TRUNCATION_PERCENTILE_LONG_NAME,
            shortName = EXTREME_OUTLIER_TRUNCATION_PERCENTILE_SHORT_NAME,
            minValue = 0.,
            maxValue = 100.
    )
    private double extremeOutlierTruncationPercentile = DEFAULT_EXTREME_OUTLIER_TRUNCATION_PERCENTILE;

    @Argument(
            doc = "Number of eigensamples to use for randomized SVD and to store in the panel of normals.  " +
                    "The number of samples retained after filtering will be used instead if it is smaller than this.",
            fullName = NUMBER_OF_EIGENSAMPLES_LONG_NAME,
            shortName = NUMBER_OF_EIGENSAMPLES_SHORT_NAME,
            minValue = 1
    )
    private int numberOfEigensamples = DEFAULT_NUMBER_OF_EIGENSAMPLES;

    @Argument(
            doc = "Output file for the genomic-interval weights (given by the inverse variance of the denoised copy ratio)." +
                    "Output is given as a plain-text file.  " +
                    "By default, " + INTERVAL_WEIGHTS_FILE_SUFFIX + " is appended to the panel-of-normals filename.",
            shortName = INTERVAL_WEIGHTS_SHORT_NAME,
            fullName = INTERVAL_WEIGHTS_LONG_NAME,
            optional = true
    )
    private File outputIntervalWeightsFile = null;

    @Override
//    protected void runPipeline(final JavaSparkContext ctx) {
    protected Object doWork() {
        if (!new HDF5Library().load(null)) {  //Note: passing null means using the default temp dir.
            throw new UserException.HardwareFeatureException("Cannot load the required HDF5 library. " +
                    "HDF5 is currently supported on x86-64 architecture and Linux or OSX systems.");
        }

        //validate parameters and parse optional parameters
        validateArguments();
        //TODO check for GC annotations
        if (outputIntervalWeightsFile == null) {
            outputIntervalWeightsFile = new File(outputPanelOfNormalsFile + INTERVAL_WEIGHTS_FILE_SUFFIX);
        }
        final List<Target> intervals = getIntervalsFromFirstReadCountFile();

        //validate input read-count files (i.e., check intervals and that only integer counts are contained) and aggregate as a RealMatrix
        logger.info("Validing and aggregating input read-count files...");
        final int numSamples = inputReadCountFiles.size();
        final int numIntervals = intervals.size();
        final RealMatrix readCountsMatrix = new Array2DRowRealMatrix(numSamples, numIntervals);
        final ListIterator<File> inputReadCountFilesIterator = inputReadCountFiles.listIterator();
        while (inputReadCountFilesIterator.hasNext()) {
            final int sampleIndex = inputReadCountFilesIterator.nextIndex();
            final File inputReadCountFile = inputReadCountFilesIterator.next();
            logger.info(String.format("Aggregating read-count file %s (%d / %d)", inputReadCountFile, sampleIndex + 1, numSamples));
            final ReadCountCollection readCounts;
            try {
                readCounts = ReadCountCollectionUtils.parse(inputReadCountFile);
            } catch (final IOException e) {
                throw new UserException.CouldNotReadInputFile(inputReadCountFile);
            }
            SVDDenoisingUtils.validateReadCounts(readCounts);
            Utils.validateArg(readCounts.targets().equals(intervals),
                    String.format("Intervals for read-count file %s do not match those in other read-count files.", inputReadCountFile));
            readCountsMatrix.setRow(sampleIndex, readCounts.getColumn(0));
        }

//        //create the PoN
//        logger.info("Creating the panel of normals...");
//        HDF5PCACoveragePoNCreationUtils.create(outputPanelOfNormalsFile, getCommandLine(), coverageMatrix, gcAnnotatedIntervals,
//                minimumIntervalMedianPercentile, maximumZerosInSamplePercentage, maximumZerosInIntervalPercentage,
//                extremeSampleMedianPercentile, extremeOutlierTruncationPercentile, ctx);
//
//        //output a copy of the interval weights to file
//        logger.info("Writing interval-weights file to " + outputIntervalWeightsFile + "...");
//        writeIntervalWeightsFile(outputPanelOfNormalsFile, outputIntervalWeightsFile);

        logger.info("Panel of normals successfully created.");
        return "SUCCESS";
    }

    private List<Target> getIntervalsFromFirstReadCountFile() {
        final File firstReadCountFile = inputReadCountFiles.get(0);
        logger.info(String.format("Retrieving intervals from first read-count file (%s)...", firstReadCountFile));
        final ReadCountCollection firstSampleReadCounts;
        try {
            firstSampleReadCounts = ReadCountCollectionUtils.parse(firstReadCountFile);
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(firstReadCountFile);
        }
        return firstSampleReadCounts.targets();
    }

    private void validateArguments() {
        inputReadCountFiles.forEach(IOUtils::canReadFile);
        Utils.validateArg(numberOfEigensamples <= inputReadCountFiles.size(),
                String.format("Number of eigensamples cannot be greater than the number of input samples (%d).", inputReadCountFiles.size()));
    }

//    /**
//     * Read interval variances from an HDF5 PoN file and write the corresponding weights
//     * to a file that can be read in by R CBS.
//     */
//    private static void writeIntervalWeightsFile(final File ponFile, final File outputFile) {
//        IOUtils.canReadFile(ponFile);
//        try (final HDF5File file = new HDF5File(ponFile, HDF5File.OpenMode.READ_ONLY)) {
//            final HDF5PCACoveragePoN pon = new HDF5PCACoveragePoN(file);
//            final double[] intervalWeights = DoubleStream.of(pon.getIntervalVariances()).map(v -> 1 / v).toArray();
//            ParamUtils.writeValuesToFile(intervalWeights, outputFile);
//        }
//    }
}
