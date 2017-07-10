package org.broadinstitute.hellbender.tools.copynumber.legacy;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5Library;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.pon.coverage.pca.HDF5PCACoveragePoN;
import org.broadinstitute.hellbender.tools.pon.coverage.pca.HDF5PCACoveragePoNCreationUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.OptionalInt;
import java.util.stream.DoubleStream;

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
 * gatk-launch --javaOptions "-Xmx4g" CreatePanelOfNormalsSpark \
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
        summary = "Create a panel of normals (PoN) for copy-ratio denoising given the read counts for samples in the panel.",
        oneLineSummary = "Create a panel of normals for copy-ratio denoising",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public class CreatePanelOfNormalsSpark extends SparkCommandLineProgram {
    //parameter names
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

    private static final String INFER_NUMBER_OF_EIGENSAMPLES = "auto";
    private static final String DEFAULT_NUMBER_OF_EIGENSAMPLES = INFER_NUMBER_OF_EIGENSAMPLES;
    private static final String INTERVAL_WEIGHTS_FILE_SUFFIX = ".interval_weights.txt";

    @Argument(
            doc = "Input integer read-count files for all samples in the panel of normals.  " +
                    "Intervals must be identical for all samples.  " +
                    "Duplicate samples are not removed.",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME
    )
    private List<File> inputReadCountFiles = new ArrayList<>();

    @Argument(
            doc = "Output file name for the panel of normals.  " +
                    "Output is given in the HDF5 file format.  Contents can be viewed with the hdfview program.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputPanelOfNormalsFile;

    @Argument(
            doc = "Genomic intervals with a median coverage (across samples) below this percentile are filtered out.  " +
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
            doc = "Samples with a median of genomic-interval-median-normalized coverage " +
                    "below this percentile or above the complementary percentile are filtered out.  " +
                    "(This is the fourth filter applied.)",
            fullName = EXTREME_SAMPLE_MEDIAN_PERCENTILE_LONG_NAME,
            shortName = EXTREME_SAMPLE_MEDIAN_PERCENTILE_SHORT_NAME,
            minValue = 0.,
            maxValue = 100.
    )
    private double extremeSampleMedianPercentile = DEFAULT_EXTREME_SAMPLE_MEDIAN_PERCENTILE;

    @Argument(
            doc = "Genomic-interval-median-normalized coverages " +
                    "below this percentile or above the complementary percentile are set to the corresponding percentile value.  " +
                    "(This is applied after all filters.)",
            fullName = EXTREME_OUTLIER_TRUNCATION_PERCENTILE_LONG_NAME,
            shortName = EXTREME_OUTLIER_TRUNCATION_PERCENTILE_SHORT_NAME,
            minValue = 0.,
            maxValue = 100.
    )
    private double extremeOutlierTruncationPercentile = DEFAULT_EXTREME_OUTLIER_TRUNCATION_PERCENTILE;

    @Argument(
            doc = "Number of eigensamples to use for denoising.  " +
                    "By default, the number will be automatically inferred using Jolliffe's rule " +
                    "(eigensamples that explain 70% of the total variance will be retained).",
            fullName = NUMBER_OF_EIGENSAMPLES_LONG_NAME,
            shortName = NUMBER_OF_EIGENSAMPLES_SHORT_NAME,
            optional = true
    )
    private String numberOfEigensamplesString = DEFAULT_NUMBER_OF_EIGENSAMPLES;

    @Argument(
            doc = "Output file for the genomic-interval weights (given by the inverse variance of the denoised copy ratio)." +
                    "Output is given as a plain-text file.  " +
                    "By default, " + INTERVAL_WEIGHTS_FILE_SUFFIX + " is appended to the PoN file name.",
            shortName = INTERVAL_WEIGHTS_SHORT_NAME,
            fullName = INTERVAL_WEIGHTS_LONG_NAME,
            optional = true
    )
    private File outputIntervalWeightsFile = null;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        if (!new HDF5Library().load(null)) {  //Note: passing null means using the default temp dir.
            throw new UserException.HardwareFeatureException("Cannot load the required HDF5 library. " +
                    "HDF5 is currently supported on x86-64 architecture and Linux or OSX systems.");
        }

        //parse optional parameters
        if (outputIntervalWeightsFile == null) {
            outputIntervalWeightsFile = new File(outputPanelOfNormalsFile + INTERVAL_WEIGHTS_FILE_SUFFIX);
        }
        final OptionalInt numberOfEigensamples = parseNumberOfEigensamples(numberOfEigensamplesString);

        //parse input read-count files

        //perform optional GC correction

        //create the PoN
        logger.info("Creating the panel of normals...");
        HDF5PCACoveragePoNCreationUtils.create(outputPanelOfNormalsFile, coverageMatrix,
                minimumIntervalMedianPercentile, maximumZerosInSamplePercentage, maximumZerosInIntervalPercentage,
                extremeSampleMedianPercentile, extremeOutlierTruncationPercentile, ctx);

        //output a copy of the interval weights to file
        logger.info("Writing interval-weights file to " + outputIntervalWeightsFile + "...");
        writeIntervalWeightsFile(outputPanelOfNormalsFile, outputIntervalWeightsFile);

        logger.info("Panel of normals successfully created.");
    }

    /**
     * Composes the preferred number of eigenvalues optional given the user input and checks that this is
     * not greater than the number of input samples, if necessary.
     *
     * @return an empty optional if the user elected to use the automatic/inferred value,
     * otherwise a strictly positive integer.
     */
    private OptionalInt parseNumberOfEigensamples(final String numberOfEigensamplesString) {
        if (numberOfEigensamplesString.equalsIgnoreCase(INFER_NUMBER_OF_EIGENSAMPLES)) {
            return OptionalInt.empty();
        } else {
            try {
                final int result = Integer.parseInt(numberOfEigensamplesString);
                Utils.validate(result > 0, NUMBER_OF_EIGENSAMPLES_LONG_NAME + " must be positive.");
                Utils.validate(result <= inputReadCountFiles.size(),
                        String.format("Number of eigensamples cannot be greater than the number of input samples (%d).", inputReadCountFiles.size()));
                return OptionalInt.of(result);
            } catch (final NumberFormatException ex) {
                throw new IllegalArgumentException(NUMBER_OF_EIGENSAMPLES_LONG_NAME + " must be either '" + INFER_NUMBER_OF_EIGENSAMPLES + "' or an integer value");
            }
        }
    }

    /**
     * Read interval variances from an HDF5 PoN file and write the corresponding weights
     * to a file that can be read in by R CBS.
     */
    private static void writeIntervalWeightsFile(final File ponFile, final File outputFile) {
        IOUtils.canReadFile(ponFile);
        try (final HDF5File file = new HDF5File(ponFile, HDF5File.OpenMode.READ_ONLY)) {
            final HDF5PCACoveragePoN pon = new HDF5PCACoveragePoN(file);
            final double[] intervalWeights = DoubleStream.of(pon.getIntervalVariances()).map(v -> 1 / v).toArray();
            ParamUtils.writeValuesToFile(intervalWeights, outputFile);
        }
    }
}
