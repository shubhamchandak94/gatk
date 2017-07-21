package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.rsvd;

import htsjdk.samtools.util.Lazy;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.mllib.linalg.Matrix;
import org.apache.spark.mllib.linalg.SingularValueDecomposition;
import org.apache.spark.mllib.linalg.distributed.RowMatrix;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.spark.SparkConverter;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * TODO
 *
 * HDF5 File backed coverage panel of normals data structure.
 *
 * Several attributes are stored transposed (in other words, the rows and columns are interchanged).
 * This dodges a very slow write time in HDF5, since we usually have many more rows (targets) than columns (samples),
 * and HDF5 writes matrices with few rows and many columns much faster than matrices with many rows and few columns.
 *
 * The following are stored as transposes:
 *
 * <ul>
 *  <li>Normalized Counts</li>
 *  <li>Log-Normalized Counts</li>
 *  <li>Reduced Panel Counts</li>
 *</ul>
 *
 * In these cases, the samples are the rows and the targets are the columns.  No transposing is performed for
 * pseudoinverses since they already have dimensions of samples x targets.
 *
 * This is only for storage.  When saving/loading the above attributes, the transposing is handled transparently.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class HDF5RandomizedSVDReadCountPanelOfNormals implements SVDReadCountPanelOfNormals {
    private static final Logger logger = LogManager.getLogger(HDF5RandomizedSVDReadCountPanelOfNormals.class);

    private static final int NUM_SLICES_FOR_SPARK_MATRIX_CONVERSION = 100;
    private static final double EPSILON = 1E-9;

    private final HDF5File file;

    /**
     * The version number is a double where the integer part is the
     * major version and the decimal part is the minor version.
     * The minor version should be only a single digit.
     */
    private static final double CURRENT_PON_VERSION = 7.0;
    private static final String PON_VERSION_STRING_FORMAT = "%.1f";

    private static final String VERSION_PATH = "/version/value";    //note that full path names must include a top-level group name ("version" here)
    private static final String COMMAND_LINE_PATH = "/command_line/value";
    private static final String ORIGINAL_DATA_GROUP_NAME = "/original_data";
    private static final String ORIGINAL_READ_COUNTS_PATH = ORIGINAL_DATA_GROUP_NAME + "/transposed_read_counts";
    private static final String ORIGINAL_SAMPLE_FILENAMES_PATH = ORIGINAL_DATA_GROUP_NAME + "/sample_filenames";
    private static final String ORIGINAL_INTERVALS_PATH = ORIGINAL_DATA_GROUP_NAME + "/intervals";
    private static final String ORIGINAL_INTERVAL_GC_CONTENT_PATH = ORIGINAL_DATA_GROUP_NAME + "/interval_gc_content";

    private static final String PANEL_GROUP_NAME = "/panel";
    private static final String PANEL_SAMPLE_FILENAMES_PATH = PANEL_GROUP_NAME + "/sample_filenames";
    private static final String PANEL_INTERVALS_PATH = PANEL_GROUP_NAME + "/intervals";
    private static final String PANEL_INTERVAL_FRACTIONAL_MEDIANS_PATH = PANEL_GROUP_NAME + "/interval_fractional_medians";
    private static final String PANEL_SINGULAR_VALUES_PATH = PANEL_GROUP_NAME + "/singular_values";
    private static final String PANEL_LEFT_SINGULAR_PATH = PANEL_GROUP_NAME + "/transposed_left_singular";

    private final Lazy<List<Locatable>> originalIntervals;
    private final Lazy<List<Locatable>> panelIntervals;

    /**
     * <p>DEV NOTE:  If you are adding attributes that are neither RealMatrix nor a primitive,
     * you must follow the pattern in the constructor (i.e. the Lazy loading pattern).
     * Otherwise, some operations will hang.</p>
     */
    private HDF5RandomizedSVDReadCountPanelOfNormals(final HDF5File file) {
        Utils.nonNull(file);
        this.file = file;
        originalIntervals = new Lazy<>(() -> IntervalHelper.readIntervals(file, ORIGINAL_INTERVALS_PATH));
        panelIntervals = new Lazy<>(() -> IntervalHelper.readIntervals(file, PANEL_INTERVALS_PATH));
    }
    
    @Override
    public double getVersion() {
        if (!file.isPresent(VERSION_PATH)) {    //the version path may be different in older PoNs
            throw new UserException.BadInput(String.format("The panel of normals is out of date and incompatible.  " +
                    "Please use a panel of normals that was created by CreateReadCountPanelOfNormals and is version " +
                    PON_VERSION_STRING_FORMAT + ".", CURRENT_PON_VERSION));
        }
        return file.readDouble(VERSION_PATH);
    }

    @Override
    public int getNumEigensamples() {
        return file.readStringArray(PANEL_SAMPLE_FILENAMES_PATH).length;
    }

    @Override
    public double[][] getOriginalReadCounts() {
        return new Array2DRowRealMatrix(file.readDoubleMatrix(ORIGINAL_READ_COUNTS_PATH), false).transpose().getData();
    }

    @Override
    public List<Locatable> getOriginalIntervals() {
        return originalIntervals.get();
    }

    @Override
    public double[] getOriginalIntervalGCContent() {
        if (!file.isPresent(ORIGINAL_INTERVAL_GC_CONTENT_PATH)) {
            return null;
        }
        return file.readDoubleArray(ORIGINAL_INTERVAL_GC_CONTENT_PATH);
    }

    @Override
    public List<Locatable> getPanelIntervals() {
        return panelIntervals.get();
    }

    @Override
    public double[] getPanelIntervalFractionalMedians() {
        return file.readDoubleArray(PANEL_INTERVAL_FRACTIONAL_MEDIANS_PATH);
    }

    @Override
    public double[] getSingularValues() {
        return file.readDoubleArray(PANEL_SINGULAR_VALUES_PATH);
    }

    @Override
    public double[][] getLeftSingular() {
        return new Array2DRowRealMatrix(file.readDoubleMatrix(PANEL_LEFT_SINGULAR_PATH), false).transpose().getData();
    }

    /**
     * Create an interface to an HDF5 file.  A version check is performed and a warning message logged if the
     * version number is not up to date.
     */
    public static HDF5RandomizedSVDReadCountPanelOfNormals read(final HDF5File file) {
        final HDF5RandomizedSVDReadCountPanelOfNormals pon = new HDF5RandomizedSVDReadCountPanelOfNormals(file);
        if (pon.getVersion() < CURRENT_PON_VERSION) {
            throw new UserException.BadInput(String.format("The version of the specified panel of normals (%f) is older than the current version (%f).",
                    pon.getVersion(), CURRENT_PON_VERSION));
        }
        return pon;
    }

    /**
     * Create the panel of normals and write it to an HDF5 file.  All inputs are assumed to be valid.
     * The dimensions of {@code originalReadCounts} should be samples x intervals.
     * To reduce memory footprint, {@code originalReadCounts} is modified in place.
     * If {@code intervalGCContent} is null, GC-bias correction will not be performed.
     */
    public static void create(final File outFile,
                              final String commandLine,
                              final RealMatrix originalReadCounts,
                              final List<String> originalSampleFilenames,
                              final List<Locatable> originalIntervals,
                              final double[] intervalGCContent,
                              final double minimumIntervalMedianPercentile,
                              final double maximumZerosInSamplePercentage,
                              final double maximumZerosInIntervalPercentage,
                              final double extremeSampleMedianPercentile,
                              final double extremeOutlierTruncationPercentile,
                              final int numEigensamplesRequested,
                              final JavaSparkContext ctx) {
        try (final HDF5File file = new HDF5File(outFile, HDF5File.OpenMode.CREATE)) {
            logger.info("Creating " + outFile.getAbsolutePath() + "...");
            final HDF5RandomizedSVDReadCountPanelOfNormals pon = new HDF5RandomizedSVDReadCountPanelOfNormals(file);

            logger.info(String.format("Writing version number (" + PON_VERSION_STRING_FORMAT + ")...", CURRENT_PON_VERSION));
            pon.writeVersion(CURRENT_PON_VERSION);

            logger.info("Writing command line...");
            pon.writeCommandLine(commandLine);

            logger.info(String.format("Writing original read counts (transposed to %d x %d)...",
                    originalReadCounts.getColumnDimension(), originalReadCounts.getRowDimension()));
            pon.writeOriginalReadCountsPath(originalReadCounts);

            logger.info(String.format("Writing original sample filenames (%d)...", originalSampleFilenames.size()));
            pon.writeOriginalSampleFilenames(originalSampleFilenames);

            logger.info(String.format("Writing original intervals (%d)...", originalIntervals.size()));
            pon.writeOriginalIntervals(originalIntervals);

            if (intervalGCContent != null) {
                logger.info(String.format("Writing GC-content annotations for original intervals (%d)...", intervalGCContent.length));
                pon.writeOriginalIntervalGCContent(intervalGCContent);
            }

            //preprocess and standardize read counts and determine filters
            //(originalReadCounts is modified in place and a filtered submatrix is returned)
            logger.info("Preprocessing and standardizing read counts...");
            final SVDDenoisingUtils.PreprocessedStandardizedResult preprocessedStandardizedResult =
                    SVDDenoisingUtils.preprocessAndStandardize(originalReadCounts, intervalGCContent,
                    minimumIntervalMedianPercentile, maximumZerosInSamplePercentage, maximumZerosInIntervalPercentage,
                    extremeSampleMedianPercentile, extremeOutlierTruncationPercentile);

            //filter samples and intervals
            final List<String> panelSampleFilenames = IntStream.range(0, originalSampleFilenames.size())
                    .filter(sampleIndex -> !preprocessedStandardizedResult.filterSamples[sampleIndex])
                    .mapToObj(originalSampleFilenames::get).collect(Collectors.toList());
            final List<Locatable> panelIntervals = IntStream.range(0, originalIntervals.size())
                    .filter(intervalIndex -> !preprocessedStandardizedResult.filterIntervals[intervalIndex])
                    .mapToObj(originalIntervals::get).collect(Collectors.toList());

            logger.info(String.format("Writing panel sample filenames (%d)...", panelSampleFilenames.size()));
            pon.writePanelSampleFilenames(panelSampleFilenames);

            logger.info(String.format("Writing panel intervals (%d)...", panelIntervals.size()));
            pon.writePanelIntervals(panelIntervals);

            //get panel interval fractional medians (calculated as an intermediate result during preprocessing)
            final double[] panelIntervalFractionalMedians = preprocessedStandardizedResult.panelIntervalFractionalMedians;

            logger.info(String.format("Writing panel interval fractional medians (%d)...", panelIntervalFractionalMedians.length));
            pon.writePanelIntervalFractionalMedians(panelIntervalFractionalMedians);

            final int numPanelIntervals = preprocessedStandardizedResult.preprocessedStandardizedProfile.getRowDimension();
            final int numPanelSamples = preprocessedStandardizedResult.preprocessedStandardizedProfile.getColumnDimension();
            final int numEigensamples = Math.min(numEigensamplesRequested, numPanelSamples);
            if (numEigensamples < numEigensamplesRequested) {
                logger.warn(String.format("%d eigensamples were requested but only %d are available in the panel of normals...",
                        numEigensamplesRequested, numEigensamples));
            }
            logger.info(String.format("Performing SVD (truncated at %d eigensamples) of standardized counts (%d x %d)...",
                    Math.min(numEigensamples, numPanelSamples), numPanelIntervals, numPanelSamples));
            final SingularValueDecomposition<RowMatrix, Matrix> svd = SparkConverter.convertRealMatrixToSparkRowMatrix(
                    ctx, preprocessedStandardizedResult.preprocessedStandardizedProfile, NUM_SLICES_FOR_SPARK_MATRIX_CONVERSION)
                    .computeSVD(numEigensamples, true, EPSILON);
            final double[] singularValues = svd.s().toArray();    //should be in decreasing order (with corresponding matrices below)
            final double[][] leftSingular = SparkConverter.convertSparkRowMatrixToRealMatrix(svd.U(), numPanelIntervals).getData();

            logger.info(String.format("Writing singular values (%d)...", singularValues.length));
            pon.writeSingularValues(singularValues);

            logger.info(String.format("Writing left-singular matrix (transposed to %d x %d)...", leftSingular[0].length, leftSingular.length));
            pon.writeLeftSingular(leftSingular);
        } catch (final RuntimeException e) {
            //if any exceptions encountered, delete partial output and rethrow
            logger.warn(String.format("Exception encountered during creation of panel of normals.  Attempting to delete partial output in %s...",
                    outFile.getAbsolutePath()));
            IOUtils.tryDelete(outFile);
            logger.warn("Partial output deleted.");
            throw e;
        }
        logger.info(String.format("Read-count panel of normals written to %s.", outFile));
    }

    //PRIVATE WRITERS (write values to HDF5 file)
    //these are private to prevent fields from being written individually, which could leave the file in a bad state

    private void writeVersion(final double version) {
        file.makeDouble(VERSION_PATH, version);
    }

    private void writeCommandLine(final String commandLine) {
        file.makeStringArray(COMMAND_LINE_PATH, commandLine);
    }

    private void writeOriginalReadCountsPath(final RealMatrix originalReadCounts) {
        file.makeDoubleMatrix(ORIGINAL_READ_COUNTS_PATH, originalReadCounts.transpose().getData());
    }

    private void writeOriginalSampleFilenames(final List<String> originalSampleFilenames) {
        file.makeStringArray(ORIGINAL_SAMPLE_FILENAMES_PATH, originalSampleFilenames.toArray(new String[originalSampleFilenames.size()]));
    }

    private void writeOriginalIntervals(final List<Locatable> originalIntervals) {
        IntervalHelper.writeIntervals(file, ORIGINAL_INTERVALS_PATH, originalIntervals);
    }

    private void writeOriginalIntervalGCContent(final double[] originalIntervalGCContent) {
        file.makeDoubleArray(ORIGINAL_INTERVAL_GC_CONTENT_PATH, originalIntervalGCContent);
    }

    private void writePanelSampleFilenames(final List<String> panelSampleFilenames) {
        file.makeStringArray(PANEL_SAMPLE_FILENAMES_PATH, panelSampleFilenames.toArray(new String[panelSampleFilenames.size()]));
    }

    private void writePanelIntervals(final List<Locatable> panelIntervals) {
        IntervalHelper.writeIntervals(file, PANEL_INTERVALS_PATH, panelIntervals);
    }

    private void writePanelIntervalFractionalMedians(final double[] panelIntervalFractionalMedians) {
        file.makeDoubleArray(PANEL_INTERVAL_FRACTIONAL_MEDIANS_PATH, panelIntervalFractionalMedians);
    }

    private void writeSingularValues(final double[] singularValues) {
        file.makeDoubleArray(PANEL_SINGULAR_VALUES_PATH, singularValues);
    }

    private void writeLeftSingular(final double[][] leftSingular) {
        file.makeDoubleMatrix(PANEL_LEFT_SINGULAR_PATH, new Array2DRowRealMatrix(leftSingular, false).transpose().getData());
    }

    private static final class IntervalHelper {
        //writing intervals as a string matrix is expensive,
        //so we instead store a map from integer indices to contig strings and
        //store (index, start, end) in a double matrix

        private static final String INTERVAL_CONTIG_NAMES_PATH_SUFFIX = "/contig_names";
        private static final String INTERVAL_MATRIX_PATH_SUFFIX = "/contig-index_start_end";

        private enum IntervalField {
            CONTIG_INDEX(0),
            START (1),
            END (2);

            private final int index;

            IntervalField(final int index) {
                this.index = index;
            }
        }
        private static final int NUM_INTERVAL_FIELDS = IntervalField.values().length;

        private static List<Locatable> readIntervals(final HDF5File file,
                                                     final String path) {
            final String[] contigNames = file.readStringArray(path + INTERVAL_CONTIG_NAMES_PATH_SUFFIX);
            final double[][] matrix = file.readDoubleMatrix(path + INTERVAL_MATRIX_PATH_SUFFIX);
            final int numIntervals = matrix[0].length;
            return IntStream.range(0, numIntervals).boxed()
                    .map(i -> (new SimpleInterval(
                            contigNames[(int) matrix[IntervalField.CONTIG_INDEX.index][i]],
                            (int) matrix[IntervalField.START.index][i],
                            (int) matrix[IntervalField.END.index][i])))
                    .collect(Collectors.toList());
        }

        private static void writeIntervals(final HDF5File file,
                                           final String path,
                                           final List<Locatable> intervals) {
            final String[] contigNames = intervals.stream().map(Locatable::getContig).distinct().toArray(String[]::new);
            file.makeStringArray(path + INTERVAL_CONTIG_NAMES_PATH_SUFFIX, contigNames);
            final Map<String, Double> contigNamesToIndexMap = IntStream.range(0, contigNames.length).boxed()
                    .collect(Collectors.toMap(i -> contigNames[i], i -> (double) i));
            final double[][] matrix = new double[NUM_INTERVAL_FIELDS][intervals.size()];
            for (int i = 0; i < intervals.size(); i++) {
                final Locatable interval = intervals.get(i);
                matrix[IntervalField.CONTIG_INDEX.index][i] = contigNamesToIndexMap.get(interval.getContig());
                matrix[IntervalField.START.index][i] = interval.getStart();
                matrix[IntervalField.END.index][i] = interval.getEnd();
            }
            file.makeDoubleMatrix(path + INTERVAL_MATRIX_PATH_SUFFIX, matrix);
        }
    }
}
