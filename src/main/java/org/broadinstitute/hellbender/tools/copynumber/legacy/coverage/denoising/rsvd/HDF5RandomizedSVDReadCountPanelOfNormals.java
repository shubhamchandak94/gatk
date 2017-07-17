package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.rsvd;

import com.google.common.primitives.Doubles;
import htsjdk.samtools.util.Lazy;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.gcbias.GCCorrector;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MatrixSummaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.svd.SVD;
import org.broadinstitute.hellbender.utils.svd.SVDFactory;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
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

    private final HDF5File file;

    /**
     * The version number is a double where the integer part is the
     * major version and the decimal part is the minor version.
     * The minor version should be only a single digit.
     */
    private static final double CURRENT_PON_VERSION = 7.0;

    private static final String VERSION_PATH = "/version/value";    //note that full path names must include a top-level group name ("version" here)
    private static final String COMMAND_LINE_PATH = "/command_line/value";
    private static final String ORIGINAL_DATA_GROUP_NAME = "/original_data";
    private static final String ORIGINAL_READ_COUNTS_PATH = ORIGINAL_DATA_GROUP_NAME + "/read_counts";
    private static final String ORIGINAL_SAMPLE_FILENAMES_PATH = ORIGINAL_DATA_GROUP_NAME + "/sample_filenames";
    private static final String ORIGINAL_INTERVALS_PATH = ORIGINAL_DATA_GROUP_NAME + "/intervals";
    private static final String ORIGINAL_INTERVAL_GC_CONTENT_PATH = ORIGINAL_DATA_GROUP_NAME + "/interval_gc_content";

    private static final String PANEL_GROUP_NAME = "/panel";
    private static final String PANEL_SAMPLE_FILENAMES_PATH = PANEL_GROUP_NAME + "/sample_filenames";
    private static final String PANEL_INTERVALS_PATH = PANEL_GROUP_NAME + "/intervals";
    private static final String PANEL_INTERVAL_FRACTIONAL_MEDIANS_PATH = PANEL_GROUP_NAME + "/interval_fractional_medians";
    private static final String PANEL_SINGULAR_VALUES_PATH = PANEL_GROUP_NAME + "/singular_values";
    private static final String PANEL_LEFT_SINGULAR_PATH = PANEL_GROUP_NAME + "/left_singular";
    private static final String PANEL_LEFT_SINGULAR_PSEUDOINVERSE_PATH = PANEL_GROUP_NAME + "/left_singular_pseudoinverse";

    private static final String NUM_INTERVAL_COLUMNS_PATH = "/num_interval_columns/value";

    private enum IntervalColumn {
        CONTIG (0),
        START (1),
        END (2);

        private final int index;

        IntervalColumn(final int index) {
            this.index = index;
        }
    }
    private static final int NUM_INTERVAL_COLUMNS = IntervalColumn.values().length;

    private final Lazy<List<SimpleInterval>> originalIntervals;
    private final Lazy<List<SimpleInterval>> panelIntervals;

    /**
     * <p>DEV NOTE:  If you are adding attributes that are neither RealMatrix nor a primitive,
     * you must follow the pattern in the constructor (i.e. the Lazy loading pattern).
     * Otherwise, some operations will hang.</p>
     */
    private HDF5RandomizedSVDReadCountPanelOfNormals(final HDF5File file) {
        Utils.nonNull(file);
        this.file = file;
        originalIntervals = new Lazy<>(() -> readIntervals(file, ORIGINAL_INTERVALS_PATH));
        panelIntervals = new Lazy<>(() -> readIntervals(file, PANEL_INTERVALS_PATH));
    }
    
    @Override
    public double getVersion() {
        return file.readDouble(VERSION_PATH);
    }

    @Override
    public int getNumEigensamples() {
        return file.readStringArray(PANEL_SAMPLE_FILENAMES_PATH).length;
    }

    @Override
    public RealMatrix getOriginalReadCounts() {
        return new Array2DRowRealMatrix(file.readDoubleMatrix(ORIGINAL_READ_COUNTS_PATH));
    }

    @Override
    public List<SimpleInterval> getOriginalIntervals() {
        return originalIntervals.get();
    }

    @Override
    public double[] getOriginalIntervalGCContent() {
        try {
            return file.readDoubleArray(ORIGINAL_INTERVAL_GC_CONTENT_PATH);
        } catch (final IllegalArgumentException e) {
            return null;
        }
    }

    @Override
    public List<SimpleInterval> getPanelIntervals() {
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
        return file.readDoubleMatrix(PANEL_LEFT_SINGULAR_PATH);
    }

    @Override
    public double[][] getLeftSingularPseudoinverse() {
        return file.readDoubleMatrix(PANEL_LEFT_SINGULAR_PSEUDOINVERSE_PATH);
    }

    /**
     * Create an interface to an HDF5 file.  A version check is performed and a warning message logged if the
     * version number is not up to date.
     */
    public static HDF5RandomizedSVDReadCountPanelOfNormals read(final HDF5File file) {
        final HDF5RandomizedSVDReadCountPanelOfNormals pon = new HDF5RandomizedSVDReadCountPanelOfNormals(file);
        if (pon.getVersion() < CURRENT_PON_VERSION) {
            logger.warn(String.format("The version of the specified panel of normals (%f) is older than the latest version (%f).",
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
                              final List<SimpleInterval> originalIntervals,
                              final double[] intervalGCContent,
                              final double minimumIntervalMedianPercentile,
                              final double maximumZerosInSamplePercentage,
                              final double maximumZerosInIntervalPercentage,
                              final double extremeSampleMedianPercentile,
                              final double extremeOutlierTruncationPercentile,
                              final JavaSparkContext ctx) {
        try (final HDF5File file = new HDF5File(outFile, HDF5File.OpenMode.CREATE)) {
            logger.info("Creating " + outFile.getAbsolutePath() + "...");
            final HDF5RandomizedSVDReadCountPanelOfNormals pon = new HDF5RandomizedSVDReadCountPanelOfNormals(file);

            logger.info(String.format("Writing version number (%.1f)...", CURRENT_PON_VERSION));
            pon.writeVersion(CURRENT_PON_VERSION);

            logger.info("Writing command line...");
            pon.writeCommandLine(commandLine);

            logger.info(String.format("Writing original read counts (%d x %d)...", originalReadCounts.getRowDimension(), originalReadCounts.getColumnDimension()));
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
            logger.info("Preprocessing and standardizing read counts...");
            final PreprocessedStandardizedResult preprocessedStandardizedResult =
                    preprocessAndStandardize(originalReadCounts, intervalGCContent,
                    minimumIntervalMedianPercentile, maximumZerosInSamplePercentage, maximumZerosInIntervalPercentage,
                    extremeSampleMedianPercentile, extremeOutlierTruncationPercentile);

            //filter samples and intervals
            final List<String> panelSampleFilenames = IntStream.range(0, originalSampleFilenames.size())
                    .filter(i -> !preprocessedStandardizedResult.filterSamples[i])
                    .mapToObj(originalSampleFilenames::get).collect(Collectors.toList());
            final List<SimpleInterval> panelIntervals = IntStream.range(0, originalIntervals.size())
                    .filter(i -> !preprocessedStandardizedResult.filterIntervals[i])
                    .mapToObj(originalIntervals::get).collect(Collectors.toList());

            logger.info(String.format("Writing panel sample filenames (%d)...", panelSampleFilenames.size()));
            pon.writePanelSampleFilenames(panelSampleFilenames);

            logger.info(String.format("Writing panel intervals (%d)...", panelIntervals.size()));
            pon.writePanelIntervals(panelIntervals);

            //get panel interval fractional medians (calculated as an intermediate result during preprocessing)
            final double[] panelIntervalFractionalMedians = preprocessedStandardizedResult.panelIntervalFractionalMedians;

            logger.info(String.format("Writing panel interval fractional medians (%d)...", panelIntervalFractionalMedians.length));
            pon.writePanelIntervalFractionalMedians(panelIntervalFractionalMedians);

            logger.info("Performing SVD...");
            final SVD svd = SVDFactory.createSVD(preprocessedStandardizedResult.preprocessedStandardizedReadCounts, ctx);
            final double[] singularValues = svd.getSingularValues();
            final double[][] leftSingular = svd.getU().getData();
            final double[][] leftSingularPseudoinverse = svd.getPinv().getData();

            logger.info(String.format("Writing singular values (%d)...", singularValues.length));
            pon.writeSingularValues(singularValues);

            logger.info(String.format("Writing left-singular matrix (%d x %d)...", leftSingular.length, leftSingular[0].length));
            pon.writeLeftSingular(leftSingular);

            logger.info(String.format("Writing left-singular pseudoinverse (%d x %d)...", leftSingularPseudoinverse.length, leftSingularPseudoinverse[0].length));
            pon.writePanelLeftSingularPseudoinverse(leftSingularPseudoinverse);
        }
        logger.info(String.format("Read-count panel of normals written to %s.", outFile));
    }

    private static final class PreprocessedStandardizedResult {
        private final RealMatrix preprocessedStandardizedReadCounts;
        private final double[] panelIntervalFractionalMedians;
        private final boolean[] filterIntervals;
        private final boolean[] filterSamples;

        private PreprocessedStandardizedResult(final RealMatrix preprocessedStandardizedReadCounts,
                                               final double[] panelIntervalFractionalMedians,
                                               final boolean[] filterIntervals,
                                               final boolean[] filterSamples) {
            this.preprocessedStandardizedReadCounts = preprocessedStandardizedReadCounts;
            this.panelIntervalFractionalMedians = panelIntervalFractionalMedians;
            this.filterIntervals = filterIntervals;
            this.filterSamples = filterSamples;
        }
    }

    /**
     * Preprocess (i.e., filter, impute, and truncate) and standardize read counts from a panel of normals.
     * All inputs are assumed to be valid.
     * The dimensions of {@code readCounts} should be intervals x samples.
     * To reduce memory footprint, {@code readCounts} is modified in place when possible.
     * Filtering is performed by using boolean arrays to keep track of intervals and samples
     * that have been filtered at any step and masking {@code readCounts} with them appropriately.
     * If {@code intervalGCContent} is null, GC-bias correction will not be performed.
     */
    private static PreprocessedStandardizedResult preprocessAndStandardize(final RealMatrix readCounts,
                                                                           final double[] intervalGCContent,
                                                                           final double minimumIntervalMedianPercentile,
                                                                           final double maximumZerosInSamplePercentage,
                                                                           final double maximumZerosInIntervalPercentage,
                                                                           final double extremeSampleMedianPercentile,
                                                                           final double extremeOutlierTruncationPercentile) {
        final int numOriginalIntervals = readCounts.getRowDimension();
        final int numOriginalSamples = readCounts.getColumnDimension();

        final boolean[] filterIntervals = new boolean[numOriginalIntervals];
        final boolean[] filterSamples = new boolean[numOriginalSamples];

        logger.info("Transforming read counts to fractional coverage...");
        final double[] sampleSums = GATKProtectedMathUtils.columnSums(readCounts);
        readCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int intervalIndex, int sampleIndex, double value) {
                return value / sampleSums[sampleIndex];
            }
        });

        if (intervalGCContent != null) {
            logger.info("Performing explicit GC-bias correction...");
            GCCorrector.correctCoverage(readCounts, intervalGCContent);
        }

        final double[] originalIntervalMedians = MatrixSummaryUtils.getRowMedians(readCounts);
        final double minimumIntervalMedianThreshold = minimumIntervalMedianPercentile == 0. //Percentile requires arguments > 0
                ? Doubles.min(originalIntervalMedians)
                : new Percentile(minimumIntervalMedianPercentile).evaluate(originalIntervalMedians);
        IntStream.range(0, numOriginalIntervals)
                .filter(i -> originalIntervalMedians[i] < minimumIntervalMedianThreshold)
                .forEach(i -> filterIntervals[i] = true);
        logger.info(String.format("After filtering intervals with median (across samples) below %.2f percentile, %d out of %d intervals remain.",
                minimumIntervalMedianPercentile, countNumberPassingFilter(filterIntervals), numOriginalIntervals));

        logger.info("Dividing by interval medians...");
        for (int intervalIndex = 0; intervalIndex < numOriginalIntervals; intervalIndex++) {
            if (!filterIntervals[intervalIndex]) {
                for (int sampleIndex = 0; sampleIndex < numOriginalSamples; sampleIndex++) {
                    if (!filterSamples[sampleIndex]) {
                        final double value = readCounts.getEntry(intervalIndex, sampleIndex);
                        readCounts.setEntry(intervalIndex, sampleIndex, value / originalIntervalMedians[intervalIndex]);
                    }
                }
            }
        }
//        readCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
//            @Override
//            public double visit(int intervalIndex, int sampleIndex, double value) {
//                return value / originalIntervalMedians[intervalIndex];
//            }
//        });

        //remove samples with zeros
        for (int sampleIndex = 0; sampleIndex < numOriginalSamples; sampleIndex++) {
            int numZerosInSample = 0;
            for (int intervalIndex = 0; intervalIndex < numOriginalIntervals; intervalIndex++) {
                if (!filterIntervals[intervalIndex] && readCounts.getEntry(intervalIndex, sampleIndex) == 0.) {
                    numZerosInSample++;
                }
            }
            if (numZerosInSample > calculateMaximumZerosCount(numZerosInSample, maximumZerosInSamplePercentage)) {
                filterSamples[sampleIndex] = true;
            }
        }
        logger.info(String.format("After filtering samples with an amount of zero-coverage intervals above %.2f percent, %d out of %d samples remain.",
                maximumZerosInSamplePercentage, countNumberPassingFilter(filterSamples), numOriginalSamples));

        //remove intervals with zeros

        //remove extreme samples

        //impute zeros as interval median

        //truncate extreme counts

        logger.info("Dividing by sample medians and transforming to log2 space...");
        final double[] sampleMeans = MatrixSummaryUtils.getColumnMedians(readCounts);
        readCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int intervalIndex, int sampleIndex, double value) {
                return SVDDenoisingUtils.safeLog2(value / sampleMeans[sampleIndex]);
            }
        });

        logger.info("Subtracting median of sample medians...");
        final double[] sampleMedians = MatrixSummaryUtils.getColumnMedians(readCounts);
        final double medianOfSampleMedians = new Median().evaluate(sampleMedians);
        readCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int intervalIndex, int sampleIndex, double value) {
                return value - medianOfSampleMedians;
            }
        });

        logger.info("Panel read counts standardized.");

        //construct the final filtered results
        final int[] panelIntervalIndices = IntStream.range(0, numOriginalIntervals).filter(i -> !filterIntervals[i]).toArray();
        final int[] panelSamplesIndices = IntStream.range(0, numOriginalSamples).filter(i -> !filterSamples[i]).toArray();
        final RealMatrix preprocessedStandardizedReadCounts = readCounts.getSubMatrix(panelIntervalIndices, panelSamplesIndices);
        final double[] panelIntervalFractionalMedians = IntStream.range(0, numOriginalIntervals)
                .filter(i -> !filterIntervals[i]).mapToDouble(i -> originalIntervalMedians[i]).toArray();
        return new PreprocessedStandardizedResult(preprocessedStandardizedReadCounts, panelIntervalFractionalMedians, filterIntervals, filterSamples);
    }

    private static int countNumberPassingFilter(final boolean[] filter) {
        final int numPassingFilter = (int) IntStream.range(0, filter.length).filter(i -> !filter[i]).count();
        if (numPassingFilter == 0) {
            throw new UserException.BadInput("Filtering removed all samples or intervals.  Select less strict filtering criteria.");
        }
        return numPassingFilter;
    }

    private static int calculateMaximumZerosCount(final int numZeroCounts, final double percentage) {
        return (int) Math.ceil(numZeroCounts * percentage / 100.0);
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
        file.makeDoubleMatrix(ORIGINAL_READ_COUNTS_PATH, originalReadCounts.getData());
    }

    private void writeOriginalSampleFilenames(final List<String> originalSampleFilenames) {
        file.makeStringArray(ORIGINAL_SAMPLE_FILENAMES_PATH, originalSampleFilenames.toArray(new String[originalSampleFilenames.size()]));
    }

    private void writeOriginalIntervals(final List<SimpleInterval> originalIntervals) {
        writeIntervals(file, ORIGINAL_INTERVALS_PATH, originalIntervals);
    }

    private void writeOriginalIntervalGCContent(final double[] originalIntervalGCContent) {
        file.makeDoubleArray(ORIGINAL_INTERVAL_GC_CONTENT_PATH, originalIntervalGCContent);
    }

    private void writePanelSampleFilenames(final List<String> panelSampleFilenames) {
        file.makeStringArray(PANEL_SAMPLE_FILENAMES_PATH, panelSampleFilenames.toArray(new String[panelSampleFilenames.size()]));
    }

    private void writePanelIntervals(final List<SimpleInterval> panelIntervals) {
        writeIntervals(file, PANEL_INTERVALS_PATH, panelIntervals);
    }

    private void writePanelIntervalFractionalMedians(final double[] panelIntervalFractionalMedians) {
        file.makeDoubleArray(PANEL_INTERVAL_FRACTIONAL_MEDIANS_PATH, panelIntervalFractionalMedians);
    }

    private void writeSingularValues(final double[] singularValues) {
        file.makeDoubleArray(PANEL_SINGULAR_VALUES_PATH, singularValues);
    }

    private void writeLeftSingular(final double[][] leftSingular) {
        file.makeDoubleMatrix(PANEL_LEFT_SINGULAR_PATH, leftSingular);
    }

    private void writePanelLeftSingularPseudoinverse(final double[][] leftSingularPseudoinverse) {
        file.makeDoubleMatrix(PANEL_LEFT_SINGULAR_PSEUDOINVERSE_PATH, leftSingularPseudoinverse);
    }

    private static List<SimpleInterval> readIntervals(final HDF5File file, final String path) {
        final String[][] values = file.readStringMatrix(path, NUM_INTERVAL_COLUMNS_PATH);
        final List<SimpleInterval> result = new ArrayList<>(values.length);
        for (final String[] row : values) {
            result.add(new SimpleInterval(
                    row[IntervalColumn.CONTIG.index],
                    Integer.parseInt(row[IntervalColumn.START.index]),
                    Integer.parseInt(row[IntervalColumn.END.index])));
        }
        return result;
    }

    private static void writeIntervals(final HDF5File file, final String path, final List<SimpleInterval> intervals) {
        final String[][] values = new String[intervals.size()][NUM_INTERVAL_COLUMNS];
        for (int i = 0; i < intervals.size(); i++) {
            final SimpleInterval interval = intervals.get(i);
            values[i][IntervalColumn.CONTIG.index] = interval.getContig();
            values[i][IntervalColumn.START.index] = String.valueOf(interval.getStart());
            values[i][IntervalColumn.END.index] = String.valueOf(interval.getEnd());
        }
        file.makeStringMatrix(path, values, NUM_INTERVAL_COLUMNS_PATH);
    }
}
