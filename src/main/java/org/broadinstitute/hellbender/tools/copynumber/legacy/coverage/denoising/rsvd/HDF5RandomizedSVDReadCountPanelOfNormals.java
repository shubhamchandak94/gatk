package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.rsvd;

import htsjdk.samtools.util.Lazy;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

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
     */
    private static final double CURRENT_PON_VERSION = 7.0;

    private static final String VERSION_PATH = "/version/values";
    private static final String COMMAND_LINE_PATH = "/command_line";
    private static final String ORIGINAL_DATA_GROUP_NAME = "/original_data";
    private static final String ORIGINAL_READ_COUNTS_PATH = ORIGINAL_DATA_GROUP_NAME + "/read_counts";
    private static final String ORIGINAL_SAMPLE_FILENAMES_PATH = ORIGINAL_DATA_GROUP_NAME + "/sample_filenames";
    private static final String ORIGINAL_INTERVALS_PATH = ORIGINAL_DATA_GROUP_NAME + "/intervals";

    private static final String PANEL_GROUP_NAME = "/panel";
    private static final String PANEL_SAMPLE_FILENAMES_PATH = PANEL_GROUP_NAME + "/sample_filenames";
    private static final String PANEL_INTERVALS_PATH = PANEL_GROUP_NAME + "/intervals";
    private static final String PANEL_INTERVAL_FRACTIONAL_MEDIANS_PATH = PANEL_GROUP_NAME + "/interval_fractional_medians";
    private static final String PANEL_SINGULAR_VALUES_PATH = PANEL_GROUP_NAME + "/singular_values";
    private static final String PANEL_RIGHT_SINGULAR_PATH = PANEL_GROUP_NAME + "/right_singular";
    private static final String PANEL_RIGHT_SINGULAR_PSEUDOINVERSE_PATH = PANEL_GROUP_NAME + "/right_singular_pseudoinverse";

    private static final String NUM_COLUMNS_PATH = "/num_columns";

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

    private Lazy<List<SimpleInterval>> originalIntervals;
    private Lazy<List<SimpleInterval>> panelIntervals;

    /**
     * <p>DEV NOTE:  If you are adding attributes that are neither RealMatrix nor a primitive,
     * you must follow the pattern in the constructor (i.e. the Lazy loading pattern).
     * Otherwise, some operations will hang.</p>
     */
    private HDF5RandomizedSVDReadCountPanelOfNormals(final HDF5File file) {
        Utils.nonNull(file);
        this.file = file;
        originalIntervals  = new Lazy<>(() -> readIntervals(file, ORIGINAL_INTERVALS_PATH));
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
    public List<SimpleInterval> getPanelIntervals() {
        return panelIntervals.get();
    }

    @Override
    public double[] getPanelIntervalFractionalMedians() {
        return file.readDoubleArray(PANEL_INTERVAL_FRACTIONAL_MEDIANS_PATH);
    }

    @Override
    public double[] getSingularValues() {
        return file.readDoubleArray(PANEL_SINGULAR_VALUES_PATH)
    }

    @Override
    public double[][] getRightSingular() {
        return file.readDoubleMatrix(PANEL_RIGHT_SINGULAR_PATH);
    }

    @Override
    public double[][] getRightSingularPseudoinverse() {
        return file.readDoubleMatrix(PANEL_RIGHT_SINGULAR_PSEUDOINVERSE_PATH);
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
     * Create the panel of normals and write it to an HDF5 file.
     */
    public static void create(final File outFile,
                              final String commandLine,
                              final RealMatrix originalReadCounts,
                              final List<String> originalSampleFilenames,
                              final List<SimpleInterval> originalIntervals) {
        Utils.nonNull(outFile);
        //TODO validate all inputs

        try (final HDF5File file = new HDF5File(outFile, HDF5File.OpenMode.CREATE)) {
            logger.info("Creating " + outFile.getAbsolutePath() + "...");
            final HDF5RandomizedSVDReadCountPanelOfNormals pon = new HDF5RandomizedSVDReadCountPanelOfNormals(file);

            logger.info(String.format("Writing version number (%f)...", CURRENT_PON_VERSION));
            pon.setVersion(CURRENT_PON_VERSION);

            logger.info("Writing command line...");
            pon.setCommandLine(commandLine);

            logger.info("Writing original read counts ( x )...");
            pon.setOriginalReadCountsPath(originalReadCounts);

            logger.info(String.format("Writing original sample filenames (%d)...", originalSampleFilenames.size()));
            pon.setOriginalSampleFilenames(originalSampleFilenames);

            logger.info(String.format("Writing original intervals (%d)...", originalIntervals.size()));
            pon.setOriginalIntervals(originalIntervals);

            logger.info(String.format("Writing panel sample filenames (%d)...", panelSampleFilenames.size()));
            pon.setPanelSampleFilenames(panelSampleFilenames);

            logger.info(String.format("Writing panel intervals (%d)...", panelIntervals.size()));
            pon.setPanelIntervals(panelIntervals);

            logger.info(String.format("Writing panel interval fractional medians (%d)...", panelIntervalFractionalMedians.length));
            pon.setPanelIntervalFractionalMedians(panelIntervalFractionalMedians);

            logger.info(String.format("Writing singular values (%d)...", singularValues.length));
            pon.setSingularValues(singularValues);

            logger.info(String.format("Writing right-singular matrix (%d x %d)...", rightSingular.length, rightSingular[0].length));
            pon.setRightSingular(rightSingular);

            logger.info(String.format("Writing right-singular pseudoinverse (%d x %d)...", rightSingularPseudoinverse.length, rightSingularPseudoinverse[0].length));
            pon.setPanelRightSingularPseudoinverse(rightSingularPseudoinverse);
        }
        logger.info(String.format("Read-count panel of normals written to %s.", outFile));
    }

    //PRIVATE SETTERS (write values to HDF5 file and set internally held, lazily loaded fields to them)
    //these are private to prevent fields from being written individually, which could leave the file in a bad state

    private void setVersion(final double version) {
        file.makeDouble(VERSION_PATH, version);
    }

    private void setCommandLine(final String commandLine) {
        file.makeStringArray(COMMAND_LINE_PATH, commandLine);
    }

    private void setOriginalReadCountsPath(final RealMatrix originalReadCounts) {
        file.makeDoubleMatrix(ORIGINAL_READ_COUNTS_PATH, originalReadCounts.getData());
    }

    private void setOriginalSampleFilenames(final List<String> originalSampleFilenames) {
        file.makeStringArray(ORIGINAL_SAMPLE_FILENAMES_PATH, originalSampleFilenames.toArray(new String[originalSampleFilenames.size()]));
    }

    private void setOriginalIntervals(final List<SimpleInterval> originalIntervals) {
        writeIntervals(file, ORIGINAL_INTERVALS_PATH, originalIntervals);
        this.originalIntervals = new Lazy<>(() -> readIntervals(file, ORIGINAL_INTERVALS_PATH));
    }

    private void setPanelSampleFilenames(final List<String> panelSampleFilenames) {
        file.makeStringArray(PANEL_SAMPLE_FILENAMES_PATH, panelSampleFilenames.toArray(new String[panelSampleFilenames.size()]));
    }

    private void setPanelIntervals(final List<SimpleInterval> panelIntervals) {
        writeIntervals(file, PANEL_INTERVALS_PATH, panelIntervals);
        this.panelIntervals = new Lazy<>(() -> readIntervals(file, PANEL_INTERVALS_PATH));
    }

    private void setPanelIntervalFractionalMedians(final double[] panelIntervalFractionalMedians) {
        file.makeDoubleArray(PANEL_INTERVAL_FRACTIONAL_MEDIANS_PATH, panelIntervalFractionalMedians);
    }

    private void setSingularValues(final double[] singularValues) {
        file.makeDoubleArray(PANEL_SINGULAR_VALUES_PATH, singularValues);
    }

    private void setRightSingular(final double[][] rightSingular) {
        file.makeDoubleMatrix(PANEL_RIGHT_SINGULAR_PATH, rightSingular);
    }

    private void setPanelRightSingularPseudoinverse(final double[][] rightSingularPseudoinverse) {
        file.makeDoubleMatrix(PANEL_RIGHT_SINGULAR_PSEUDOINVERSE_PATH, rightSingularPseudoinverse);
    }

    //HDF5 array/matrix read/write methods

    private static List<SimpleInterval> readIntervals(final HDF5File file, final String path) {
        final String[][] values = readStringMatrix(file, path);
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
        writeStringMatrix(file, path, values);
    }

    private static String[][] readStringMatrix(final HDF5File file, final String path) {
        return file.readStringMatrix(path, path + NUM_COLUMNS_PATH);
    }

    private static void writeStringMatrix(final HDF5File file, final String path, final String[][] values) {
        file.makeStringMatrix(path, values, path + NUM_COLUMNS_PATH);
    }
}
