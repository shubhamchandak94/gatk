package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.rsvd;

import com.google.common.annotations.VisibleForTesting;
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
 * Several attributes are stored transposed (in other words, the rows and columns are interchanged).
 * This dodges a very slow write time in HDF5, since we usually have many more rows (targets) than columns (samples),
 * and HDF5 writes matrices with few rows and many columns much faster than matrices with many rows and few columns.
 *
 * The following are stored as transposes:
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class HDF5RandomizedSVDReadCountPanelOfNormals implements SVDReadCountPanelOfNormals {
    private static final Logger logger = LogManager.getLogger(HDF5RandomizedSVDReadCountPanelOfNormals.class);

    private static final int CHUNK_DIVISOR = 16;    //limits number of intervals to 16777215
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
    private static final String ORIGINAL_READ_COUNTS_PATH = ORIGINAL_DATA_GROUP_NAME + "/read_counts_samples_by_intervals";
    private static final String ORIGINAL_SAMPLE_FILENAMES_PATH = ORIGINAL_DATA_GROUP_NAME + "/sample_filenames";
    private static final String ORIGINAL_INTERVALS_PATH = ORIGINAL_DATA_GROUP_NAME + "/intervals";
    private static final String ORIGINAL_INTERVAL_GC_CONTENT_PATH = ORIGINAL_DATA_GROUP_NAME + "/interval_gc_content";

    private static final String PANEL_GROUP_NAME = "/panel";
    private static final String PANEL_SAMPLE_FILENAMES_PATH = PANEL_GROUP_NAME + "/sample_filenames";
    private static final String PANEL_INTERVALS_PATH = PANEL_GROUP_NAME + "/intervals";
    private static final String PANEL_INTERVAL_FRACTIONAL_MEDIANS_PATH = PANEL_GROUP_NAME + "/interval_fractional_medians";
    private static final String PANEL_SINGULAR_VALUES_PATH = PANEL_GROUP_NAME + "/singular_values";
    private static final String PANEL_EIGENSAMPLE_VECTORS_PATH = PANEL_GROUP_NAME + "/transposed_eigensamples_samples_by_intervals";

    private final Lazy<List<Locatable>> originalIntervals;
    private final Lazy<List<Locatable>> panelIntervals;

    /**
     * <p>DEV NOTE:  If you are adding attributes that are neither RealMatrix nor a primitive,
     * you must follow the pattern in the constructor (i.e. the Lazy loading pattern).
     * Otherwise, some operations will hang.</p>
     */
    private HDF5RandomizedSVDReadCountPanelOfNormals(final HDF5File file) {
        Utils.nonNull(file);
        IOUtils.canReadFile(file.getFile());
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
        return HDF5ChunkedDoubleMatrixHelper.readMatrix(file, ORIGINAL_READ_COUNTS_PATH);
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
    public double[][] getEigensampleVectors() {
        return new Array2DRowRealMatrix(
                HDF5ChunkedDoubleMatrixHelper.readMatrix(file, PANEL_EIGENSAMPLE_VECTORS_PATH), false)
                .transpose().getData();
    }

    /**
     * Create an interface to an HDF5 file.  A version check is performed and a warning message logged if the
     * version number is not up to date.
     */
    public static HDF5RandomizedSVDReadCountPanelOfNormals read(final HDF5File file) {
        Utils.nonNull(file);
        IOUtils.canReadFile(file.getFile());
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

            logger.info(String.format("Writing original read counts (%d x %d)...",
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
                    SVDDenoisingUtils.preprocessAndStandardizePanel(originalReadCounts, intervalGCContent,
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

            final int numPanelSamples = preprocessedStandardizedResult.preprocessedStandardizedProfile.getRowDimension();
            final int numPanelIntervals = preprocessedStandardizedResult.preprocessedStandardizedProfile.getColumnDimension();
            final int numEigensamples = Math.min(numEigensamplesRequested, numPanelSamples);
            if (numEigensamples < numEigensamplesRequested) {
                logger.warn(String.format("%d eigensamples were requested but only %d are available in the panel of normals...",
                        numEigensamplesRequested, numEigensamples));
            }
            logger.info(String.format("Performing SVD (truncated at %d eigensamples) of standardized counts (transposed to %d x %d)...",
                    Math.min(numEigensamples, numPanelSamples), numPanelIntervals, numPanelSamples));
            final SingularValueDecomposition<RowMatrix, Matrix> svd = SparkConverter.convertRealMatrixToSparkRowMatrix(
                    ctx, preprocessedStandardizedResult.preprocessedStandardizedProfile.transpose(), NUM_SLICES_FOR_SPARK_MATRIX_CONVERSION)
                    .computeSVD(numEigensamples, true, EPSILON);
            final double[] singularValues = svd.s().toArray();    //should be in decreasing order (with corresponding matrices below)
            final double[][] eigensampleVectors = SparkConverter.convertSparkRowMatrixToRealMatrix(svd.U(), numPanelIntervals).getData();

            logger.info(String.format("Writing singular values (%d)...", singularValues.length));
            pon.writeSingularValues(singularValues);

            logger.info(String.format("Writing eigensample vectors (transposed to %d x %d)...", eigensampleVectors.length, eigensampleVectors[0].length));
            pon.writeEigensampleVectors(eigensampleVectors);
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
        HDF5ChunkedDoubleMatrixHelper.writeMatrix(file, ORIGINAL_READ_COUNTS_PATH, originalReadCounts.getData(), CHUNK_DIVISOR);
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

    private void writeEigensampleVectors(final double[][] eigensampleVectors) {
        HDF5ChunkedDoubleMatrixHelper.writeMatrix(file, PANEL_EIGENSAMPLE_VECTORS_PATH,
                new Array2DRowRealMatrix(eigensampleVectors, false).transpose().getData(),
                CHUNK_DIVISOR);
    }

    private static final class IntervalHelper {
        //writing intervals as a string matrix is expensive,
        //so we instead store a map from integer indices to contig strings and
        //store (index, start, end) in a double matrix

        private static final String INTERVAL_CONTIG_NAMES_SUB_PATH = "/indexed_contig_names";
        private static final String INTERVAL_MATRIX_SUB_PATH = "/transposed_index_start_end";

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
            final String[] contigNames = file.readStringArray(path + INTERVAL_CONTIG_NAMES_SUB_PATH);
            final double[][] matrix = file.readDoubleMatrix(path + INTERVAL_MATRIX_SUB_PATH);
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
            file.makeStringArray(path + INTERVAL_CONTIG_NAMES_SUB_PATH, contigNames);
            final Map<String, Double> contigNamesToIndexMap = IntStream.range(0, contigNames.length).boxed()
                    .collect(Collectors.toMap(i -> contigNames[i], i -> (double) i));
            final double[][] matrix = new double[NUM_INTERVAL_FIELDS][intervals.size()];
            for (int i = 0; i < intervals.size(); i++) {
                final Locatable interval = intervals.get(i);
                matrix[IntervalField.CONTIG_INDEX.index][i] = contigNamesToIndexMap.get(interval.getContig());
                matrix[IntervalField.START.index][i] = interval.getStart();
                matrix[IntervalField.END.index][i] = interval.getEnd();
            }
            file.makeDoubleMatrix(path + INTERVAL_MATRIX_SUB_PATH, matrix);
        }
    }

    @VisibleForTesting
    static final class HDF5ChunkedDoubleMatrixHelper {

        private static final String NUMBER_OF_ROWS_SUB_PATH = "/num_rows";
        private static final String NUMBER_OF_COLUMNS_SUB_PATH = "/num_columns";
        private static final String NUMBER_OF_CHUNKS_SUB_PATH = "/num_chunks";
        private static final String CHUNK_INDEX_PATH_SUFFIX = "/chunk_";

        private static final int MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX = Integer.MAX_VALUE / Byte.SIZE;

        @VisibleForTesting
        static double[][] readMatrix(final HDF5File file,
                                     final String path) {
            Utils.nonNull(file);
            IOUtils.canReadFile(file.getFile());
            Utils.nonNull(path);

            final String numRowsPath = path + NUMBER_OF_ROWS_SUB_PATH;
            final String numColumnsPath = path + NUMBER_OF_COLUMNS_SUB_PATH;
            final String numChunksPath = path + NUMBER_OF_CHUNKS_SUB_PATH;
            Utils.validateArg(file.isPresent(numRowsPath) && file.isPresent(numColumnsPath) && file.isPresent(numChunksPath),
                    String.format("HDF5 file %s does not contain a chunked matrix in path %s.", file.getFile().getAbsolutePath(), path));

            final int numRows = (int) file.readDouble(numRowsPath);
            final int numColumns = (int) file.readDouble(numColumnsPath);
            final int numChunks = (int) file.readDouble(numChunksPath);

            final double[][] fullMatrix = new double[numRows][numColumns];
            int numRowsRead = 0;
            for (int chunkIndex = 0; chunkIndex < numChunks; chunkIndex++) {
                final double[][] matrixChunk = file.readDoubleMatrix(path + CHUNK_INDEX_PATH_SUFFIX + chunkIndex);
                if (numRowsRead + matrixChunk.length > numRows) {
                    throw new UserException.BadInput("Matrix chunk contains too many rows.");
                }
                if (matrixChunk[0].length != numColumns) {
                    throw new UserException.BadInput("Matrix chunk does not contain expected number of columns.");
                }
                System.arraycopy(matrixChunk, 0, fullMatrix, numRowsRead, matrixChunk.length);
                numRowsRead += matrixChunk.length;
            }
            if (numRowsRead != numRows) {
                throw new UserException.BadInput("Matrix chunks do not contain expected total number of rows.");
            }
            return fullMatrix;
        }

        /**
         * @param chunkDivisor  The maximum number of values in each chunk
         *                      is given by {@code MAX_NUM_VALUES_PER_HDF5_MATRIX} / {@code chunkDivisor},
         *                      so increasing this number will reduce heap usage when writing chunks,
         *                      which requires subarrays to be copied.  However, since a single row is not allowed
         *                      to be split across multiple chunks, the number of columns must be less
         *                      than the maximum number of values in each chunk.  For example,
         *                      {@code chunkDivisor} = 8 allows for 16777215 columns.
         */
        @VisibleForTesting
        static void writeMatrix(final HDF5File file,
                                final String path,
                                final double[][] matrix,
                                final int chunkDivisor) {
            Utils.nonNull(file);
            IOUtils.canReadFile(file.getFile());
            Utils.nonNull(path);
            Utils.nonNull(matrix);
            Utils.validateArg(chunkDivisor > 0, "Chunk divisor must be positive.");
            final int maxNumValuesPerChunk = MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX / chunkDivisor;
            final long numRows = matrix.length;
            Utils.validateArg(numRows > 0, "Matrix must contain at least one row.");
            final long numColumns = matrix[0].length;
            Utils.validateArg(numColumns > 0, "Matrix must contain at least one column.");
            Utils.validateArg(numColumns <= maxNumValuesPerChunk,
                    String.format("Number of columns (%d) exceeds the maximum number of values allowed per chunk (%d).",
                            numColumns, maxNumValuesPerChunk));

            final int numRowsPerFilledChunk = (int) (maxNumValuesPerChunk / numColumns);
            final int numFilledChunks = numRowsPerFilledChunk == 0 ? 0 : (int) numRows / numRowsPerFilledChunk;
            final boolean needPartialChunk = numFilledChunks == 0 || numRows % numRowsPerFilledChunk != 0;

            logger.debug("Number of values in matrix / maximum number allowed for HDF5 matrix: " + (double) numRows * numColumns / MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX);
            logger.debug("Maximum number of values per chunk: " + maxNumValuesPerChunk);
            logger.debug("Number of filled chunks: " + numFilledChunks);
            logger.debug("Number of rows per filled chunk: " + numRowsPerFilledChunk);
            logger.debug("Partial chunk needed: " + needPartialChunk);

            final String numRowsPath = path + NUMBER_OF_ROWS_SUB_PATH;
            final String numColumnsPath = path + NUMBER_OF_COLUMNS_SUB_PATH;
            final String numChunksPath = path + NUMBER_OF_CHUNKS_SUB_PATH;
            file.makeDouble(numRowsPath, numRows);
            file.makeDouble(numColumnsPath, numColumns);
            file.makeDouble(numChunksPath, needPartialChunk ? numFilledChunks + 1 : numFilledChunks);

            //TODO we could add makeDoubleMatrix(path, matrix, startRow, endRow, startCol, endCol) method to avoid copying
            int numRowsWritten = 0;
            for (int chunkIndex = 0; chunkIndex < numFilledChunks; chunkIndex++) {
                final double[][] matrixChunk = new double[numRowsPerFilledChunk][(int) numColumns];
                System.arraycopy(matrix, numRowsWritten, matrixChunk, 0, numRowsPerFilledChunk);
                file.makeDoubleMatrix(path + CHUNK_INDEX_PATH_SUFFIX + chunkIndex, matrixChunk);    //write filled chunks
                numRowsWritten += numRowsPerFilledChunk;
            }
            if (needPartialChunk) {
                final int numRowsPartialChunk = (int) numRows - numRowsWritten;
                logger.debug("Number of rows in partial chunk: " + numRowsPartialChunk);
                final double[][] matrixChunk = new double[numRowsPartialChunk][(int) numColumns];
                System.arraycopy(matrix, numRowsWritten, matrixChunk, 0, numRowsPartialChunk);
                file.makeDoubleMatrix(path + CHUNK_INDEX_PATH_SUFFIX + numFilledChunks, matrixChunk);    //write final partially filled chunk
            }
        }
    }
}
