package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.rsvd;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.legacy.CreateReadCountPanelOfNormals;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.gcbias.GCCorrector;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MatrixSummaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Utility class for package-private methods for performing SVD-based denoising and related operations.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SVDDenoisingUtils {
    private static final Logger logger = LogManager.getLogger(SVDDenoisingUtils.class);

    private static final double EPSILON = 1E-9;
    private static final double INV_LN2 = GATKProtectedMathUtils.INV_LOG_2;
    private static final double LN2_EPSILON = Math.log(EPSILON) * INV_LN2;

    private SVDDenoisingUtils() {}

    //TODO remove this method once ReadCountCollection is refactored to only store single sample, non-negative integer counts
    public static void validateReadCounts(final ReadCountCollection readCountCollection) {
        Utils.nonNull(readCountCollection);
        if (readCountCollection.columnNames().size() != 1) {
            throw new UserException.BadInput("Read-count file must contain counts for only a single sample.");
        }
        if (readCountCollection.targets().isEmpty()) {
            throw new UserException.BadInput("Read-count file must contain counts for at least one genomic interval.");
        }
        final double[] readCounts = readCountCollection.counts().getColumn(0);
        if (!IntStream.range(0, readCounts.length).allMatch(i -> (readCounts[i] >= 0) && ((int) readCounts[i] == readCounts[i]))) {
            throw new UserException.BadInput("Read-count file must contain non-negative integer counts.");
        }
    }

    static final class PreprocessedStandardizedResult {
        final RealMatrix preprocessedStandardizedReadCounts;
        final double[] panelIntervalFractionalMedians;
        final boolean[] filterIntervals;
        final boolean[] filterSamples;

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
    static PreprocessedStandardizedResult preprocessAndStandardize(final RealMatrix readCounts,
                                                                   final double[] intervalGCContent,
                                                                   final double minimumIntervalMedianPercentile,
                                                                   final double maximumZerosInSamplePercentage,
                                                                   final double maximumZerosInIntervalPercentage,
                                                                   final double extremeSampleMedianPercentile,
                                                                   final double extremeOutlierTruncationPercentile) {
        logger.info("Transforming read counts to fractional coverage...");
        transformToFractionalCoverage(readCounts);

        if (intervalGCContent != null) {
            logger.info("Performing explicit GC-bias correction...");
            GCCorrector.correctCoverage(readCounts, intervalGCContent);
        }

        final PreprocessedStandardizedResult preprocessedStandardizedResult = preprocessPanel(readCounts,
                minimumIntervalMedianPercentile, maximumZerosInSamplePercentage, maximumZerosInIntervalPercentage,
                extremeSampleMedianPercentile, extremeOutlierTruncationPercentile);

        logger.info("Dividing by sample medians and transforming to log2 space...");
        divideBySampleMedianAndTransformToLog2(preprocessedStandardizedResult.preprocessedStandardizedReadCounts);

        logger.info("Subtracting median of sample medians...");
        final double[] sampleLog2Medians = MatrixSummaryUtils.getColumnMedians(preprocessedStandardizedResult.preprocessedStandardizedReadCounts);
        final double medianOfSampleMedians = new Median().evaluate(sampleLog2Medians);
        preprocessedStandardizedResult.preprocessedStandardizedReadCounts
                .walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int intervalIndex, int sampleIndex, double value) {
                return value - medianOfSampleMedians;
            }
        });

        logger.info("Panel read counts standardized.");

        return preprocessedStandardizedResult;
    }

    /**
     * Perform SVD-based denoising of integer read counts for a single sample using a panel of normals.
     * Only the eigensamples (which are sorted by singular value in decreasing order) specified by
     * {@code numEigensamples} are used to denoise.
     */
    static SVDDenoisedCopyRatioResult denoise(final SVDReadCountPanelOfNormals panelOfNormals,
                                              final ReadCountCollection readCounts,
                                              final int numEigensamples) {
        Utils.nonNull(panelOfNormals);
        validateReadCounts(readCounts);
        ParamUtils.isPositive(numEigensamples, "Number of eigensamples to use for denoising must be positive.");
        Utils.validateArg(numEigensamples <= panelOfNormals.getNumEigensamples(),
                "Number of eigensamples to use for denoising is greater than the number available in the panel of normals.");

        logger.info("Validating sample intervals against original intervals used to build panel of normals...");
        Utils.validateArg(panelOfNormals.getOriginalIntervals().equals(readCounts.targets().stream().map(Target::getInterval).collect(Collectors.toList())),
                "Sample intervals must be identical to the original intervals used to build the panel of normals.");

        logger.info("Standardizing sample read counts...");
        final RealMatrix standardizedCounts = preprocessAndStandardizeSample(panelOfNormals, readCounts.counts());

        logger.info(String.format("Using %d out of %d eigensamples to denoise...", numEigensamples, panelOfNormals.getNumEigensamples()));

        logger.info("Subtracting projection onto space spanned by eigensamples...");
        final RealMatrix denoisedCounts = subtractProjection(standardizedCounts, panelOfNormals.getLeftSingular(), numEigensamples);

        logger.info("Sample denoised.");

        //construct the result
        //TODO clean this up once Targets are removed
        final Set<SimpleInterval> panelIntervals = new HashSet<>(panelOfNormals.getPanelIntervals());
        final List<Target> targets = readCounts.targets().stream().filter(t -> panelIntervals.contains(t.getInterval())).collect(Collectors.toList());
        final ReadCountCollection standardizedProfile = new ReadCountCollection(targets, readCounts.columnNames(), standardizedCounts);
        final ReadCountCollection denoisedProfile = new ReadCountCollection(targets, readCounts.columnNames(), denoisedCounts);

        return new SVDDenoisedCopyRatioResult(standardizedProfile, denoisedProfile);
    }

    /**
     * TODO
     */
    private static PreprocessedStandardizedResult preprocessPanel(final RealMatrix transformedReadCounts,
                                                                  final double minimumIntervalMedianPercentile,
                                                                  final double maximumZerosInSamplePercentage,
                                                                  final double maximumZerosInIntervalPercentage,
                                                                  final double extremeSampleMedianPercentile,
                                                                  final double extremeOutlierTruncationPercentile) {
        final int numOriginalIntervals = transformedReadCounts.getRowDimension();
        final int numOriginalSamples = transformedReadCounts.getColumnDimension();

        final boolean[] filterIntervals = new boolean[numOriginalIntervals];
        final boolean[] filterSamples = new boolean[numOriginalSamples];

        final double[] originalIntervalMedians = MatrixSummaryUtils.getRowMedians(transformedReadCounts);

        //filter intervals by fractional median
        if (minimumIntervalMedianPercentile == 0.) {
            logger.info(String.format("A value of 0 was provided for argument %s, so the corresponding filtering step will be skipped...",
                    CreateReadCountPanelOfNormals.MINIMUM_INTERVAL_MEDIAN_PERCENTILE_LONG_NAME));
        } else {
            logger.info(String.format("Filtering intervals with median (across samples) below the %.2f percentile...", minimumIntervalMedianPercentile));
            //calculate percentile
            final double minimumIntervalMedianThreshold = new Percentile(minimumIntervalMedianPercentile).evaluate(originalIntervalMedians);
            //filter intervals
            IntStream.range(0, numOriginalIntervals)
                    .filter(intervalIndex -> originalIntervalMedians[intervalIndex] < minimumIntervalMedianThreshold)
                    .forEach(intervalIndex -> filterIntervals[intervalIndex] = true);
            logger.info(String.format("After filtering, %d out of %d intervals remain...", countNumberPassingFilter(filterIntervals), numOriginalIntervals));
        }

        logger.info("Dividing by interval medians...");
        IntStream.range(0, numOriginalIntervals)
                .filter(intervalIndex -> !filterIntervals[intervalIndex])
                .forEach(intervalIndex -> IntStream.range(0, numOriginalSamples).filter(sampleIndex -> !filterSamples[sampleIndex]).forEach(sampleIndex -> {
                    final double value = transformedReadCounts.getEntry(intervalIndex, sampleIndex);
                    transformedReadCounts.setEntry(intervalIndex, sampleIndex, value / originalIntervalMedians[intervalIndex]);
                }));

        //filter samples by percentage of zero-coverage intervals not already filtered
        if (maximumZerosInSamplePercentage == 100.) {
            logger.info(String.format("A value of 100 was provided for argument %s, so the corresponding filtering step will be skipped...",
                    CreateReadCountPanelOfNormals.MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE_LONG_NAME));
        } else {
            logger.info(String.format("Filtering samples with a fraction of zero-coverage intervals above %.2f percent...", maximumZerosInSamplePercentage));
            IntStream.range(0, numOriginalSamples)
                    .filter(sampleIndex -> !filterSamples[sampleIndex])
                    .forEach(sampleIndex -> {
                        final int numZerosInSample = (int) IntStream.range(0, numOriginalIntervals)
                                .filter(intervalIndex -> !filterIntervals[intervalIndex] && transformedReadCounts.getEntry(intervalIndex, sampleIndex) == 0.)
                                .count();
                        if (numZerosInSample > calculateMaximumZerosCount(numZerosInSample, maximumZerosInSamplePercentage)) {
                            filterSamples[sampleIndex] = true;
                        }
                    });
            logger.info(String.format("After filtering, %d out of %d samples remain...", countNumberPassingFilter(filterSamples), numOriginalSamples));
        }

        //filter intervals by percentage of zero-coverage samples not already filtered
        if (maximumZerosInIntervalPercentage == 100.) {
            logger.info(String.format("A value of 100 was provided for argument %s, so the corresponding filtering step will be skipped...",
                    CreateReadCountPanelOfNormals.MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE_LONG_NAME));
        } else {
            logger.info(String.format("Filtering intervals with a fraction of zero-coverage samples above %.2f percent...", maximumZerosInIntervalPercentage));
            IntStream.range(0, numOriginalIntervals)
                    .filter(intervalIndex -> !filterIntervals[intervalIndex])
                    .forEach(intervalIndex -> {
                        final int numZerosInInterval = (int) IntStream.range(0, numOriginalSamples)
                                .filter(sampleIndex -> !filterSamples[sampleIndex] && transformedReadCounts.getEntry(sampleIndex, sampleIndex) == 0.)
                                .count();
                        if (numZerosInInterval > calculateMaximumZerosCount(numZerosInInterval, maximumZerosInIntervalPercentage)) {
                            filterIntervals[intervalIndex] = true;
                        }
                    });
            logger.info(String.format("After filtering, %d out of %d intervals remain...", countNumberPassingFilter(filterIntervals), numOriginalIntervals));
        }

        //filter samples with extreme medians
        if (extremeSampleMedianPercentile == 0.) {
            logger.info(String.format("A value of 0 was provided for argument %s, so the corresponding filtering step will be skipped...",
                    CreateReadCountPanelOfNormals.EXTREME_SAMPLE_MEDIAN_PERCENTILE_LONG_NAME));
        } else {
            logger.info(String.format("Filtering samples with a median (across intervals) below the %.2f percentile or above the %.2f percentile...",
                    extremeSampleMedianPercentile, 100. - extremeSampleMedianPercentile));
            //calculate the medians for all samples (which, although unnecessary, makes bookkeeping easier) across intervals not already filtered
            final double[] sampleMedians = IntStream.range(0, numOriginalSamples)
                    .mapToDouble(sampleIndex -> new Median().evaluate(IntStream.range(0, numOriginalIntervals)
                            .filter(intervalIndex -> !filterIntervals[intervalIndex])
                            .mapToDouble(intervalIndex -> transformedReadCounts.getEntry(intervalIndex, sampleIndex))
                            .toArray()))
                    .toArray();
            //calculate percentiles
            final double minimumSampleMedianThreshold = new Percentile(extremeSampleMedianPercentile).evaluate(sampleMedians);
            final double maximumSampleMedianThreshold = new Percentile(100. - extremeSampleMedianPercentile).evaluate(sampleMedians);
            //filter samples
            IntStream.range(0, numOriginalSamples)
                    .filter(sampleIndex -> sampleMedians[sampleIndex] < minimumSampleMedianThreshold || sampleMedians[sampleIndex] > maximumSampleMedianThreshold)
                    .forEach(sampleIndex -> filterSamples[sampleIndex] = true);
            logger.info(String.format("After filtering, %d out of %d samples remain...", countNumberPassingFilter(filterSamples), numOriginalSamples));
        }

        //construct the filtered results as a new matrix, which will be modified in place from this point on
        final int[] panelIntervalIndices = IntStream.range(0, numOriginalIntervals).filter(intervalIndex -> !filterIntervals[intervalIndex]).toArray();
        final int[] panelSampleIndices = IntStream.range(0, numOriginalSamples).filter(sampleIndex -> !filterSamples[sampleIndex]).toArray();
        final RealMatrix preprocessedStandardizedReadCounts = transformedReadCounts.getSubMatrix(panelIntervalIndices, panelSampleIndices);
        final double[] panelIntervalFractionalMedians = IntStream.range(0, numOriginalIntervals)
                .filter(intervalIndex -> !filterIntervals[intervalIndex])
                .mapToDouble(intervalIndex -> originalIntervalMedians[intervalIndex]).toArray();

        //impute zeros as median of non-zero values in interval
        //TODO make this optional
        final int numPanelIntervals = panelIntervalIndices.length;
        final double[] intervalNonZeroMedians = IntStream.range(0, numPanelIntervals)
                .mapToObj(intervalIndex -> Arrays.stream(preprocessedStandardizedReadCounts.getRow(intervalIndex)).filter(value -> value > 0.).toArray())
                .mapToDouble(nonZeroValues -> new Median().evaluate(nonZeroValues))
                .toArray();
        final int[] numImputed = {0};  //needs to be effectively final to be used inside visitor
        preprocessedStandardizedReadCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int intervalIndex, int sampleIndex, double value) {
                if (value == 0.) {
                    numImputed[0]++;
                    return intervalNonZeroMedians[intervalIndex];
                }
                return value;
            }
        });
        logger.info(String.format("%d zero-coverage values were imputed to the median of the non-zero values in the corresponding interval...",
                numImputed[0]));

        //truncate extreme values to the corresponding percentile
        if (extremeOutlierTruncationPercentile == 0.) {
            logger.info(String.format("A value of 0 was provided for argument %s, so the corresponding truncation step will be skipped...",
                    CreateReadCountPanelOfNormals.EXTREME_OUTLIER_TRUNCATION_PERCENTILE_LONG_NAME));
        } else {
            final double[] values = Doubles.concat(preprocessedStandardizedReadCounts.getData());
            final double minimumOutlierTruncationThreshold = new Percentile(extremeOutlierTruncationPercentile).evaluate(values);
            final double maximumOutlierTruncationThreshold = new Percentile(100. - extremeOutlierTruncationPercentile).evaluate(values);
            final int[] numTruncated = {0};  //needs to be effectively final to be used inside visitor
            preprocessedStandardizedReadCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
                @Override
                public double visit(int intervalIndex, int sampleIndex, double value) {
                    if (value < minimumOutlierTruncationThreshold) {
                        numTruncated[0]++;
                        return minimumOutlierTruncationThreshold;
                    }
                    if (value > maximumOutlierTruncationThreshold) {
                        numTruncated[0]++;
                        return maximumOutlierTruncationThreshold;
                    }
                    return value;
                }
            });
            logger.info(String.format("%d values below the %.2f percentile or above the %.2f percentile were truncated to the corresponding value...",
                    numTruncated[0], extremeOutlierTruncationPercentile, 100. - extremeOutlierTruncationPercentile));
        }
        return new PreprocessedStandardizedResult(transformedReadCounts, panelIntervalFractionalMedians,
                filterIntervals, filterSamples);
    }

    /**
     * Standardize read counts for a single sample, using interval fractional medians from a panel of normals.
     * The original {@code readCounts} has dimensions intervals x 1 and is not modified.
     */
    private static RealMatrix preprocessAndStandardizeSample(final SVDReadCountPanelOfNormals panelOfNormals,
                                                             final RealMatrix readCounts) {
        RealMatrix result = readCounts.copy();

        logger.info("Transforming read counts to fractional coverage...");
        transformToFractionalCoverage(result);

        if (panelOfNormals.getOriginalIntervalGCContent() != null) {
            logger.info("GC-content annotations found in the panel of normals; performing explicit GC-bias correction...");
            GCCorrector.correctCoverage(result, panelOfNormals.getOriginalIntervalGCContent());
        }

        logger.info("Subsetting sample intervals to post-filter panel intervals...");
        final Set<SimpleInterval> panelIntervals = new HashSet<>(panelOfNormals.getPanelIntervals());
        final int[] subsetIntervalIndices = IntStream.range(0, panelOfNormals.getOriginalIntervals().size())
                .filter(i -> panelIntervals.contains(panelOfNormals.getOriginalIntervals().get(i)))
                .toArray();
        result = result.getSubMatrix(subsetIntervalIndices, new int[]{0});

        logger.info("Dividing by interval medians from the panel of normals...");
        final double[] intervalMedians = panelOfNormals.getPanelIntervalFractionalMedians();
        result.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int intervalIndex, int sampleIndex, double value) {
                return value / intervalMedians[intervalIndex];
            }
        });

        logger.info("Dividing by sample median and transforming to log2 space...");
        divideBySampleMedianAndTransformToLog2(result);

        logger.info("Subtracting sample median...");
        final double[] sampleLog2Medians = MatrixSummaryUtils.getColumnMedians(result);
        result.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int intervalIndex, int sampleIndex, double value) {
                return value - sampleLog2Medians[sampleIndex];
            }
        });

        logger.info("Sample read counts standardized.");

        return result;
    }

    /**
     * Given standardized read counts specified by a column vector S (dimensions {@code M x 1})
     * and left-singular vectors U (dimensions {@code M x K}),
     * returns s - U<sub>k</sub> U<sub>k</sub><sup>T</sup> s,
     * where U<sub>k</sub> contains the first {@code numEigensamples}.
     */
    private static RealMatrix subtractProjection(final RealMatrix standardizedProfile,
                                                 final double[][] leftSingular,
                                                 final int numEigensamples) {
        final int numIntervals = leftSingular.length;
        final int numAllEigensamples = leftSingular[0].length;

        logger.info("Distributing the standardized read counts...");

        logger.info("Composing left-singular matrix for the requested number of eigensamples and transposing them...");
        final RealMatrix leftSingularTruncatedMatrix = numEigensamples == numAllEigensamples
                ? new Array2DRowRealMatrix(leftSingular)
                : new Array2DRowRealMatrix(leftSingular).getSubMatrix(0, numIntervals - 1, 0, numEigensamples - 1);

        logger.info("Computing projection of transpose...");
        final RealMatrix projectionTranspose = standardizedProfile.transpose()
                .multiply(leftSingularTruncatedMatrix)
                .multiply(leftSingularTruncatedMatrix.transpose());

        logger.info("Subtracting projection...");
        return standardizedProfile.subtract(projectionTranspose.transpose());
    }

    private static int countNumberPassingFilter(final boolean[] filter) {
        final int numPassingFilter = (int) IntStream.range(0, filter.length).filter(i -> !filter[i]).count();
        if (numPassingFilter == 0) {
            throw new UserException.BadInput("Filtering removed all samples or intervals.  Select less strict filtering criteria.");
        }
        return numPassingFilter;
    }

    private static void transformToFractionalCoverage(final RealMatrix m) {
        final double[] sampleSums = GATKProtectedMathUtils.columnSums(m);
        m.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int intervalIndex, int sampleIndex, double value) {
                return value / sampleSums[sampleIndex];
            }
        });
    }

    private static void divideBySampleMedianAndTransformToLog2(final RealMatrix m) {
        final double[] sampleMedians = MatrixSummaryUtils.getColumnMedians(m);
        m.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int intervalIndex, int sampleIndex, double value) {
                return safeLog2(value / sampleMedians[sampleIndex]);
            }
        });
    }

    private static int calculateMaximumZerosCount(final int numZeroCounts, final double percentage) {
        return (int) Math.ceil(numZeroCounts * percentage / 100.0);
    }

    private static double safeLog2(final double x) {
        return x < EPSILON ? LN2_EPSILON : Math.log(x) * INV_LN2;
    }
}
