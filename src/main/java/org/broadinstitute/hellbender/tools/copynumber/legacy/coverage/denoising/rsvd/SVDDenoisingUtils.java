package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.rsvd;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.gcbias.GCCorrector;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MatrixSummaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

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

    public static final double EPSILON = 1E-30;
    private static final double INV_LN2 = GATKProtectedMathUtils.INV_LOG_2;
    private static final double LN2_EPSILON = Math.log(EPSILON) * INV_LN2;

    private static final int NUM_SLICES_SPARK = 50;

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
        final RealMatrix denoisedCounts = subtractProjection(standardizedCounts,
                panelOfNormals.getLeftSingular(), panelOfNormals.getLeftSingularPseudoinverse(), numEigensamples);

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
     * Standardize read counts for a single sample, using interval fractional medians from a panel of normals.
     * The original {@code readCounts} has dimensions intervals x 1 and is not modified.
     */
    private static RealMatrix preprocessAndStandardizeSample(final SVDReadCountPanelOfNormals panelOfNormals,
                                                             final RealMatrix readCounts) {
        RealMatrix result = readCounts.copy();

        logger.info("Transforming read counts to fractional coverage...");
        final double[] sampleSums = GATKProtectedMathUtils.columnSums(readCounts);
        result.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int intervalIndex, int sampleIndex, double value) {
                return value / sampleSums[sampleIndex];
            }
        });

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
        final double[] sampleMedians = MatrixSummaryUtils.getColumnMedians(result);
        result.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int intervalIndex, int sampleIndex, double value) {
                return safeLog2(value / sampleMedians[sampleIndex]);
            }
        });

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
     * Given standardized read counts specified by a column vector S (dimensions {@code M x 1}),
     * left-singular vectors U (dimensions {@code M x K}), and the pseudoinverse U<sup>+</sup> (dimensions {@code K x M}),
     * returns s - U U<sup>+</sup> s.
     */
    private static RealMatrix subtractProjection(final RealMatrix standardizedProfile,
                                                 final double[][] leftSingular,
                                                 final double[][] leftSingularPseudoinverse,
                                                 final int numEigensamples) {
        final int numIntervals = leftSingular.length;

        logger.info("Distributing the standardized read counts...");

        logger.info("Composing left-singular matrix and pseudoinverse for the requested number of eigensamples and transposing them...");
        final RealMatrix leftSingularTruncatedTransposedMatrix = new Array2DRowRealMatrix(leftSingular)
                .getSubMatrix(0, numIntervals - 1, 0, numEigensamples - 1)
                .transpose();
        final RealMatrix leftSingularPseudoinverseTruncatedTransposedMatrix = new Array2DRowRealMatrix(leftSingularPseudoinverse)
                .getSubMatrix(0, numEigensamples - 1, 0, numIntervals - 1)
                .transpose();

        logger.info("Computing projection of transpose...");
        final RealMatrix projectionTranspose = standardizedProfile.transpose()
                .multiply(leftSingularPseudoinverseTruncatedTransposedMatrix)
                .multiply(leftSingularTruncatedTransposedMatrix);

        logger.info("Subtracting projection...");
        return standardizedProfile.subtract(projectionTranspose.transpose());
    }

    static double safeLog2(final double x) {
        return x < EPSILON ? LN2_EPSILON : Math.log(x) * INV_LN2;
    }
}
