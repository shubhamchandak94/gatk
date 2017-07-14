package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.rsvd;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.mllib.linalg.DenseMatrix;
import org.apache.spark.mllib.linalg.Matrix;
import org.apache.spark.mllib.linalg.distributed.RowMatrix;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.DenoisedCopyRatioResult;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.spark.SparkConverter;

import java.util.HashSet;
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
    static DenoisedCopyRatioResult denoise(final SVDReadCountPanelOfNormals panelOfNormals,
                                           final ReadCountCollection readCounts,
                                           final int numEigensamples,
                                           final JavaSparkContext ctx) {
        Utils.nonNull(panelOfNormals);
        validateReadCounts(readCounts);
        Utils.nonNull(ctx);
        ParamUtils.isPositive(numEigensamples, "Number of eigensamples to use for denoising must be positive.");
        Utils.validateArg(numEigensamples <= panelOfNormals.getNumEigensamples(),
                "Number of eigensamples to use for denoising is greater than the number available in the panel of normals.");

        logger.info("Validating sample intervals against original intervals used to build panel of normals...");
        Utils.validateArg(panelOfNormals.getOriginalIntervals().equals(readCounts.targets()),
                "Sample intervals must be identical to the original intervals used to build the panel of normals.");

        logger.info("Subsetting sample intervals to post-filter panel intervals...");
        final RealMatrix subsetReadCounts = readCounts.subsetTargets(new HashSet<>(panelOfNormals.getPanelIntervals())).counts();

        logger.info("Standardizing sample read counts...");
        final RealMatrix standardizedCounts = standardizeSample(subsetReadCounts, panelOfNormals.getPanelIntervalFractionalMedians());

        logger.info(String.format("Using %d out of %d eigensamples to denoise...", numEigensamples, panelOfNormals.getNumEigensamples()));
        logger.info("Subtracting projection onto space spanned by eigensamples...");
        final RealMatrix denoisedCounts = subtractProjection(standardizedCounts,
                panelOfNormals.getRightSingular(), panelOfNormals.getRightSingularPseudoinverse(), numEigensamples, ctx);

        logger.info("Sample denoised.");

        //construct the result
        final ReadCountCollection standardizedProfile = new ReadCountCollection(readCounts.targets(), readCounts.columnNames(), standardizedCounts);
        final ReadCountCollection denoisedProfile = new ReadCountCollection(readCounts.targets(), readCounts.columnNames(), denoisedCounts);

        return new DenoisedCopyRatioResult(standardizedProfile, denoisedProfile);
    }

    /**
     * Standardize read counts from a panel of normals.  This is performed in place to minimize memory footprint.
     */
    static void standardizePanel(final RealMatrix readCounts) {
        //TODO add filtering, etc. and extract

        logger.info("Transforming read counts to fractional coverage...");
        final double[] sampleSums = GATKProtectedMathUtils.columnSums(readCounts);
        readCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int intervalIndex, final int sampleIndex, final double value) {
                return value / sampleSums[sampleIndex];
            }
        });

        //use interval fractional medians from panel of normals if provided, else calculate from provided read counts
        logger.info("Dividing by interval medians...");
        final double[] intervalMedians = IntStream.range(0, readCounts.getColumnDimension())
                        .mapToDouble(intervalIndex -> new Median().evaluate(readCounts.getColumn(intervalIndex))).toArray();
        readCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int intervalIndex, final int sampleIndex, final double value) {
                return value / intervalMedians[intervalIndex];
            }
        });

        logger.info("Dividing by sample means and transforming to log2 space...");
        final double[] sampleMeans = GATKProtectedMathUtils.columnMeans(readCounts);
        readCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int intervalIndex, final int sampleIndex, final double value) {
                return safeLog2(value / sampleMeans[sampleIndex]);
            }
        });

        logger.info("Subtracting sample medians...");
        final double[] sampleMedians = IntStream.range(0, readCounts.getRowDimension())
                .mapToDouble(sampleIndex -> new Median().evaluate(readCounts.getColumn(sampleIndex))).toArray();
        readCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int intervalIndex, final int sampleIndex, final double value) {
                return value - sampleMedians[sampleIndex];
            }
        });

        logger.info("Panel read counts standardized.");
    }

    /**
     * Standardize read counts for a single sample, using interval fractional medians from a panel of normals.
     */
    private static RealMatrix standardizeSample(final RealMatrix readCounts, final double[] intervalFractionalMedians) {
        final RealMatrix result = readCounts.copy();

        logger.info("Transforming read counts to fractional coverage...");
        final double[] sampleSums = GATKProtectedMathUtils.columnSums(readCounts);
        readCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int intervalIndex, final int sampleIndex, final double value) {
                return value / sampleSums[sampleIndex];
            }
        });

        //use interval fractional medians from panel of normals if provided, else calculate from provided read counts
        logger.info("Dividing by interval medians from the panel of normals...");
        readCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int intervalIndex, final int sampleIndex, final double value) {
                return value / intervalFractionalMedians[intervalIndex];
            }
        });

        logger.info("Dividing by sample mean and transforming to log2 space...");
        final double[] sampleMeans = GATKProtectedMathUtils.columnMeans(readCounts);
        readCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int intervalIndex, final int sampleIndex, final double value) {
                return safeLog2(value / sampleMeans[sampleIndex]);
            }
        });

        logger.info("Subtracting sample median...");
        final double[] sampleMedians = IntStream.range(0, readCounts.getRowDimension())
                .mapToDouble(sampleIndex -> new Median().evaluate(readCounts.getColumn(sampleIndex))).toArray();
        readCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int intervalIndex, final int sampleIndex, final double value) {
                return value - sampleMedians[sampleIndex];
            }
        });

        logger.info("Sample read counts standardized.");

        return result;
    }

    /**
     * Given standardized read counts specified by a column vector S (dimensions {@code M x 1}),
     * right-singular vectors V (dimensions {@code M x K}), and the pseudoinverse V<sup>+</sup> (dimensions {@code K x M}),
     * returns s - V V<sup>+</sup> s.
     */
    private static RealMatrix subtractProjection(final RealMatrix standardizedProfile,
                                                 final double[][] rightSingular,
                                                 final double[][] rightSingularPseudoinverse,
                                                 final int numEigensamples,
                                                 final JavaSparkContext ctx) {
        final int numIntervals = rightSingular.length;

        logger.info("Distributing the standardized read counts...");
        final RowMatrix standardizedProfileDistributedMatrix = SparkConverter.convertRealMatrixToSparkRowMatrix(
                ctx, standardizedProfile.transpose(), NUM_SLICES_SPARK);

        logger.info("Composing right-singular matrix and pseudoinverse for the requested number of eigensamples and transposing them...");
        final RealMatrix rightSingularTruncatedMatrix = new Array2DRowRealMatrix(rightSingular)
                .getSubMatrix(0, numIntervals, 0, numEigensamples);
        final Matrix rightSingularTruncatedTransposedLocalMatrix = new DenseMatrix(numIntervals, numEigensamples,
                Doubles.concat(rightSingularTruncatedMatrix.getData()))
                .transpose();
        final RealMatrix rightSingularPseudoinverseTruncatedMatrix = new Array2DRowRealMatrix(rightSingularPseudoinverse)
                .getSubMatrix(0, numEigensamples, 0, numIntervals);
        final Matrix rightSingularPseudoinverseTransposeLocalMatrix = new DenseMatrix(numEigensamples, numIntervals,
                Doubles.concat(rightSingularPseudoinverseTruncatedMatrix.getData()), true)
                .transpose();

        logger.info("Computing projection of transpose...");
        final RowMatrix projectionTransposeDistributedMatrix = standardizedProfileDistributedMatrix
                .multiply(rightSingularPseudoinverseTransposeLocalMatrix)
                .multiply(rightSingularTruncatedTransposedLocalMatrix);

        logger.info("Converting and transposing result...");
        final int numRows = standardizedProfile.getRowDimension();
        final RealMatrix projection = SparkConverter.convertSparkRowMatrixToRealMatrix(projectionTransposeDistributedMatrix, numRows)
                .transpose();

        logger.info("Subtracting projection...");
        return standardizedProfile.subtract(projection);
    }

    private static double safeLog2(final double x) {
        return x < EPSILON ? LN2_EPSILON : Math.log(x) * INV_LN2;
    }
}
