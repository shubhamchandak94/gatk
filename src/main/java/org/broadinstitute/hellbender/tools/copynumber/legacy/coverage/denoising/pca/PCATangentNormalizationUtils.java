package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.pca;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.mllib.linalg.DenseMatrix;
import org.apache.spark.mllib.linalg.Matrix;
import org.apache.spark.mllib.linalg.distributed.RowMatrix;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.CaseToPoNTargetMapper;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.spark.SparkConverter;

import java.util.List;
import java.util.stream.IntStream;

/**
 * Utility class for package-private methods for performing tangent normalization (and related operations).
 *
 * Currently, only supports tangent normalization in the reduced hyperplane, not the logNormal hyperplane
 */
public final class PCATangentNormalizationUtils {
    private static final Logger logger = LogManager.getLogger(PCATangentNormalizationUtils.class);

    /**
     * Minimum target normalized and column centered count possible.
     *
     * <p>
     *     It must be small yet greater than 0 to avoid -Inf problems in the calculations.
     * </p>
     */
    public static final double EPSILON = 1E-9;
    private static final double INV_LN2 = GATKProtectedMathUtils.INV_LOG_2;
    private static final double LOG_2_EPSILON = Math.log(EPSILON) * INV_LN2;

    private static final int TN_NUM_SLICES_SPARK = 50;

    private PCATangentNormalizationUtils() {}

    /**
     * Target-factor normalizes a {@link RealMatrix} in-place given target factors..
     */
    static void factorNormalize(final RealMatrix input, final double[] targetFactors) {
        Utils.nonNull(input, "Input matrix cannot be null.");
        Utils.nonNull(targetFactors, "Target factors cannot be null.");
        Utils.validateArg(targetFactors.length == input.getRowDimension(),
                "Number of target factors does not correspond to the number of rows.");
        // Divide all counts by the target factor for the row.
        input.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int row, final int column, final double value) {
                return value / targetFactors[row];
            }
        });
    }

    /**
     *  Do the full tangent normalization process given a {@link PCACoveragePoN} and a proportional-coverage profile.
     *
     *  This includes:
     *   <ul><li>normalization by target factors (optional)</li>
     *   <li>projection of the normalized coverage profile into the hyperplane from the PoN</li>
     *   </ul>
     *
     * @param pon -- never {@code null}
     * @param profile -- never {@code null}.  Must contain data for at least one sample.
     * @param ctx spark context.  Use {@code null} if no context is available
     * @return never {@code null}
     */
    static PCATangentNormalizationResult tangentNormalize(final PCACoveragePoN pon,
                                                          final ReadCountCollection profile,
                                                          final JavaSparkContext ctx) {
        Utils.nonNull(pon, "PoN cannot be null.");
        Utils.nonNull(profile, "Proportional coverages cannot be null.");
        ParamUtils.isPositive(profile.columnNames().size(), "Column names cannot be an empty list.");

        final ReadCountCollection factorNormalizedCoverage = mapTargetsToPoNAndFactorNormalize(profile, pon);

        return tangentNormalize(factorNormalizedCoverage, pon.getPanelTargetNames(), pon.getReducedPanelCounts(), pon.getReducedPanelPInverseCounts(), ctx);
    }

    /**
     * Returns a target-factor-normalized {@link ReadCountCollection} given a {@link PCACoveragePoN}..
     */
    private static ReadCountCollection mapTargetsToPoNAndFactorNormalize(final ReadCountCollection input, final PCACoveragePoN pon) {
        final CaseToPoNTargetMapper targetMapper = new CaseToPoNTargetMapper(input.targets(), pon.getTargetNames());
        final RealMatrix inputCounts = targetMapper.fromCaseToPoNCounts(input.counts());
        factorNormalize(inputCounts, pon.getAllIntervalProportionalMedians());   //factor normalize in-place
        return targetMapper.fromPoNtoCaseCountCollection(inputCounts, input.columnNames());
    }

    /**
     * Tangent normalize given the raw PoN data.
     *
     *  Ahat^T = (C^T P^T) A^T
     *  Therefore, C^T is the RowMatrix
     *
     *  pinv: P
     *  panel: A
     *  projection: Ahat
     *  cases: C
     *  betahat: C^T P^T
     *  tangentNormalizedCounts: C - Ahat
     */
    private static PCATangentNormalizationResult tangentNormalize(final ReadCountCollection targetFactorNormalizedCounts,
                                                                  final List<String> panelTargetNames,
                                                                  final RealMatrix reducedPanelCounts,
                                                                  final RealMatrix reducedPanelPInvCounts,
                                                                  final JavaSparkContext ctx) {
        final CaseToPoNTargetMapper targetMapper = new CaseToPoNTargetMapper(targetFactorNormalizedCounts.targets(), panelTargetNames);

        // The input counts with rows (targets) sorted so that they match the PoN's order.
        final RealMatrix tangentNormalizationRawInputCounts = targetMapper.fromCaseToPoNCounts(targetFactorNormalizedCounts.counts());

        // We prepare the counts for tangent normalization.
        final RealMatrix tangentNormalizationInputCounts = composeTangentNormalizationInputMatrix(tangentNormalizationRawInputCounts);

        // Make the C^T a distributed matrix (RowMatrix)
        final RowMatrix caseTDistMat = SparkConverter.convertRealMatrixToSparkRowMatrix(
                ctx, tangentNormalizationInputCounts.transpose(), TN_NUM_SLICES_SPARK);

        // Spark local matrices (transposed)
        final Matrix pinvTLocalMat = new DenseMatrix(
                reducedPanelPInvCounts.getRowDimension(), reducedPanelPInvCounts.getColumnDimension(),
                Doubles.concat(reducedPanelPInvCounts.getData()), true).transpose();
        final Matrix panelTLocalMat = new DenseMatrix(
                reducedPanelCounts.getRowDimension(), reducedPanelCounts.getColumnDimension(),
                Doubles.concat(reducedPanelCounts.getData()), true).transpose();

        // Calculate the projection transpose in a distributed matrix, then convert to Apache Commons matrix (not transposed)
        final RowMatrix betahatDistMat = caseTDistMat.multiply(pinvTLocalMat);
        final RowMatrix projectionTDistMat = betahatDistMat.multiply(panelTLocalMat);
        final RealMatrix projection = SparkConverter.convertSparkRowMatrixToRealMatrix(
                projectionTDistMat, tangentNormalizationInputCounts.transpose().getRowDimension()).transpose();

        // Subtract the projection from the cases
        final RealMatrix tangentNormalizedCounts = tangentNormalizationInputCounts.subtract(projection);

        // Construct the result object and return it with the correct targets.
        final ReadCountCollection tangentNormalized = targetMapper.fromPoNtoCaseCountCollection(
                tangentNormalizedCounts, targetFactorNormalizedCounts.columnNames());
        final ReadCountCollection preTangentNormalized = targetMapper.fromPoNtoCaseCountCollection(
                tangentNormalizationInputCounts, targetFactorNormalizedCounts.columnNames());

        return new PCATangentNormalizationResult(tangentNormalized, preTangentNormalized);
    }

    /**
     * Prepares the data to perform tangent normalization.
     * <p>
     * This is done by count group or column:
     *   <ol>
     *     </li>we divide counts by the column mean,</li>
     *     </li>then we transform value to their log_2,</li>
     *     </li>and finally we center them around the median.</li>
     *   </ol>
     * </p>
     *
     * @param matrix input matrix.
     * @return never {@code null}.
     */
    private static RealMatrix composeTangentNormalizationInputMatrix(final RealMatrix matrix) {
        final RealMatrix result = matrix.copy();

        // step 1: divide by column means and log_2 transform
        final double[] columnMeans = GATKProtectedMathUtils.columnMeans(matrix);
        result.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int row, final int column, final double value) {
                return truncatedLog2(value / columnMeans[column]);
            }
        });

        // step 2: subtract column medians
        final double[] columnMedians = IntStream.range(0, matrix.getColumnDimension())
                .mapToDouble(c -> new Median().evaluate(result.getColumn(c))).toArray();
        result.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int row, final int column, final double value) {
                return value - columnMedians[column];
            }
        });

        return result;
    }

    private static double truncatedLog2(final double x) {
        return x < EPSILON ? LOG_2_EPSILON : Math.log(x) * INV_LN2;
    }
}
