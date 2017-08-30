package org.broadinstitute.hellbender.tools.copynumber.utils.segmentation;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.function.BiFunction;

/**
 * Segment data (i.e., find multiple changepoints) using a method based on the kernel-segmentation algorithm
 * described in <a href="https://hal.inria.fr/hal-01413230/document">https://hal.inria.fr/hal-01413230/document</a>, which gives a framework to quickly calculate the cost of
 * a segment given a low-rank approximation to a specified kernel.  However, unlike the algorithm described there,
 * which seeks to minimize a global segmentation cost, the method implemented here instead finds candidate changepoints
 * based on the local costs within windows of various sizes.
 *
 * Given <i>N</i> data points to segment, the basic steps of the method are:
 *
 * <ol>
 *     1) Select <i>C<sub>max</sub></i>, the maximum number of changepoints to discover.
 * </ol>
 * <ol>
 *     2) Select a kernel (linear for sensitivity to changes in the distribution mean,
 *     Gaussian with a specified variance <i>σ<sup>2</sup></i> for multimodal data, etc.)
 *     and a subsample of <i>p</i> points to approximate it using singular value decomposition.
 * </ol>
 * <ol>
 *     3) Select window sizes <i>w<sub>j</sub></i> for which to compute local costs at each point.
 *     To be precise, we compute the cost of a changepoint at the point with index <i>i</i>,
 *     assuming adjacent segments containing the points with indices <i>[i - w<sub>j</sub> + 1, i]</i>
 *     and <i>[i + 1, i + w<sub>j</sub>]</i>.
 * </ol>
 * <ol>
 *     4) For each of these cost functions, find (up to) the  <i>C<sub>max</sub></i> most significant local minima.
 *     The problem of finding local minima of a noisy function can be solved by using topological persistence
 *     (e.g., <a href="https://people.mpi-inf.mpg.de/~weinkauf/notes/persistence1d.html">https://people.mpi-inf.mpg.de/~weinkauf/notes/persistence1d.html</a>
 *     and <a href="http://www2.iap.fr/users/sousbie/web/html/indexd3dd.html?post/Persistence-and-simplification>http://www2.iap.fr/users/sousbie/web/html/indexd3dd.html?post/Persistence-and-simplification</a>).
 * </ol>
 * <ol>
 *     5) These sets of local minima from all window sizes together provide the pool of candidate changepoints
 *     (some of which may overlap exactly or approximately).  We perform backwards selection using the global segmentation cost.
 *     That is, we compute the global segmentation cost given all the candidate changepoints,
 *     calculate the cost change for removing each of the changepoints individually,
 *     remove the changepoint with the minimum cost change, and repeat.
 *     This gives the global cost as a function of the number of changepoints <i>C</i>.
 * </ol>
 * <ol>
 *     6) Add a penalty <i>A * C + B * C * log N / C</i> to the global cost and find the minimum to determine the
 *     number of changepoints, where <i>A</i> and <i>B</i> are specified penalty factors.
 *
 * </ol>
 *
 * See discussion at <a href="https://github.com/broadinstitute/gatk/issues/2858#issuecomment-324125586">https://github.com/broadinstitute/gatk/issues/2858#issuecomment-324125586</a>
 * and accompanying plots for more detail.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class KernelSegmenter<T> {
    private static final int RANDOM_SEED = 1216;

    private final List<T> data;

    public KernelSegmenter(final List<T> data) {
        this.data = Collections.unmodifiableList(new ArrayList<T>(data));
    }

    public List<Integer> findChangepoints(final int maxNumChangepoints,
                                          final BiFunction<T, T, Double> kernel,
                                          final int kernelApproximationDimension,
                                          final List<Integer> windowSizes,
                                          final double numChangepointsPenaltyLinearFactor,
                                          final double numChangepointsPenaltyLogLinearFactor) {
        ParamUtils.isPositiveOrZero(maxNumChangepoints, "Maximum number of changepoints must be non-negative.");
        ParamUtils.isPositive(kernelApproximationDimension, "Dimension of kernel approximation must be positive.");
        Utils.validateArg(windowSizes.stream().allMatch(ws -> ws > 0), "Window sizes must all be positive.");
        Utils.validateArg(new HashSet<>(windowSizes).size() == windowSizes.size(), "Window sizes must all be unique.");
        Utils.validateArg(numChangepointsPenaltyLinearFactor == 0. || numChangepointsPenaltyLinearFactor >= 1.,
                "Linear factor for the penalty on the number of changepoints per chromosome must be either zero or greater than or equal to 1.");
        Utils.validateArg(numChangepointsPenaltyLogLinearFactor == 0. || numChangepointsPenaltyLogLinearFactor >= 1.,
                "Log-linear factor for the penalty on the number of changepoints per chromosome must be either zero or greater than or equal to 1.");

        if (maxNumChangepoints == 0) {
            return Collections.emptyList();
        }

        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));
        final RealMatrix reducedObservationMatrix = computeReducedObservationMatrix(rng, data, kernel, kernelApproximationDimension);
        final double[] kernelApproximationDiagonal = computeKernelApproximationDiagonal(reducedObservationMatrix);
        final List<Integer> changepointCandidates = findChangepointCandidates(
                data, reducedObservationMatrix, kernelApproximationDiagonal, maxNumChangepoints, windowSizes);
        return selectChangepoints(
                changepointCandidates, maxNumChangepoints, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor,
                reducedObservationMatrix, kernelApproximationDiagonal);
    }

    private final class Segment {
        private final int start;    //index of start point, inclusive
        private final int end;      //index of end point, inclusive
        private final double cost;

        private Segment(final int start,
                        final int end,
                        final double cost) {
            this.start = start;
            this.end = end;
            this.cost = cost;
        }

        private Segment(final int start,
                        final int end,
                        final RealMatrix reducedObservationMatrix,
                        final double[] kernelApproximationDiagonal) {
            this(start, end, computeSegmentCost(start, end, reducedObservationMatrix, kernelApproximationDiagonal));
        }
    }

    private static <T> RealMatrix computeReducedObservationMatrix(final RandomGenerator rng,
                                                                  final List<T> data,
                                                                  final BiFunction<T, T, Double> kernel,
                                                                  final double kernelApproximationDimension) {

    }

    private static double[] computeKernelApproximationDiagonal(final RealMatrix reducedObservationMatrix) {

    }

    private static <T> List<Integer> findChangepointCandidates(final List<T> data,
                                                               final RealMatrix reducedObservationMatrix,
                                                               final double[] kernelApproximationDiagonal,
                                                               final int maxNumChangepoints,
                                                               final List<Integer> windowSizes) {
        //warn if window sizes too large
    }

    private static List<Integer> selectChangepoints(final List<Integer> changepointCandidates,
                                                    final int maxNumChangepoints,
                                                    final double numChangepointsPenaltyLinearFactor,
                                                    final double numChangepointsPenaltyLogLinearFactor,
                                                    final RealMatrix reducedObservationMatrix,
                                                    final double[] kernelApproximationDiagonal) {
    }

    /**
     * Compute the cost of a segment.  This is defined by Eq. 11 of <a href="https://hal.inria.fr/hal-01413230/document">https://hal.inria.fr/hal-01413230/document</a>
     * (except we use the low-rank approximation to the kernel, as described in Sec. 3.2, ibid).
     * Various recursion relations are used to compute costs iteratively.
     * @param start inclusive start index of segment
     * @param end   inclusive end index of segment
     * @param reducedObservationMatrix      N x p matrix of projected observations, where N is the number of data points
     *                                      and p is the dimension of the low-rank approximation of the kernel matrix;
     *                                      this is the Z matrix described in the text preceding Eq. 14, ibid
     * @param kernelApproximationDiagonal   N diagonal terms of the low-rank approximation to the kernel matrix
     */
    private static double computeSegmentCost(final int start,
                                             final int end,
                                             final RealMatrix reducedObservationMatrix,
                                             final double[] kernelApproximationDiagonal) {

    }

    private static double[] computeWindowCosts(final RealMatrix reducedObservationMatrix,
                                               final double[] kernelApproximationDiagonal,
                                               final int windowSize) {

    }
}
