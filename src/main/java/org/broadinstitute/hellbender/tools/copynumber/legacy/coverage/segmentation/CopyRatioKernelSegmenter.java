package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * TODO
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyRatioKernelSegmenter {
    //we use a linear kernel for segmentation so that we are sensitive to changes in the mean of the copy-ratio distribution
    private static final BiFunction<Double, Double, Double> linearKernel = (x, y) -> x * y;

    private final Map<String, List<Double>> denoisedCopyRatiosPerChromosome;

    public CopyRatioKernelSegmenter(final ReadCountCollection denoisedCopyRatioProfile) {
        Utils.nonNull(denoisedCopyRatioProfile);
        denoisedCopyRatiosPerChromosome = IntStream.range(0, denoisedCopyRatioProfile.targets().size()).boxed()
                .map(i -> new ImmutablePair<>(
                        denoisedCopyRatioProfile.targets().get(i).getInterval().getContig(),
                        denoisedCopyRatioProfile.getRow(0)[i]))
                .collect(Collectors.groupingBy(Pair::getKey, Collectors.mapping(Pair::getValue, Collectors.toList())));
    }

    public List<SimpleInterval> findSegments(final int maxNumChangepointsPerChromosome,
                                             final double kernelVariance,
                                             final int kernelApproximationDimension,
                                             final List<Integer> windowSizes,
                                             final double numChangepointsPenaltyLinearFactor,
                                             final double numChangepointsPenaltyLogLinearFactor) {
        ParamUtils.isPositiveOrZero(maxNumChangepointsPerChromosome, "Maximum number of changepoints must be non-negative.");
        ParamUtils.isPositiveOrZero(kernelVariance, "Variance of Gaussian kernel must be non-negative (if zero, a linear kernel will be used).");
        ParamUtils.isPositive(kernelApproximationDimension, "Dimension of kernel approximation must be positive.");
        Utils.validateArg(windowSizes.stream().allMatch(ws -> ws > 0), "Window sizes must all be positive.");
        Utils.validateArg(new HashSet<>(windowSizes).size() == windowSizes.size(), "Window sizes must all be unique.");
        Utils.validateArg(numChangepointsPenaltyLinearFactor == 0. || numChangepointsPenaltyLinearFactor >= 1.,
                "Linear factor for the penalty on the number of changepoints per chromosome must be either zero or greater than or equal to 1.");
        Utils.validateArg(numChangepointsPenaltyLogLinearFactor == 0. || numChangepointsPenaltyLogLinearFactor >= 1.,
                "Log-linear factor for the penalty on the number of changepoints per chromosome must be either zero or greater than or equal to 1.");

        //loop over chromosomes and find changepoints
        final Map<String, List<Integer>> changepointsPerChromosome = new HashMap<>();
        for (final String chromosome : denoisedCopyRatiosPerChromosome.keySet()) {
            final List<Integer> changepoints = new KernelSegmenter<Double>(denoisedCopyRatiosPerChromosome.get(chromosome))
                .findChangepoints(maxNumChangepointsPerChromosome, linearKernel, kernelApproximationDimension,
                        windowSizes, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor);
        }
        return new ArrayList<>();
    }
}
