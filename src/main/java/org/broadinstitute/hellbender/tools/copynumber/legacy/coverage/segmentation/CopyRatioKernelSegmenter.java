package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * TODO
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyRatioKernelSegmenter {
    //Gaussian kernel for a specified variance; if variance is zero, use a linear kernel
    private static final Function<Double, BiFunction<Double, Double, Double>> kernel =
            variance -> variance == 0.
                    ? (x, y) -> x * y
                    : (x, y) -> Math.exp(-(x - y) * (x - y) / variance);

    private final List<SimpleInterval> intervals;
    private final Map<String, List<Double>> denoisedCopyRatiosPerChromosome;    //in log2 space

    /**
     * @param denoisedCopyRatioProfile  in log2 space
     */
    public CopyRatioKernelSegmenter(final ReadCountCollection denoisedCopyRatioProfile) {
        Utils.nonNull(denoisedCopyRatioProfile);
        intervals = denoisedCopyRatioProfile.targets().stream().map(Target::getInterval).collect(Collectors.toList());
        denoisedCopyRatiosPerChromosome = IntStream.range(0, denoisedCopyRatioProfile.targets().size()).boxed()
                .map(i -> new ImmutablePair<>(
                        intervals.get(i).getContig(),
                        denoisedCopyRatioProfile.getRow(0)[i]))
                .collect(Collectors.groupingBy(Pair::getKey, Collectors.mapping(Pair::getValue, Collectors.toList())));
    }

    public CopyRatioSegmentationResult findSegmentation(final int maxNumChangepointsPerChromosome,
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

        //loop over chromosomes, find changepoints, and create copy-ratio segments
        final List<CopyRatioSegmentationResult.CopyRatioSegment> segments = new ArrayList<>();
        int numDataPointsProcessed = 0;
        for (final String chromosome : denoisedCopyRatiosPerChromosome.keySet()) {
            final List<Double> denoisedCopyRatiosInChromosome = denoisedCopyRatiosPerChromosome.get(chromosome);

            final List<Integer> changepoints = new KernelSegmenter<>(denoisedCopyRatiosInChromosome)
                .findChangepoints(maxNumChangepointsPerChromosome, kernel.apply(kernelVariance), kernelApproximationDimension,
                        windowSizes, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor);

            if (!changepoints.contains(denoisedCopyRatiosInChromosome.size())) {
                changepoints.add(denoisedCopyRatiosInChromosome.size());
            }
            for (final int changepoint : changepoints) {
                final int start = intervals.get(numDataPointsProcessed).getStart();
                final int end = intervals.get(numDataPointsProcessed + changepoint).getEnd();
                final List<Double> denoisedCopyRatiosInSegment = denoisedCopyRatiosInChromosome.subList(
                        numDataPointsProcessed, numDataPointsProcessed + changepoint + 1);
                segments.add(new CopyRatioSegmentationResult.CopyRatioSegment(
                        new SimpleInterval(chromosome, start, end),
                        denoisedCopyRatiosInSegment));
                numDataPointsProcessed += denoisedCopyRatiosInSegment.size();
            }
        }
        return new CopyRatioSegmentationResult(segments);
    }
}
