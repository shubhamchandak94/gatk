package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation;

import com.google.common.primitives.Doubles;
import javafx.util.Pair;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * TODO
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyRatioKernelSegmenter {
    private final List<SimpleInterval> intervals;
    private final Map<String, List<Double>> denoisedCopyRatiosPerChromosome;

    public CopyRatioKernelSegmenter(final ReadCountCollection denoisedCopyRatioProfile) {
        Utils.nonNull(denoisedCopyRatioProfile);
        intervals = denoisedCopyRatioProfile.targets().stream().map(Target::getInterval).collect(Collectors.toList());
//        denoisedCopyRatiosPerChromosome = IntStream.range(0, intervals.size())
//                .collect(Collectors.groupingBy()); //TODO
        denoisedCopyRatiosPerChromosome = new HashMap<>();
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

        //loop over chromosomes (optionally?)
//        final List<Integer> changepoints = new KernelSegmenter<Double>(denoisedCopyRatios)
//                .findChangepoints(maxNumChangepointsPerChromosome, kernel, kernelApproximationDimension,
//                        windowSizes, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor);
        return new ArrayList<>();
    }
}
