package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation;

import com.google.common.primitives.Doubles;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * TODO
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyRatioKernelSegmenter {
    private final List<SimpleInterval> intervals;
    private final List<Double> denoisedCopyRatios;

    public CopyRatioKernelSegmenter(final ReadCountCollection denoisedCopyRatioProfile) {
        Utils.nonNull(denoisedCopyRatioProfile);
        intervals = denoisedCopyRatioProfile.targets().stream().map(Target::getInterval).collect(Collectors.toList());
        denoisedCopyRatios = Doubles.asList(denoisedCopyRatioProfile.counts().getColumn(0));
    }

    public List<SimpleInterval> findSegments(final int maxNumSegmentsPerChromosome,
                                             final double kernelVariance,
                                             final int kernelApproximationDimension,
                                             final List<Integer> windowSizes,
                                             final double numChangepointsPenaltyLinearFactor,
                                             final double numChangepointsPenaltyLogLinearFactor) {
        Utils.validateArg(numChangepointsPenaltyLinearFactor == 0. || numChangepointsPenaltyLinearFactor >= 1.,
                "Linear factor for the penalty on the number of changepoints per chromosome must be either zero or greater than or equal to 1.");
        Utils.validateArg(numChangepointsPenaltyLogLinearFactor == 0. || numChangepointsPenaltyLogLinearFactor >= 1.,
                "Log-linear factor for the penalty on the number of changepoints per chromosome must be either zero or greater than or equal to 1.");

        return new ArrayList<>();
    }
}
