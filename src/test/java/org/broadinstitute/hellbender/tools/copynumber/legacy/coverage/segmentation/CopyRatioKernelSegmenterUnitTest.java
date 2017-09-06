package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation;

import com.google.common.collect.ImmutableSet;
import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenterUnitTest;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.testng.Assert.*;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class CopyRatioKernelSegmenterUnitTest {
    private static final int RANDOM_SEED = 1;   //reset seed before each simulated test case

    /**
     * Generates same Gaussian test data as {@link KernelSegmenterUnitTest#dataKernelSegmenter()}
     * but introduces further segments by placing data on different chromosomes
     */
    @DataProvider(name = "dataCopyRatioKernelSegmenter")
    public Object[][] dataCopyRatioKernelSegmenter() {
        final int numPoints = 1000;

        final Random rng = new Random(RANDOM_SEED);
        rng.setSeed(RANDOM_SEED);
        final List<Double> dataGaussian = IntStream.range(0, numPoints).boxed()
                .map(i -> Math.abs(i / 100 - 5) + 0.1 * rng.nextGaussian())
                .collect(Collectors.toList());  //changepoints: 299, 699, 99, 899, 399, 199, 599, 499, 799

        return new Object[][]{
                {denoisedCopyRatioProfile, segmentsExpected}
        };
    }

    @Test(dataProvider = "dataCopyRatioKernelSegmenter")
    public void testCopyRatioKernelSegmenter(final ReadCountCollection denoisedCopyRatioProfile,
                                             final CopyRatioSegmentationResult segmentsExpected) {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
        final int maxNumChangepointsPerChromosome = 25;
        final double kernelVariance = 0.;
        final int kernelApproximationDimension = 20;
        final List<Integer> windowSizes = Arrays.asList(8, 16, 32, 64);
        final double numChangepointsPenaltyLinearFactor = 2.;
        final double numChangepointsPenaltyLogLinearFactor = 2.;

        final CopyRatioSegmentationResult segments = new CopyRatioKernelSegmenter(denoisedCopyRatioProfile)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVariance, kernelApproximationDimension,
                        ImmutableSet.copyOf(windowSizes).asList(),
                        numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor);
        Assert.assertEquals(segments, segmentsExpected);
    }

}