package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenterUnitTest;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class CopyRatioKernelSegmenterUnitTest {
    private static final int RANDOM_SEED = 1;   //reset seed before each simulated test case

    /**
     * Generates same Gaussian test data as {@link KernelSegmenterUnitTest#dataKernelSegmenter()},
     * but introduces further segments by placing data on different chromosomes.
     */
    @DataProvider(name = "dataCopyRatioKernelSegmenter")
    public Object[][] dataCopyRatioKernelSegmenter() {
        final int numPoints = 1000;

        final Random rng = new Random(RANDOM_SEED);
        rng.setSeed(RANDOM_SEED);
        final List<Double> dataGaussian = IntStream.range(0, numPoints).boxed()
                .map(i -> Math.abs(i / 100 - 5) + 0.1 * rng.nextGaussian())
                .collect(Collectors.toList());              //changepoints at 99, 199, 299, 399, 499, 599, 699, 799, 899

        final List<SimpleInterval> intervals = IntStream.range(0, numPoints).boxed()
                .map(i -> new SimpleInterval(
                        Integer.toString(i / 250 + 1),  //start a new chromosome every 250 points, which adds additional changepoints
                        (i % 250) * 10 + 1,
                        (i % 250) * 10 + 10))         //intervals for copy-ratio data points have length = 10
                .collect(Collectors.toList());

        final ReadCountCollection denoisedCopyRatioProfile = new ReadCountCollection(
                intervals.stream().map(Target::new).collect(Collectors.toList()),
                Collections.singletonList("testSampleName"),
                new Array2DRowRealMatrix(dataGaussian.stream().mapToDouble(Double::doubleValue).toArray())
        );

        final CopyRatioSegmentationResult segmentsExpected =
                new CopyRatioSegmentationResult(Arrays.asList(
                        new CopyRatioSegmentationResult.CopyRatioSegment(new SimpleInterval("1", 1, 1000), dataGaussian.subList(0, 100)),
                        new CopyRatioSegmentationResult.CopyRatioSegment(new SimpleInterval("1", 1001, 2000), dataGaussian.subList(100, 200)),
                        new CopyRatioSegmentationResult.CopyRatioSegment(new SimpleInterval("1", 2001, 2500), dataGaussian.subList(200, 250)),
                        new CopyRatioSegmentationResult.CopyRatioSegment(new SimpleInterval("2", 1, 500), dataGaussian.subList(250, 300)),
                        new CopyRatioSegmentationResult.CopyRatioSegment(new SimpleInterval("2", 501, 1500), dataGaussian.subList(300, 400)),
                        new CopyRatioSegmentationResult.CopyRatioSegment(new SimpleInterval("2", 1501, 2500), dataGaussian.subList(400, 500)),
                        new CopyRatioSegmentationResult.CopyRatioSegment(new SimpleInterval("3", 1, 1000), dataGaussian.subList(500, 600)),
                        new CopyRatioSegmentationResult.CopyRatioSegment(new SimpleInterval("3", 1001, 2000), dataGaussian.subList(600, 700)),
                        new CopyRatioSegmentationResult.CopyRatioSegment(new SimpleInterval("3", 2001, 2500), dataGaussian.subList(700, 750)),
                        new CopyRatioSegmentationResult.CopyRatioSegment(new SimpleInterval("4", 1, 500), dataGaussian.subList(750, 800)),
                        new CopyRatioSegmentationResult.CopyRatioSegment(new SimpleInterval("4", 501, 1500), dataGaussian.subList(800, 900)),
                        new CopyRatioSegmentationResult.CopyRatioSegment(new SimpleInterval("4", 1501, 2500), dataGaussian.subList(900, 1000))));

        return new Object[][]{
                {denoisedCopyRatioProfile, segmentsExpected}
        };
    }

    @Test(dataProvider = "dataCopyRatioKernelSegmenter")
    public void testCopyRatioKernelSegmenter(final ReadCountCollection denoisedCopyRatioProfile,
                                             final CopyRatioSegmentationResult segmentsExpected) {
        final int maxNumChangepointsPerChromosome = 25;
        final double kernelVariance = 0.;
        final int kernelApproximationDimension = 20;
        final List<Integer> windowSizes = Arrays.asList(8, 16, 32, 64);
        final double numChangepointsPenaltyLinearFactor = 2.;
        final double numChangepointsPenaltyLogLinearFactor = 2.;

        final CopyRatioSegmentationResult segments = new CopyRatioKernelSegmenter(denoisedCopyRatioProfile)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVariance, kernelApproximationDimension,
                        windowSizes, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor);
        Assert.assertEquals(segments, segmentsExpected);
    }

}