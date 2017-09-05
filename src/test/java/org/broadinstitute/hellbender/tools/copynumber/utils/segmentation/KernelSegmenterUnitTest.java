package org.broadinstitute.hellbender.tools.copynumber.utils.segmentation;

import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
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
public class KernelSegmenterUnitTest extends BaseTest {
    private static final int RANDOM_SEED = 1;   //reset seed before each simulated test case

    @DataProvider(name = "dataKernelSegmenter")
    public Object[][] dataKernelSegmenter() {
        final int numPoints = 1000;

        final BiFunction<Double, Double, Double> linearKernel = (x, y) -> x * y;
        final BiFunction<Double, Double, Double> gaussianKernel = (x, y) -> Math.exp(-(x - y) * (x - y));

        final Random rng = new Random(RANDOM_SEED);
        rng.setSeed(RANDOM_SEED);
        final List<Double> dataGaussian = IntStream.range(0, numPoints).boxed()
                .map(i -> Math.abs(i / 100 - 5) + 0.1 * rng.nextGaussian())
                .collect(Collectors.toList());
        final List<Integer> changepointsExpectedGaussian = Arrays.asList(
                99, 199, 299, 399, 499, 599, 699, 799, 899, 999);

        return new Object[][]{
                {dataGaussian, linearKernel, changepointsExpectedGaussian}
        };
    }

    @Test(dataProvider = "dataKernelSegmenter")
    public void testKernelSegmenter(final List<Double> data,
                                    final BiFunction<Double, Double, Double> kernel,
                                    final List<Integer> changepointsExpected) {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
        final int maxNumChangepoints = 25;
        final int kernelApproximationDimension = 20;
        final List<Integer> windowSizes = Arrays.asList(8, 16, 32, 64);
        final double numChangepointsPenaltyLinearFactor = 2.;
        final double numChangepointsPenaltyLogLinearFactor = 2.;

        System.out.println(data);
        final KernelSegmenter<Double> segmenter = new KernelSegmenter<>(data);
        Assert.assertEquals(
                segmenter.findChangepoints(maxNumChangepoints, kernel, kernelApproximationDimension, windowSizes,
                        numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor),
                changepointsExpected);
    }
}