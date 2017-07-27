package org.broadinstitute.hellbender.tools.copynumber.temporary;

import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

public class HDF5ReadCountCollectionUnitTest extends BaseTest {

    @Test
    public void basicTest() {
        final File outputFile = createTempFile("HDF5ReadCountCollection", ".cov");
        final List<Target> newTargets = new ArrayList<>();
        newTargets.add(new Target(new SimpleInterval("1", 1000, 2000)));
        newTargets.add(new Target(new SimpleInterval("1", 5000, 6000)));

        final List<String> sampleNames = new ArrayList<>();
        sampleNames.add("SAMPLE1");

        final double[][] newTargetValues = {{2, 10}};

        // The output file already exists at this point, since it is a temp file.
        HDF5ReadCountCollection.write(outputFile, newTargets, newTargetValues, sampleNames);

        final ReadCountCollection rcc = ReadCountCollectionUtils.parseHdf5AsDouble(outputFile);
        Assert.assertEquals(rcc.columnNames(), sampleNames);
        Assert.assertEquals(rcc.targets(), newTargets);
        Assert.assertFalse(rcc.targets() == newTargets);
        Assert.assertFalse(rcc.columnNames() == sampleNames);

        Assert.assertTrue(IntStream.range(0, newTargetValues.length).allMatch(i -> newTargetValues[i][0] == rcc.counts().getData()[0][i]));

        Assert.assertEquals(rcc.counts().transpose().getData().length, newTargetValues.length);
        Assert.assertEquals(rcc.counts().transpose().getData()[0].length, newTargetValues[0].length);
        Assert.assertFalse(rcc.counts().transpose().getData() == newTargetValues);
    }
}
