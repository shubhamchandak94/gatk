package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Created by tsato on 8/1/17.
 */
public class ContextDependentArtifactFilterIntegrationTest extends CommandLineProgramTest {

    @Test
    public void test() {
        final File table = createTempFile("contamination", ".table");
        final String[] args = {
                "-R", v37_chr17_1Mb_Reference,
                "-I", NA12878_chr17_1k_BAM,
                "-O", table.getAbsolutePath()
        };


        runCommandLine(args);
        final int yo = 30;
    }

    /***
     * Add these tests in the future to make sure that we count correctly. The data verified via manual IGV review.
     *
     * position, altDepth, altF1R2Depth
     * 20:21196089-21196089, 9, 5
     * 20:46952142-46952142, 1, 1
     * 20:53355622-53355622, 9, 2
     */
}