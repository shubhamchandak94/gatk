package org.broadinstitute.hellbender.tools.copynumber.legacy;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.junit.Test;

import static org.testng.Assert.*;

/**
 * Created by slee on 9/6/17.
 */
public class ModelSegmentsIntegrationTest extends CommandLineProgramTest {
    @Test
    public void testTCGA() {
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-01A-01D-1265-08-gc-corrected.tn.tsv",
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-01A-01D-1265-08-gc-corrected.seg",
                "-" + StandardArgumentDefinitions.VERBOSITY_NAME, "DEBUG"
        };
        runCommandLine(arguments);
    }
}