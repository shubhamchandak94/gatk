package org.broadinstitute.hellbender.cmdline;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;


public final class PicardProgramExecutorIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DATA_PATH = new File(getTestDataDir(), "picard/fastq/NormalizeFasta");

    @Override
    public String getTestedClassName() {
        return picard.reference.NormalizeFasta.class.getSimpleName();
    }

    @Test
    public void testPicardNormalizeFasta() throws IOException {
        final File input = new File(TEST_DATA_PATH, "testfasta.fasta");
        final File expectedFile = new File(TEST_DATA_PATH, "testFASTA_WS_4.fasta");
        final File outfile = createTempFile("normalized", ".fasta");

        final String[] args = {
                "-I", input.getAbsolutePath(),
                "-O", outfile.getAbsolutePath(),
                "--TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE", "TRUE",
                "--LINE_LENGTH", "5",
        };

        Assert.assertEquals(runCommandLine(args), 0);

        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile);
    }

    @Test
    public void testPicardNormalizeFastaWithBadArgs() throws IOException {
        final File input = new File(TEST_DATA_PATH, "testfasta.fasta");
        final File outfile = createTempFile("normalized", ".fasta");

        // use GATK-style lower case names, which are rejected by Picard
        final String[] args = {
                "--input", input.getAbsolutePath(),
                "--output", outfile.getAbsolutePath(),
                "--TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE", "TRUE",
                "--LINE_LENGTH", "5",
        };

        Assert.assertNotEquals(runCommandLine(args), 0);
    }

}
