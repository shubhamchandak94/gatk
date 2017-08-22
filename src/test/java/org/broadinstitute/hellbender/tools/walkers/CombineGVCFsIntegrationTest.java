package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.vcf.VCFCodec;
import org.apache.commons.codec.digest.DigestUtils;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.broadinstitute.hellbender.utils.runtime.ProcessSettings;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;

public class CombineGVCFsIntegrationTest extends CommandLineProgramTest {
    private static final List<String> NO_EXTRA_ARGS = Collections.emptyList();
    private static final String BASE_PAIR_EXPECTED = "gvcf.basepairResolution.gatk3.7_30_ga4f720357.output.vcf";
    private static final String b38_reference_20_21 = largeFileTestDir + "Homo_sapiens_assembly38.20.21.fasta";
    private static final String BASE_PAIR_GVCF = "gvcf.basepairResolution.gvcf";

    private static final File CEUTRIO_20_21_GATK3_4_G_VCF = new File(largeFileTestDir, "gvcfs/CEUTrio.20.21.gatk3.4.g.vcf");
    private static final String CEUTRIO_20_21_EXPECTED_VCF = "CEUTrio.20.21.gatk3.7_30_ga4f720357.expected.vcf";

    private static <T> void assertForEachElementInLists(final List<T> actual, final List<T> expected, final BiConsumer<T, T> assertion) {
        Assert.assertEquals(actual.size(), expected.size(), "different number of elements in lists:\n"
                + actual.stream().map(Object::toString).collect(Collectors.joining("\n","actual:\n","\n"))
                +  expected.stream().map(Object::toString).collect(Collectors.joining("\n","expected:\n","\n")));
        for (int i = 0; i < actual.size(); i++) {
            assertion.accept(actual.get(i), expected.get(i));
        }
    }

    @DataProvider(name = "gvcfsToGenotype")
    public Object[][] gvcfsToGenotype() {
        return new Object[][]{
                //combine not supported yet, see https://github.com/broadinstitute/gatk/issues/2429 and https://github.com/broadinstitute/gatk/issues/2584
                //{"combine.single.sample.pipeline.1.vcf", null, Arrays.asList("-V", getTestFile("combine.single.sample.pipeline.2.vcf").toString() , "-V", getTestFile("combine.single.sample.pipeline.3.vcf").toString()), b37_reference_20_21},
                {new File[]{getTestFile("leadingDeletion.g.vcf")}, getTestFile("leadingDeletionRestrictToStartExpected.vcf"), Arrays.asList("-L", "20:69512-69513"), b37_reference_20_21}
        };
    }


    /*
This test is useful for testing changes in GATK4 versus different versions of GATK3.
To use, set GATK3_PATH to point to a particular version of gatk, then enable this test and run.

It will cache the gatk3 outputs in a folder called gatk3results, it does it's best to avoid reusing bad results by
comparing the md5 of the gatk3 path, input file path, reference, and commandline, but it doesn't know about internal changes to files
You have to manually delete the cache if you make changes inside the input files.

The expected outputs based on gatk3's results will be put in a folder called expectedResults.
These are overwritten during each test, the names are based on the name of the existing expected output file.

This method should be removed after GenotypeGVCFs has been completely validated against GATK3.
 */
    @Test(dataProvider = "gvcfsToGenotype")
    public void compareToGATK3(File[] inputs, File outputFile, List<String> extraArgs, String reference) throws IOException, NoSuchAlgorithmException {
        final String GATK3_PATH = "/Users/emeryj/hellbender/gsa-unstable/target/package/GenomeAnalysisTK.jar";
        final String params = GATK3_PATH + inputs[0].getAbsolutePath() + extraArgs.stream().collect(Collectors.joining()) + reference;
        final String md5 = DigestUtils.md5Hex(params);
        final File gatk3ResultsDir = new File("gatk3results");
        if(! gatk3ResultsDir.exists()){
            Assert.assertTrue(gatk3ResultsDir.mkdir());
        }
        final File gatk3Result = new File(gatk3ResultsDir, md5 + ".vcf");
        if (!gatk3Result.exists()) {
            List<String> gatk3Command = new ArrayList<>(
                    Arrays.asList("java", "-jar", GATK3_PATH, "-T", "CombineGVCFs"));
            for (File f: inputs) {
                gatk3Command.add("-V");
                gatk3Command.add(f.getAbsolutePath());
            }
            gatk3Command.add("-o");
            gatk3Command.add(gatk3Result.getAbsolutePath());
            gatk3Command.add("-R");
            gatk3Command.add(reference);
            gatk3Command.addAll(extraArgs);

            runProcess(new ProcessController(), gatk3Command.toArray(new String[gatk3Command.size()]));
        } else {
            System.out.println("Found precomputed gatk3Result");
        }
        final Path expectedResultsDir = Paths.get("expectedResults");
        if ( !Files.exists(expectedResultsDir)) {
            Files.createDirectory(expectedResultsDir);
        }
        Files.copy(gatk3Result.toPath(), expectedResultsDir.resolve(outputFile.getName()), StandardCopyOption.REPLACE_EXISTING);

        assertGenotypesMatch(Arrays.asList(inputs), gatk3Result, extraArgs, reference);
        assertVariantContextsMatch(Arrays.asList(inputs), gatk3Result, extraArgs, reference);
    }

    private static void runProcess(ProcessController processController, String[] command) {
        final ProcessSettings prs = new ProcessSettings(command);
        prs.getStderrSettings().printStandard(true);
        prs.getStdoutSettings().printStandard(true);
        final ProcessOutput output = processController.exec(prs);
        Assert.assertEquals(output.getExitValue(), 0, "Process exited with non-zero value. Command: "+ Arrays.toString(command) + "\n");
    }



    private void assertVariantContextsMatch(List<File> inputs, File expected, List<String> extraArgs, String reference) throws IOException {
        runCombineGVCFSandAssertSomething(inputs, expected, extraArgs, (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a, e, Collections.emptyList()), reference);
    }

    private void assertGenotypesMatch(List<File> inputs, File expected, List<String> additionalArguments, String reference) throws IOException {
        runCombineGVCFSandAssertSomething(inputs, expected, additionalArguments, VariantContextTestUtils::assertVariantContextsHaveSameGenotypes,
                reference);
    }

    private void runCombineGVCFSandAssertSomething(List<File> inputs, File expected, List<String> additionalArguments, BiConsumer<VariantContext, VariantContext> assertion, String reference) throws IOException {
        final File output = createTempFile("genotypegvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(reference))
                .addOutput(output);
        for (File input: inputs) {
            args.addArgument("V", input.getAbsolutePath());
        }

        additionalArguments.forEach(args::add);

        Utils.resetRandomGenerator();
        runCommandLine(args);

        final List<VariantContext> expectedVC = getVariantContexts(expected);
        final List<VariantContext> actualVC = getVariantContexts(output);
        assertForEachElementInLists(actualVC, expectedVC, assertion);
    }

    /**
     * Returns a list of VariantContext records from a VCF file
     *
     * @param vcfFile VCF file
     * @return list of VariantContext records
     * @throws IOException if the file does not exist or can not be opened
     */
    private static List<VariantContext> getVariantContexts(final File vcfFile) throws IOException {
        final VCFCodec codec = new VCFCodec();
        final FileInputStream s = new FileInputStream(vcfFile);
        final LineIterator lineIteratorVCF = codec.makeSourceFromStream(new PositionalBufferedStream(s));
        codec.readHeader(lineIteratorVCF);

        final List<VariantContext> VCs = new ArrayList<>();
        while (lineIteratorVCF.hasNext()) {
            final String line = lineIteratorVCF.next();
            Assert.assertFalse(line == null);
            VCs.add(codec.decode(line));
        }

        return VCs;
    }

//    @Test
//    public void testOneStartsBeforeTwoAndEndsAfterwards() throws Exception {
//        final String cmd = baseTestString(" -L 1:69485-69509");
//        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList(""));
//        spec.disableShadowBCF();
//        final File gVCF = executeTest("testOneStartsBeforeTwoAndEndsAfterwards", spec).first.get(0);
//        final List<VariantContext> allVCs = GATKVCFUtils.readVCF(gVCF).getSecond();
//
//        Assert.assertEquals(allVCs.size(), 2, "Observed: " + allVCs);
//
//        final VariantContext first = allVCs.get(0);
//        Assert.assertEquals(first.getStart(), 69491);
//        Assert.assertEquals(first.getEnd(), 69497);
//        Assert.assertEquals(first.getGenotypes().size(), 2);
//        Assert.assertTrue(first.getGenotype("NA1").isNoCall());
//        Assert.assertTrue(first.getGenotype("NA2").isNoCall());
//
//        final VariantContext second = allVCs.get(1);
//        Assert.assertEquals(second.getStart(), 69498);
//        Assert.assertEquals(second.getEnd(), 69506);
//        Assert.assertEquals(second.getGenotypes().size(), 2);
//        Assert.assertTrue(second.getGenotype("NA1").isNoCall());
//        Assert.assertTrue(second.getGenotype("NA2").isNoCall());
//    }
//
//    @Test(enabled = true)
//    public void testTetraploidRun() {
//        WalkerTestSpec spec = new WalkerTestSpec(
//                "-T CombineGVCFs -R " + b37KGReference + " -o %s --no_cmdline_in_header -V:sample1 " + privateTestDir + "tetraploid-gvcf-1.vcf" +
//                        " -V:sample2 " + privateTestDir + "tetraploid-gvcf-2.vcf" +
//                        " -V:sample3 " + privateTestDir + "tetraploid-gvcf-3.vcf" +
//                        " -L " + privateTestDir + "tetraploid-gvcfs.intervals",
//                1,
//                Arrays.asList("8472fcd43a41e7501b4709f0f1f5f432"));
//        executeTest("combineSingleSamplePipelineGVCF", spec);
//    }
//
//    @Test(enabled= true)
//    public void testMixedPloidyRun() {
//        WalkerTestSpec spec = new WalkerTestSpec(
//                "-T CombineGVCFs -R " + b37KGReference + " -o %s --no_cmdline_in_header -V:sample1 " + privateTestDir + "haploid-gvcf-1.vcf" +
//                        " -V:sample2 " + privateTestDir + "tetraploid-gvcf-2.vcf" +
//                        " -V:sample3 " + privateTestDir + "diploid-gvcf-3.vcf" +
//                        " -L " + privateTestDir + "tetraploid-gvcfs.intervals",
//                1,
//                Arrays.asList("f88b2b4d285276130cf088e7a03ca6a7"));
//        executeTest("combineSingleSamplePipelineGVCF", spec);
//    }
//
//    @Test
//    public void testTwoSpansManyBlocksInOne() throws Exception {
//        final String cmd = baseTestString(" -L 1:69512-69634");
//        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList(""));
//        spec.disableShadowBCF();
//        final File gVCF = executeTest("testTwoSpansManyBlocksInOne", spec).first.get(0);
//        final List<VariantContext> allVCs = GATKVCFUtils.readVCF(gVCF).getSecond();
//
//        Assert.assertEquals(allVCs.size(), 5);
//    }
//
//    @Test
//    public void testOneHasAltAndTwoHasNothing() throws Exception {
//        final String cmd = baseTestString(" -L 1:69511");
//        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList(""));
//        spec.disableShadowBCF();
//        final File gVCF = executeTest("testOneHasAltAndTwoHasNothing", spec).first.get(0);
//        final List<VariantContext> allVCs = GATKVCFUtils.readVCF(gVCF).getSecond();
//
//        Assert.assertEquals(allVCs.size(), 1);
//
//        final VariantContext first = allVCs.get(0);
//        Assert.assertEquals(first.getStart(), 69511);
//        Assert.assertEquals(first.getEnd(), 69511);
//        Assert.assertEquals(first.getGenotypes().size(), 2);
//    }
//
//    @Test
//    public void testOneHasAltAndTwoHasRefBlock() throws Exception {
//        final String cmd = baseTestString(" -L 1:69635");
//        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList(""));
//        spec.disableShadowBCF();
//        final File gVCF = executeTest("testOneHasAltAndTwoHasRefBlock", spec).first.get(0);
//        final List<VariantContext> allVCs = GATKVCFUtils.readVCF(gVCF).getSecond();
//
//        Assert.assertEquals(allVCs.size(), 1);
//
//        final VariantContext first = allVCs.get(0);
//        Assert.assertEquals(first.getStart(), 69635);
//        Assert.assertEquals(first.getEnd(), 69635);
//        Assert.assertEquals(first.getNAlleles(), 3);
//        Assert.assertEquals(first.getGenotypes().size(), 2);
//    }
//
//    @Test
//    public void testOneHasDeletionAndTwoHasRefBlock() throws Exception {
//        final String cmd = baseTestString(" -L 1:69772-69783");
//        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList(""));
//        spec.disableShadowBCF();
//        final File gVCF = executeTest("testOneHasDeletionAndTwoHasRefBlock", spec).first.get(0);
//        final List<VariantContext> allVCs = GATKVCFUtils.readVCF(gVCF).getSecond();
//
//        Assert.assertEquals(allVCs.size(), 3);
//
//        final VariantContext first = allVCs.get(0);
//        Assert.assertEquals(first.getStart(), 69772);
//        Assert.assertEquals(first.getEnd(), 69776);
//        Assert.assertEquals(first.getNAlleles(), 3);
//        Assert.assertEquals(first.getGenotypes().size(), 2);
//
//        final VariantContext second = allVCs.get(1);
//        Assert.assertEquals(second.getStart(), 69773);
//        Assert.assertEquals(second.getEnd(), 69774);
//        Assert.assertEquals(second.getGenotypes().size(), 2);
//
//        final VariantContext third = allVCs.get(2);
//        Assert.assertEquals(third.getStart(), 69775);
//        Assert.assertEquals(third.getEnd(), 69783);
//        Assert.assertEquals(third.getGenotypes().size(), 2);
//    }
//
//    @Test
//    public void testMD5s() throws Exception {
//        final String cmd = baseTestString(" -L 1:69485-69791");
//        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList("e1a888e8116cf59d53ad919634a18e6c"));
//        spec.disableShadowBCF();
//        executeTest("testMD5s", spec);
//    }
//
//    @Test
//    public void testBasepairResolutionOutput() throws Exception {
//        final String cmd = baseTestString(" -L 1:69485-69791 --convertToBasePairResolution");
//        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList("e7e86722a49ad9730743c4952cdbedc7"));
//        spec.disableShadowBCF();
//        executeTest("testBasepairResolutionOutput", spec);
//    }
//
//    @Test
//    public void testBreakBlocks() throws Exception {
//        final String cmd = baseTestString(" -L 1:69485-69791 --breakBandsAtMultiplesOf 5");
//        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList("bd279625ccfb3adcd39d97f07f3a236e"));
//        spec.disableShadowBCF();
//        executeTest("testBreakBlocks", spec);
//    }
//
//    @Test
//    public void testSpanningDeletions() {
//        WalkerTestSpec spec = new WalkerTestSpec(
//                "-T CombineGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
//                        " -V " + privateTestDir + "spanningDel.1.g.vcf -V " + privateTestDir + "spanningDel.2.g.vcf",
//                1,
//                Arrays.asList("b22238e1ff584a157335429309fbfc5b"));
//        spec.disableShadowBCF();
//        executeTest("testSpanningDeletions", spec);
//    }
//
//    @Test
//    public void testMultipleSpanningDeletionsForOneSample() {
//        WalkerTestSpec spec = new WalkerTestSpec(
//                "-T CombineGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
//                        " -V " + privateTestDir + "spanningDel.many.g.vcf",
//                1,
//                Arrays.asList("b828ba5b69422ce32e234ea24d2df4c7"));
//        spec.disableShadowBCF();
//        executeTest("testMultipleSpanningDeletionsForOneSample", spec);
//    }
//
//    @Test
//    public void testMultipleSpanningDeletionsForOneSampleHaploid() {
//        WalkerTestSpec spec = new WalkerTestSpec(
//                "-T CombineGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
//                        " -V " + privateTestDir + "spanningDel.many.haploid.g.vcf",
//                1,
//                Arrays.asList("e707335ebd61bbe20775f76ad9b8c20d"));
//        spec.disableShadowBCF();
//        executeTest("testMultipleSpanningDeletionsForOneSampleHaploid", spec);
//    }
//
//    @Test
//    public void testMultipleSpanningDeletionsForOneSampleTetraploid() {
//        WalkerTestSpec spec = new WalkerTestSpec(
//                "-T CombineGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
//                        " -V " + privateTestDir + "spanningDel.many.tetraploid.g.vcf",
//                1,
//                Arrays.asList("d4c22bd32d136414bfd7a6ebc5152026"));
//        spec.disableShadowBCF();
//        executeTest("testMultipleSpanningDeletionsForOneSampleTetraploid", spec);
//    }
//
//    @Test
//    public void testWrongReferenceBaseBugFix() throws Exception {
//        final String cmd = "-T CombineGVCFs -R " + b37KGReference + " -V " + (privateTestDir + "combine-gvcf-wrong-ref-input1.vcf"
//                + " -V " + (privateTestDir + "combine-gvcf-wrong-ref-input2.vcf") + " -o %s --no_cmdline_in_header");
//        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList("5fec22c9c8a0063f43c86ac86bb12e27"));
//        spec.disableShadowBCF();
//        executeTest("testWrongReferenceBaseBugFix",spec);
//
//    }
//
//    @Test
//    public void testBasepairResolutionInput() throws Exception {
//        final String cmd = "-T CombineGVCFs -R " + b37KGReference + " -o %s --no_cmdline_in_header -V " + privateTestDir + "gvcf.basepairResolution.vcf";
//        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList("a839f81014758d1bdc900d59d35dd5bc"));
//        spec.disableShadowBCF();
//        executeTest("testBasepairResolutionInput", spec);
//    }
//
//    @Test
//    public void testAlleleSpecificAnnotations() throws Exception {
//        final String cmd = "-T CombineGVCFs -R " + b37KGReference + " -o %s --no_cmdline_in_header -G Standard -G AS_Standard -V "
//                + privateTestDir + "NA12878.AS.chr20snippet.g.vcf -V " + privateTestDir + "NA12892.AS.chr20snippet.g.vcf";
//        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList("26606fbfc3e9e5a8813307227d898b9a"));
//        spec.disableShadowBCF();
//        executeTest("testAlleleSpecificAnnotations", spec);
//    }
//
//    @Test
//    public void testMissingAlleleSpecificAnnotationGroup() throws IOException {
//        final File logFile = createTempFile("testMissingAlleleSpecificAnnotationGroup.log", ".tmp");
//        final String cmd = "-T CombineGVCFs -R " + b37KGReference + " -o %s --no_cmdline_in_header -V "
//                + privateTestDir + "NA12878.AS.chr20snippet.g.vcf -V " + privateTestDir + "NA12892.AS.chr20snippet.g.vcf -log " +
//                logFile.getAbsolutePath();
//        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList(""));
//        spec.disableShadowBCF();
//        executeTest("testMissingAlleleSpecificAnnotationGroup", spec);
//        Assert.assertTrue(FileUtils.readFileToString(logFile).contains(ReferenceConfidenceVariantContextMerger.ADD_AS_STANDARD_MSG));
//    }
//
//    @Test
//    public void testASMateRankSumAnnotation() throws Exception {
//        final String cmd = "-T CombineGVCFs -R " + b37KGReference + " -o %s --no_cmdline_in_header -G Standard -G AS_Standard -A AS_MQMateRankSumTest -V "
//                + privateTestDir + "NA12878.AS.MateRankSum.chr20snippet.g.vcf -V " + privateTestDir + "NA12892.AS.MateRankSum.chr20snippet.g.vcf";
//        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList("9fd72e0f1f0f08a5aff6875225eec567"));
//        spec.disableShadowBCF();
//        executeTest("testASMateRankSumAnnotation", spec);
//    }
//
//    @Test
//    public void testASInsertSizeRankSumAnnotation() throws Exception {
//        final String cmd = "-T CombineGVCFs -R " + b37KGReference + " -o %s --no_cmdline_in_header -G Standard -G AS_Standard -A AS_InsertSizeRankSum -V "
//                + privateTestDir + "NA12878.AS.InsertSizeRankSum.chr20snippet.g.vcf -V " + privateTestDir + "NA12892.AS.InsertSizeRankSum.chr20snippet.g.vcf";
//        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList("b8e10684010702d11ee71538b1f201ac"));
//        spec.disableShadowBCF();
//        executeTest("testASInsertSizeRankSumAnnotation", spec);
//    }


}
