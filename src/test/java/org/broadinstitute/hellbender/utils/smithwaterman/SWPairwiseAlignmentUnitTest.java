package org.broadinstitute.hellbender.utils.smithwaterman;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWAlignerArguments;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;

public final class SWPairwiseAlignmentUnitTest extends BaseTest {

    @DataProvider(name = "ComplexReadAlignedToRef")
    public Object[][] makeComplexReadAlignedToRef() {
        return new Object[][] {
                {"AAAGGACTGACTG", "ACTGACTGACTG", 1, "12M"}
        };
    }

    @Test(dataProvider = "ComplexReadAlignedToRef")
    public void testReadAlignedToRefComplexAlignment(final String reference, final String read, final int expectedStart, final String expectedCigar) {
        assertAlignmentMatchesExpected(reference, read, expectedStart, expectedCigar, SWPairwiseAlignment.ORIGINAL_DEFAULT,
                                       SWPairwiseAlignment.DEFAULT_OVERHANG_STRATEGY);
    }

    @DataProvider(name = "OddNoAlignment")
    public Object[][] makeOddNoAlignment() {
        final String ref1     = "AAAGACTACTG";
        final String read1    = "AACGGACACTG";
        return new Object[][] {
                {ref1, read1, new SWAlignerArguments.Weights( 50, -100, -220, -12), 1,  "2M2I3M1D4M"},
                {ref1, read1, new SWAlignerArguments.Weights(200, -50, -300, -22), 0, "11M"}
        };
    }

    @Test(dataProvider = "OddNoAlignment")
    public void testOddNoAlignment(final String reference, final String read, final SWAlignerArguments.Weights weights,
                                   final int expectedStart, final String expectedCigar) {
        assertAlignmentMatchesExpected(reference, read, expectedStart, expectedCigar, weights, SWPairwiseAlignment.DEFAULT_OVERHANG_STRATEGY);
    }

    @Test
    public void testIndelsAtStartAndEnd() {
        final String match     = "CCCCC";
        final String reference = "AAA" + match;
        final String read      = match + "GGG";
        final int expectedStart = 3;
        final String expectedCigar = "5M3S";
        assertAlignmentMatchesExpected(reference, read, expectedStart, expectedCigar, SWPairwiseAlignment.ORIGINAL_DEFAULT, SWPairwiseAlignment.DEFAULT_OVERHANG_STRATEGY);
    }

    @Test
    public void testDegenerateAlignmentWithIndelsAtBothEnds() {
        logger.warn("testDegenerateAlignmentWithIndelsAtBothEnds");
        final String ref = "TGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGA";
        final String alt =               "ACAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA";
        final int expectedStart = 14;
        final String expectedCigar = "31M20S";
        assertAlignmentMatchesExpected(ref, alt, expectedStart, expectedCigar, SWPairwiseAlignment.STANDARD_NGS,
                                       SWPairwiseAlignment.DEFAULT_OVERHANG_STRATEGY);
    }

    @Test
    public void  testForIdenticalAlignmentsWithDifferingFlankLengths() {
        //This test is designed to ensure that the indels are correctly placed
        //if the region flanking these indels is extended by a varying amount.
        //It checks for problems caused by floating point rounding leading to different
        //paths being selected.

        //Create two versions of the same sequence with different flanking regions.
        final byte[] paddedRef="GCGTCGCAGTCTTAAGGCCCCGCCTTTTCAGACAGCTTCCGCTGGGCCTGGGCCGCTGCGGGGCGGTCACGGCCCCTTTAAGCCTGAGCCCCGCCCCCTGGCTCCCCGCCCCCTCTTCTCCCCTCCCCCAAGCCAGCACCTGGTGCCCCGGCGGGTCGTGCGGCGCGGCGCTCCGCGGTGAGCGCCTGACCCCGAGGGGGCCCGGGGCCGCGTCCCTGGGCCCTCCCCACCCTTGCGGTGGCCTCGCGGGTCCCAGGGGCGGGGCTGGAGCGGCAGCAGGGCCGGGGAGATGGGCGGTGGGGAGCGCGGGAGGGACCGGGCCGAGCCGGGGGAAGGGCTCCGGTGACT".getBytes();
        final byte[] paddedHap="GCGTCGCAGTCTTAAGGCCCCGCCTTTTCAGACAGCTTCCGCTGGGCCTGGGCCGCTGCGGGGCGGTCACGGCCCCTTTAAGCCTGAGCCCCGCCCCCTGGCTCCCCGCCCCCTCTTCTCCCCTCCCCCAAGCCAGCACCTGGTGCCCCGGCGGGTCGTGCGGCGCGGCGCTCCGCGGTGAGCGCCTGACCCCGA--GGGCC---------------GGGCCCTCCCCACCCTTGCGGTGGCCTCGCGGGTCCCAGGGGCGGGGCTGGAGCGGCAGCAGGGCCGGGGAGATGGGCGGTGGGGAGCGCGGGAGGGACCGGGCCGAGCCGGGGGAAGGGCTCCGGTGACT".replace("-", "").getBytes();
        final byte[] notPaddedRef=                                                                           "CTTTAAGCCTGAGCCCCGCCCCCTGGCTCCCCGCCCCCTCTTCTCCCCTCCCCCAAGCCAGCACCTGGTGCCCCGGCGGGTCGTGCGGCGCGGCGCTCCGCGGTGAGCGCCTGACCCCGAGGGGGCCCGGGGCCGCGTCCCTGGGCCCTCCCCACCCTTGCGGTGGCCTCGCGGGTCCCAGGGGCGGGGCTGGAGCGGCAGCAGGGCCGGGGAGATGGGCGGTGGGGAGCGCGGGAGGGA".getBytes();
        final byte[] notPaddedHap=                                                                           "CTTTAAGCCTGAGCCCCGCCCCCTGGCTCCCCGCCCCCTCTTCTCCCCTCCCCCAAGCCAGCACCTGGTGCCCCGGCGGGTCGTGCGGCGCGGCGCTCCGCGGTGAGCGCCTGACCCCGA---------GGGCC--------GGGCCCTCCCCACCCTTGCGGTGGCCTCGCGGGTCCCAGGGGCGGGGCTGGAGCGGCAGCAGGGCCGGGGAGATGGGCGGTGGGGAGCGCGGGAGGGA".replace("-", "").getBytes();
        //a simplified version of the getCigar routine in the haplotype caller to align these
        final String SW_PAD = "NNNNNNNNNN";
        final String paddedsRef = SW_PAD + new String(paddedRef) + SW_PAD;
        final String paddedsHap = SW_PAD + new String(paddedHap) + SW_PAD;
        final String notPaddedsRef = SW_PAD + new String(notPaddedRef) + SW_PAD;
        final String notpaddedsHap = SW_PAD + new String(notPaddedHap) + SW_PAD;
        final SmithWatermanAligner aligner = new SWPairwiseAlignment(CigarUtils.NEW_SW_PARAMETERS, SWPairwiseAlignment.DEFAULT_OVERHANG_STRATEGY);
        final SmithWatermanAlignment paddedAlignment = aligner.align(paddedsRef.getBytes(), paddedsHap.getBytes());
        final SmithWatermanAlignment notPaddedAlignment = aligner.align(notPaddedsRef.getBytes(), notpaddedsHap.getBytes());
        //Now verify that the two sequences have the same alignment and not match positions.
        final Cigar rawPadded = paddedAlignment.getCigar();
        final Cigar notPadded= notPaddedAlignment.getCigar();
        final List<CigarElement> paddedC=rawPadded.getCigarElements();
        final List<CigarElement> notPaddedC=notPadded.getCigarElements();
        Assert.assertEquals(paddedC.size(), notPaddedC.size());
        for(int i=0;i<notPaddedC.size();i++)
        {
            final CigarElement pc=paddedC.get(i);
            final CigarElement npc=notPaddedC.get(i);
            if(pc.getOperator()== CigarOperator.M && npc.getOperator()== CigarOperator.M)
            {
                continue;
            }
            final int l1=pc.getLength();
            final int l2=npc.getLength();
            Assert.assertEquals(l1, l2);
            Assert.assertEquals(pc.getOperator(), npc.getOperator());
        }
    }

    @DataProvider
    public Object[][] getSubstringMatchTests(){
        return new Object[][]{
                {3, "5M", SWAlignerArguments.OverhangStrategy.SOFTCLIP},
                {0, "3D5M", SWAlignerArguments.OverhangStrategy.INDEL},
                {0, "3D5M", SWAlignerArguments.OverhangStrategy.LEADING_INDEL},
                {3, "5M", SWAlignerArguments.OverhangStrategy.IGNORE}
        };
    }

    @Test(dataProvider = "getSubstringMatchTests")
    public void testSubstringMatch(int expectedStart, String expectedCigar, SWAlignerArguments.OverhangStrategy strategy) {
        final String matchingSection = "CCCCC";
        final String reference = "AAA" + matchingSection;
        final String read = matchingSection;
        assertAlignmentMatchesExpected(reference, read, expectedStart, expectedCigar, SWPairwiseAlignment.ORIGINAL_DEFAULT,
                                       strategy);
    }

    private static void assertAlignmentMatchesExpected(String reference, String read, int expectedStart, String expectedCigar, SWAlignerArguments.Weights weights, SWAlignerArguments.OverhangStrategy strategy) {
        final SWPairwiseAlignment sw = new SWPairwiseAlignment(weights, strategy);
        final SmithWatermanAlignment alignment = sw.align(reference.getBytes(), read.getBytes());
        sw.printAlignment(reference.getBytes(), read.getBytes(), alignment);
        Assert.assertEquals(alignment.getAlignmentOffset(), expectedStart);
        Assert.assertEquals(alignment.getCigar().toString(), expectedCigar);
    }

    @DataProvider
    public Object[][] getSubstringMatchLong(){
        return new Object[][]{
                {359, "7M", SWAlignerArguments.OverhangStrategy.SOFTCLIP},
                {0, "1M358D6M29D", SWAlignerArguments.OverhangStrategy.INDEL},
                {0, "1M1D6M", SWAlignerArguments.OverhangStrategy.LEADING_INDEL},
                {359, "7M", SWAlignerArguments.OverhangStrategy.IGNORE}
        };
    }

    @Test(dataProvider = "getSubstringMatchLong")
    public void testSubstringMatchLong(int expectedStart, String expectedCigar, SWAlignerArguments.OverhangStrategy strategy) {
        final String reference = "ATAGAAAATAGTTTTTGGAAATATGGGTGAAGAGACATCTCCTCTTATGGAAAAAGGGATTCTAGAATTTAACAATAAATATTCCCAACTTTCCCCAAGGCTTTAAAATCTACCTTGAAGGAGCAGCTGATGTATTTCTAGAACAGACTTAGGTGTCTTGGTGTGGCCTGTAAAGAGATACTGTCTTTCTCTTTTGAGTGTAAGAGAGAAAGGACAGTCTACTCAATAAAGAGTGCTGGGAAAACTGAATATCCACACACAGAATAATAAAACTAGATCCTATCTCTCACCATATACAAAGATCAACTCAAAACAAATTAAAGACCTAAATGTAAGACAAGAAATTATAAAACTACTAGAAAAAAACACAAGGGAAATGCTTCAGGACATTGGC";
        final String read      = "AAAAAAA";
        assertAlignmentMatchesExpected(reference, read, expectedStart, expectedCigar, SWPairwiseAlignment.ORIGINAL_DEFAULT, strategy);
    }

}
