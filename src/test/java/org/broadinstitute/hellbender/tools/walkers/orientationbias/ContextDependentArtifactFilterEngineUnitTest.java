package org.broadinstitute.hellbender.tools.walkers.orientationbias;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.broadinstitute.hellbender.tools.walkers.orientationbias.ContextDependentArtifactFilterEngine.States;

import java.util.Arrays;

/**
 * Created by tsato on 8/14/17.
 */
public class ContextDependentArtifactFilterEngineUnitTest {
    private final byte A = "A".getBytes()[0];
    private final byte C = "A".getBytes()[0];
    private final byte G = "G".getBytes()[0];
    private final byte T = "T".getBytes()[0];

    /**
     * Create a test case: 100 sites (per context). Oh, perhaps actually generate the data from appropriate distributions and
     * see if we catch the artifact
     *
     * 5 alt sites per transition. But, G -> T are all orientation biased, or something like that. Good, write this test tomorrow.
     */
    @Test
    public void testCatchBasicOrientationBias() {
        final int NUM_EXAMPLES = 100;
        final int QUOTIENT = 5;
        final int NUM_F1R2_EXAMPLES = 20;

        final int DEPTH = 100;
        final short STANTARD_ALT_DEPTH = 20;
        final short BALANCED_ALT_F1R2_DEPTH = STANTARD_ALT_DEPTH/2;
        final short ENTIRELY_F1R2 = STANTARD_ALT_DEPTH;
        final short ENTIRELY_F2R1 = STANTARD_ALT_DEPTH;

        final String refContext = "AGT";
        final Allele refAllele = Allele.create("G", true);
        final Allele altAllele = Allele.create("T", false);

        final PerContextData data = new PerContextData(refContext);
        for (int n = 0; n < NUM_EXAMPLES; n++){
            // assume 20 in 100 examples is alt, and 20% allele fraction, and alt reads are entirely F1R2
            final short altDepth = n % QUOTIENT == 0 ? STANTARD_ALT_DEPTH : 0;
            final Allele allele = n % QUOTIENT == 0 ? altAllele : refAllele;
            final short altF1R2Depth = n % QUOTIENT == 0 ? ENTIRELY_F1R2 : BALANCED_ALT_F1R2_DEPTH;

            final SimpleInterval snpLocation = new SimpleInterval("20", 1_000_000 + n, 1_000_000 + n);
            final AlignmentContext alignmentContext = new AlignmentContext(snpLocation, new ReadPileup(snpLocation));

            data.addNewExample(DEPTH, altDepth, altF1R2Depth, allele, alignmentContext);
        }

        ContextDependentArtifactFilterEngine engine = new ContextDependentArtifactFilterEngine(data);
        ContextDependentArtifactFilterEngine.Hyperparameters hyperparameters = engine.runEMAlgorithm();

        // make an enum for the states of Z...but how do you code that?
        for (int a = 0; a < ContextDependentArtifactFilterEngine.NUM_ALLELES; a++){
            // TODO: how should we code for the case when there are 0 alts of that type? Can pi for that allele be 0?
            // Assert.assertEquals(MathUtils.sum(hyperparameters.getPi()[a]), 1.0);
        }

        System.out.println(Arrays.toString(engine.effectiveCounts));
        Assert.assertEquals(engine.effectiveCounts[States.F1R2.ordinal()], (double) NUM_F1R2_EXAMPLES);
        // START HERE, we get something like 78 - but we want 80. For such a simple example, we gotta get this right
        Assert.assertEquals(engine.effectiveCounts[States.BALANCED_HOM_REF.ordinal()], (double) NUM_EXAMPLES - NUM_F1R2_EXAMPLES);

        // Given that the allele is G, all fo the examples were home ref, so p(z = hom ref|a = G) must equal 1.0
        Assert.assertEquals(hyperparameters.getPi()[BaseUtils.simpleBaseToBaseIndex(G)][States.BALANCED_HOM_REF.ordinal()], 1.0);

        // Given that the allele is T, all of the examples were F1R2 artifacts, so we want p(z = F1R2|a = T) = 1.0 and
        // p(z = F2R1|a = T) = 0.0
        Assert.assertEquals(hyperparameters.getPi()[BaseUtils.simpleBaseToBaseIndex(T)][States.F1R2.ordinal()], 1.0);
        Assert.assertEquals(hyperparameters.getPi()[BaseUtils.simpleBaseToBaseIndex(T)][States.F2R1.ordinal()], 0.0);

        // We expect the model to learn the correct allele fraction given z = F1R2
        Assert.assertEquals(hyperparameters.getF()[States.F1R2.ordinal()], (double) ENTIRELY_F1R2/DEPTH);
    }
}