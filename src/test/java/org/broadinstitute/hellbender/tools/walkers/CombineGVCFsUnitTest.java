package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.junit.Test;
import org.mockito.ArgumentCaptor;
import org.mockito.Mockito;

import java.util.Arrays;
import java.util.Collections;

import static org.mockito.Mockito.when;
import static org.testng.Assert.*;

/**
 * Created by emeryj on 8/25/17.
 */
public class CombineGVCFsUnitTest {
    public static final Allele REF = Allele.create("A", true);
    public static final Allele ALT = Allele.create("C", false);

    //TODO unit test for gvcf blockign

    //TODO unit test for how the end previous states works?

    //TODO fill this out
    @Test
    public void testClearAccumulatedReads() {
//        CombineGVCFs combinegvcfs = Mockito.mock(CombineGVCFs.class);
//        VariantContext[] VCs = new VariantContext[]{
//        new VariantContextBuilder("test", "1", 100, 105, Collections.singleton(REF)).make(),
//        new VariantContextBuilder("test", "1", 100, 110, Arrays.asList(REF, ALT)).make(),
//        new VariantContextBuilder("test", "1", 100, 110, Arrays.asList(REF, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE)).make()};
//        CombineGVCFs.OverallState state = new CombineGVCFs.OverallState();
//        state.VCs = Arrays.asList( VCs);
//        when(combinegvcfs.currentOverallState = new CombineGVCFs.OverallState());
    }

}