package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.annotator.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.testng.Assert.*;

/**
 * Created by emeryj on 8/11/17.
 */
public class AS_RMSMappingQualityUnitTest extends ReducibleAnnotationBaseTest {

    @Override
    protected List<String> getAnnotationsToUse() {
        return Collections.singletonList(AS_RMSMappingQuality.class.getSimpleName());
    }

    @Override
    protected String getRawKey() {
        return GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY;
    }

    @Test
    public void testFinalizeAnnotations() throws Exception {
        final List<String> annotationGroupsToUse = Collections.emptyList();
        final List<String> annotationsToUse = Collections.singletonList(AS_RMSMappingQuality.class.getSimpleName());
        final List<String> annotationsToExclude = Collections.emptyList();
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding, features);
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");
        final Genotype genotype = new GenotypeBuilder("sample2", Arrays.asList(refAllele, altAllele))
                .AD(new int[]{8,9}).make();

        final VariantContext vc = new VariantContextBuilder(ArtificialAnnotationUtils.makeVC()).attribute(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY, "285.00|385.00").genotypes(genotype).make();
        final VariantContext result = vae.finalizeAnnotations(vc, vc);
        Assert.assertNull(result.getAttribute(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY));
        Assert.assertNotNull(result.getAttribute(GATKVCFConstants.AS_RMS_MAPPING_QUALITY_KEY));
        Assert.assertEquals(1,1);
    }

}