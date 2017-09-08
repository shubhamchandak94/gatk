package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;

import java.io.File;

/**
 *
 * The inference phase of the orientation bias filter
 * Created by tsato on 8/8/17.
 */
@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = VariantProgramGroup.class
)
public class ContextDependentArtifactFilterInferenceEngine extends VariantWalker {

    @Override
    public void onTraversalStart(){
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext ){

    }

    @Override
    public Object onTraversalSuccess(){
        return null;
    }

}
