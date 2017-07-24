package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;

import static org.testng.Assert.*;

/**
 * Created by emeryj on 8/11/17.
 */
public class AS_MappingQualityRankSumTestUnitTest extends ReducibleAnnotationBaseTest {
    @Override
    protected List<String> getAnnotationsToUse() {
        return Collections.singletonList(AS_MappingQualityRankSumTest.class.getSimpleName());
    }

    @Override
    protected String getRawKey() {
        return GATKVCFConstants.AS_RAW_MAP_QUAL_RANK_SUM_KEY;
    }
}