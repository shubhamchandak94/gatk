package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.walkers.annotator.ReadPosRankSumTest;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

public class AS_ReadPosRankSumTestUnitTest extends ReducibleAnnotationBaseTest {

    @DataProvider(name = "dataIsUsableRead")
    private Object[][] dataIsUsableRead(){
        return new Object[][]{
                {"20M6D2M", 10, 1, 0, true},
                {"20M6D2M", 10, 1, 27, true},
                {"20M6D2M", 10, 1, 29, false},
                {"1I20M1S", 10, 1, 0, true},
                {"1I20M1S", 0, 1, 0, false},
                {"1I20M1S", QualityUtils.MAPPING_QUALITY_UNAVAILABLE, 1, 0, false},
                {"1I20M1S", 10, 1, 22, false},
                {"1I20M1S", 10, 21, 42, false},
                {"1I20M1H", 10, 1, 21, false},
        };
    }

    @Test(dataProvider = "dataIsUsableRead")
    public void testIsUsableRead(final String cigarString, final int mappingQuality, final int start, final int refLoc, final boolean isUsable ) {
        final AS_ReadPosRankSumTest as_readPosRankSumTest = new AS_ReadPosRankSumTest();
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(cigar);
        read.setMappingQuality(mappingQuality);
        read.setPosition("1", start);
        Assert.assertEquals(as_readPosRankSumTest.isUsableRead(read, refLoc), isUsable);
    }

    @Override
    protected List<String> getAnnotationsToUse() {
        return Collections.singletonList(AS_ReadPosRankSumTest.class.getSimpleName());
    }

    @Override
    protected String getRawKey() {
        return GATKVCFConstants.AS_RAW_READ_POS_RANK_SUM_KEY;
    }

}
