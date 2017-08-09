package org.broadinstitute.hellbender.tools.walkers.orientationbias;

import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by tsato on 7/26/17.
 */

/***
 * This tools is the learning phase of the orientation filter.
 * Inference phase will likely feature variant context and what not.
 */
@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = VariantProgramGroup.class
)
public class ContextDependentArtifactFilter extends LocusWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "", common = false, optional = true)
    private File OUTPUT_FILE = null;

    @Argument(fullName = "", shortName = "", doc = "", optional = true)
    private File gnomad = null;


    @Argument(fullName = "", shortName = "", doc = "", optional = true)
    static final int DEFAULT_INITIAL_LIST_SIZE = 37_000_000/64; // by default we assume that that all 64 reference 3-mers are equally likely

    @Argument(fullName = "", shortName = "", doc = "", optional = true)
    static final int MINIMUM_MEDIAN_MQ_THRESHOLD = 20;

    private ContextDependentArtifactFilterEngine engine;

    public static Map<String, PerContextData> contextDependentDataMap;

    @Override
    public boolean requiresReference(){
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return ReadUtils.makeStandardReadFilters();
    }

    @Override
    public void onTraversalStart(){
        contextDependentDataMap = new HashMap<>();
        List<String> all3mers = SequenceUtil.generateAllKmers(3).stream().map(String::new).collect(Collectors.toList());

        for (final String refContext : all3mers){
            contextDependentDataMap.put(refContext,
                        new PerContextData(refContext));
        }
    }

    @Override
    public void apply(final AlignmentContext alignmentContext, final ReferenceContext referenceContext, final FeatureContext featureContext){
        // referenceContext always comes withe window of single base, so
        // manually expand the window and get the 3-mer for now.
        // TODO: implement getBasesInInterval() in referenceContext. Maybe simplify to getKmer(int k)?
        referenceContext.setWindow(1, 1);
        final String reference3mer = new String(referenceContext.getBases());
        assert reference3mer.length() == 3 : "kmer must have length 3";

        final ReadPileup pileup = alignmentContext.getBasePileup();
        final int[] baseCounts = pileup.getBaseCounts();

        // R in the docs
        final int depth = (int) MathUtils.sum(baseCounts);

        final byte refBase = reference3mer.getBytes()[1];

        /*** Enter Heuristic Land ***/

        // skip INDELs

        // skip MQ=0 loci
        List<Integer> mappingQualities = new ArrayList<>(pileup.size());

        // there is no shortcut or a standard API for converting an int[] to List<Integer> (we don't want List<int[]>)
        // so we must convert int[] to a List<Integer> with a for loop
        for (final int mq : pileup.getMappingQuals()) {
            mappingQualities.add(mq);
        }

        final double medianMQ = MathUtils.median(mappingQualities);

        if (medianMQ < MINIMUM_MEDIAN_MQ_THRESHOLD) {
            return;
        }

        final Optional<Byte> altBase = findAltBaseFromBaseCounts(baseCounts, refBase);
        final boolean isVariantSite = altBase.isPresent();

        /*** Exit Heuristic Land ***/

        // m in the docs
        final short altDepth = isVariantSite ? (short) baseCounts[BaseUtils.simpleBaseToBaseIndex(altBase.get())] : 0;

        // x in the docs
        final short altF1R2Depth = isVariantSite ? (short) pileup.getReads().stream().
                filter(r -> r.getBase(0) == altBase.get() && ReadUtils.isF2R1(r)).count() : 0;

        // FIXME: choose the correct allele
        final Allele allele = isVariantSite ? Allele.create(altBase.get(), false) : Allele.create(refBase, true);
        contextDependentDataMap.get(reference3mer).addNewSample(depth, altDepth, altF1R2Depth, allele, alignmentContext);
        return;
    }

    // FIXME: write tests
    private Optional<Byte> findAltBaseFromBaseCounts(final int[] baseCounts, final byte refBase) {
        final int[] baseCountsCopy = Arrays.copyOf(baseCounts, baseCounts.length);
        final long numObservedBases = Arrays.stream(baseCounts).filter(c -> c != 0).count();
        // FIXME: must handle hom var case when all the reads are alt
        if (numObservedBases == 1){
            return Optional.empty();
        }

        // now that we know there are multiple bases observed at the locus,
        // find the max out of the bases that are not ref
        // FIXME: also impose a minimum alt allele count and perhaps allele fraction (we're in the heuristic land anyway)
        baseCountsCopy[BaseUtils.simpleBaseToBaseIndex(refBase)] = 0;
        return Optional.of(BaseUtils.baseIndexToSimpleBase(MathUtils.argmax(baseCountsCopy)));
    }

    @Override
    public Object onTraversalSuccess() {
        engine = new ContextDependentArtifactFilterEngine();
        return null;
    }

    private List<String> makeAllPossible3Mers(){
        // TODO: this method needs to be improved
        final List<String> allPossible3Mers = new ArrayList<>(64); // 4^3 = 64
        final char[] possibleAlleles = "ACGT".toCharArray();
        for (int i = 0; i < possibleAlleles.length; i++){
            for (int j = 0; j < possibleAlleles.length; j++){
                for (int k = 0; k < possibleAlleles.length; k++){
                    allPossible3Mers.add(new StringBuilder().append(possibleAlleles[i]).append(possibleAlleles[j]).append(possibleAlleles[k])
                            .toString());
                }
            }
        }

        assert allPossible3Mers.size() == 64 : "there must be 64 kmers";
        return allPossible3Mers;
    }




}
