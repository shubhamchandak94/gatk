package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by tsato on 7/26/17.
 */

/***
 * This tools is the learning phase of the orientation filter.
 * Inference phase will likely feature variant context and what not.
 */
@CommandLineProgramProperties(
        summary = "yoooo",
        oneLineSummary = "hooooo",
        programGroup = VariantProgramGroup.class
)

public class ContextDependentArtifactFilter extends LocusWalker {
    @Argument(fullName = "", shortName = "", doc = "", optional = true)
    private File gnomad = null;

    @Argument(fullName = "", shortName = "", doc = "", optional = true)
    static final int DEFAULT_INITIAL_LIST_SIZE = 37_000_000/64; // by default we assume that that all 64 reference 3-mers are equally likely

    @Argument(fullName = "", shortName = "", doc = "exclude reads below this quality from pileup", optional = true)
    static final int MINIMUM_MEDIAN_MQ = 20;

    @Argument(fullName = "", shortName = "", doc = "exclude bases below this quality from pileup", optional = true)
    static final int MINIMUM_BASE_QUALITY = 10;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "a tab-seprated table of hyperparameters", optional = false)
    static File output = null;

    public static Map<String, PerContextData> contextDependentDataMap;

    public static final List<String> ALL_3_MERS = SequenceUtil.generateAllKmers(3).stream().map(String::new).collect(Collectors.toList());

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

        for (final String refContext : ALL_3_MERS){
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

        final ReadPileup pileup = alignmentContext.filter(pe -> pe.getQual() > MINIMUM_BASE_QUALITY).getBasePileup();

        // FIXME; this is not ideal. AlignmentContext should come filtered and not reach here if it's empty
        if (pileup.size() == 0){
            return;
        }

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

        if (medianMQ < MINIMUM_MEDIAN_MQ) {
            return;
        }

        final Optional<Byte> altBase = findAltBaseFromBaseCounts(baseCounts, refBase);
        final boolean isVariantSite = altBase.isPresent();

        /*** Exit Heuristic Land ***/

        // m in the docs
        final short altDepth = isVariantSite ? (short) baseCounts[BaseUtils.simpleBaseToBaseIndex(altBase.get())] : 0;

        // x in the docs
        final short altF1R2Depth = isVariantSite ? (short) pileup.getNumberOfElements(pe -> pe.getBase() == altBase.get() && ! ReadUtils.isF2R1(pe.getRead())) : 0;
        assert altDepth >= altF1R2Depth : String.format("altDepth >= altF1R2Depth but got %d, %d", altDepth, altF1R2Depth);


        // FIXME: choose the correct allele
        final Allele allele = isVariantSite ? Allele.create(altBase.get(), false) : Allele.create(refBase, true);
        contextDependentDataMap.get(reference3mer).addNewExample(depth, altDepth, altF1R2Depth, allele, alignmentContext);
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
        List<Hyperparameters> hyperparameterEstimates = new ArrayList<>();
        // remember we run EM separately for each of 4^3 = 64 ref contexts


        // debug
        for (final String refContext : ALL_3_MERS){
            PerContextData contextData = contextDependentDataMap.get(refContext);
            ContextDependentArtifactFilterEngine engine = new ContextDependentArtifactFilterEngine(contextData);
            final int numSites = contextData.getNumLoci();
            Hyperparameters hyperparameters = engine.runEMAlgorithm();
            hyperparameterEstimates.add(hyperparameters);


//            List<Integer> activeIndices = IntStream.range(0, numSites).filter(i -> contextData.getAltDepths().get(i) > 0)
//                    .boxed().collect(Collectors.toList());
//            List<Short> activeAltDepths = activeIndices.stream().map(i -> contextData.altDepths.get(i)).collect(Collectors.toList());
//            List<Short> activeAltF1R2Depths = activeIndices.stream().map(i -> contextData.altF1R2Depths.get(i)).collect(Collectors.toList());
//            List<String> activePositions= activeIndices.stream().map(i -> contextData.positions.get(i)).collect(Collectors.toList());
//            int bogus = 3;
        }

        Hyperparameters.writeHyperparameters(hyperparameterEstimates, output);

        /**
         *  internal state check - must go to a separate test at some point
         *  20:21196089-21196089, 9, 5 // ACG;
         *  20:46952142-46952142, 1, 1 // ACG
         *  20:53355622-53355622, 9, 2 // AAC
         */

        return "SUCCESS";
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
