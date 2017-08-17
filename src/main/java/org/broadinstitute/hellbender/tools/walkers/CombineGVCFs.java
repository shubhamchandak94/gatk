package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.GenotypeUtils;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VcfUtils;

import java.util.*;

/**
 * Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file
 *
 * <p>
 * CombineGVCFs is meant to be used for hierarchical merging of gVCFs that will eventually be input into GenotypeGVCFs.
 * One would use this tool when needing to genotype too large a number of individual gVCFs; instead of passing them
 * all in to GenotypeGVCFs, one would first use CombineGVCFs on smaller batches of samples and then pass these combined
 * gVCFs to GenotypeGVCFs.</p>
 *
 * <h3>Input</h3>
 * <p>
 * Two or more Haplotype Caller gVCFs to combine.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A combined multisample gVCF.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T CombineGVCFs \
 *   -R reference.fasta \
 *   --variant sample1.g.vcf \
 *   --variant sample2.g.vcf \
 *   -o cohort.g.vcf
 * </pre>
 *
 * <h3>Caveats</h3>
 * <p>Only gVCF files produced by HaplotypeCaller (or CombineGVCFs) can be used as input for this tool. Some other
 * programs produce files that they call gVCFs but those lack some important information (accurate genotype likelihoods
 * for every position) that GenotypeGVCFs requires for its operation.</p>
 * <p>If the gVCF files contain allele specific annotations, add -G Standard -G AS_Standard to the command line.</p>
 *
 */
@CommandLineProgramProperties(summary = "TODO", oneLineSummary = "TODO", programGroup = VariantProgramGroup.class)
@DocumentedFeature
public final class CombineGVCFs extends MultiVariantWalker {

    /**
     * Which annotations to recompute for the combined output VCF file.
     */
    @Advanced
    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to recompute.  The single value 'none' removes the default annotations", optional=true)
    protected List<String> annotationsToUse = new ArrayList<>(Arrays.asList(new String[]{"AS_RMSMappingQuality"}));

    /**
     * Which groups of annotations to add to the output VCF file. The single value 'none' removes the default group. See
     * the VariantAnnotator -list argument to view available groups. Note that this usage is not recommended because
     * it obscures the specific requirements of individual annotations. Any requirements that are not met (e.g. failing
     * to provide a pedigree file for a pedigree-based annotation) may cause the run to fail.
     */
    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", optional=true)
    protected String[] annotationGroupsToUse = { StandardAnnotation.class.getSimpleName() };


    /**
     * The rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate. Note that dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();
    public FeatureInput<VariantContext> getDbsnpRodBinding() { return dbsnp.dbsnp; }
    public List<RodBinding<VariantContext>> getCompRodBindings() { return Collections.emptyList(); }
    public RodBinding<VariantContext> getSnpEffRodBinding() { return null; }
    public List<RodBinding<VariantContext>> getResourceRodBindings() { return Collections.emptyList(); }
    public boolean alwaysAppendDbsnpId() { return false; }

    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {

    }

    protected final class PositionalState {
        final List<VariantContext> VCs;
        final Set<String> samples = new HashSet<>();
        final byte[] refBases;
        final GenomeLoc loc;
        public PositionalState(final List<VariantContext> VCs, final byte[] refBases, final GenomeLoc loc) {
            this.VCs = VCs;
            for(final VariantContext vc : VCs){
                samples.addAll(vc.getSampleNames());
            }
            this.refBases = refBases;
            this.loc = loc;
        }
    }

    protected final class OverallState {
        final LinkedList<VariantContext> VCs = new LinkedList<>();
        final Set<String> samples = new HashSet<>();
        GenomeLoc prevPos = null;
        byte refAfterPrevPos;

        public OverallState() {}
    }

    @Argument(fullName="convertToBasePairResolution", shortName="bpResolution", doc = "If specified, convert banded gVCFs to all-sites gVCFs", required=false)
    protected boolean USE_BP_RESOLUTION = false;

    /**
     * To reduce file sizes our gVCFs group similar reference positions into bands.  However, there are cases when users will want to know that no bands
     * span across a given genomic position (e.g. when scatter-gathering jobs across a compute farm).  The option below enables users to break bands at
     * pre-defined positions.  For example, a value of 10,000 would mean that we would ensure that no bands span across chr1:10000, chr1:20000, etc.
     *
     * Note that the --convertToBasePairResolution argument is just a special case of this argument with a value of 1.
     */
    @Argument(fullName="breakBandsAtMultiplesOf", shortName="breakBandsAtMultiplesOf", doc = "If > 0, reference bands will be broken up at genomic positions that are multiples of this number", required=false)
    protected int multipleAtWhichToBreakBands = 0;

    private GenomeLocParser genomeLocParser;

    @Override
    public void onTraversalStart() {
        // take care of the VCF headers
        final Map<String, VCFHeader> vcfHeaders = GATKVCFUtils.getVCFHeadersFromRods(getToolkit());
        final Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), true);
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags

        final Set<String> samples = VcfUtils.getSortedSampleSet(vcfHeaders, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        final VCFHeader vcfHeader = new VCFHeader(headerLines, samples);
        vcfWriter.writeHeader(vcfHeader);

        // collect the actual rod bindings into a list for use later
        for ( final FeatureDataSource<VariantContext> variantCollection : variantCollections )
            variants.addAll(variantCollection.getRodBindings());

        genomeLocParser = getToolkit().getGenomeLocParser();

        // create the annotation engine
        annotationEngine = new VariantAnnotatorEngine(Arrays.asList(annotationGroupsToUse), annotationsToUse, Collections.<String>emptyList(), this, getToolkit());

        //now that we have all the VCF headers, initialize the annotations (this is particularly important to turn off RankSumTest dithering in integration tests)
        annotationEngine.invokeAnnotationInitializationMethods(headerLines);

        // optimization to prevent mods when we always just want to break bands
        if ( multipleAtWhichToBreakBands == 1 )
            USE_BP_RESOLUTION = true;
    }

    public PositionalState map(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return null;

        final GenomeLoc loc = ref.getLocus();
        return new PositionalState(tracker.getValues(variants, loc), ref.getBases(), loc);
    }

    public OverallState reduceInit() {
        return new OverallState();
    }

    public OverallState reduce(final PositionalState startingStates, final OverallState previousState) {
        if ( startingStates == null )
            return previousState;

        if ( !startingStates.VCs.isEmpty() ) {
            if ( ! okayToSkipThisSite(startingStates, previousState) )
                endPreviousStates(previousState, startingStates.loc.incPos(-1), startingStates, false);
            previousState.VCs.addAll(startingStates.VCs);
            for(final VariantContext vc : previousState.VCs){
                previousState.samples.addAll(vc.getSampleNames());
            }

        }

        if ( breakBand(startingStates.loc) || containsEndingContext(previousState.VCs, startingStates.loc.getStart()) ) {
            endPreviousStates(previousState, startingStates.loc, startingStates, true);
        }

        return previousState;
    }


    private void setupVCFWriter(VCFHeader inputVCFHeader, SampleList samples) {
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());
        headerLines.addAll(getDefaultToolVCFHeaderLines());

        // Remove GCVFBlocks
        headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().startsWith(GVCF_BLOCK));
asd
        headerLines.addAll(annotationEngine.getVCFAnnotationDescriptions());
        headerLines.addAll(genotypingEngine.getAppropriateVCFInfoHeaders());

        // add headers for annotations added by this tool
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags
        if ( dbsnp.dbsnp != null  ) {
            VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DBSNP_KEY);
        }

        vcfWriter = createVCFWriter(outputFile);

        final Set<String> sampleNameSet = samples.asSetOfSamples();
        final VCFHeader vcfHeader = new VCFHeader(headerLines, new TreeSet<>(sampleNameSet));
        vcfWriter.writeHeader(vcfHeader);
    }

    /**
     * Should we break bands at the given position?
     *
     * @param loc  the genomic location to evaluate against
     *
     * @return true if we should ensure that bands should be broken at the given position, false otherwise
     */
    private boolean breakBand(final GenomeLoc loc) {
        return USE_BP_RESOLUTION ||
                (loc != null && multipleAtWhichToBreakBands > 0 && (loc.getStart()+1) % multipleAtWhichToBreakBands == 0);  // add +1 to the loc because we want to break BEFORE this base
    }

    /**
     * Is it okay to skip the given position?
     *
     * @param startingStates  state information for this position
     * @param previousState   state information for the last position for which we created a VariantContext
     * @return true if it is okay to skip this position, false otherwise
     */
    private boolean okayToSkipThisSite(final PositionalState startingStates, final OverallState previousState) {
        final int thisPos = startingStates.loc.getStart();
        final GenomeLoc lastPosRun = previousState.prevPos;
        Set<String> intersection = new HashSet<String>(startingStates.samples);
        intersection.retainAll(previousState.samples);

        //if there's a starting VC with a sample that's already in a current VC, don't skip this position
        return lastPosRun != null && thisPos == lastPosRun.getStart() + 1 && intersection.isEmpty();
    }

    /**
     * Does the given list of VariantContexts contain any whose context ends at the given position?
     *
     * @param VCs  list of VariantContexts
     * @param pos  the position to check against
     * @return true if there are one or more VCs that end at pos, false otherwise
     */
    private boolean containsEndingContext(final List<VariantContext> VCs, final int pos) {
        if ( VCs == null ) throw new IllegalArgumentException("The list of VariantContexts cannot be null");

        for ( final VariantContext vc : VCs ) {
            if ( isEndingContext(vc, pos) )
                return true;
        }
        return false;
    }

    /**
     * Does the given variant context end (in terms of reference blocks, not necessarily formally) at the given position.
     * Note that for the purposes of this method/tool, deletions are considered to be single base events (as opposed to
     * reference blocks), hence the check for the number of alleles (because we know there will always be a <NON_REF> allele).
     *
     * @param vc   the variant context
     * @param pos  the position to query against
     * @return true if this variant context "ends" at this position, false otherwise
     */
    private boolean isEndingContext(final VariantContext vc, final int pos) {
        return vc.getNAlleles() > 2 || vc.getEnd() == pos;
    }

    /**
     * Disrupt the VariantContexts so that they all stop at the given pos, write them out, and put the remainder back in the list.
     * @param state   the previous state with list of active VariantContexts
     * @param pos   the position for the starting VCs
     * @param startingStates the state for the starting VCs
     * @param atCurrentPosition  indicates whether we output a variant at the current position, independent of VCF start/end, i.e. in BP resolution mode
     */
    private void endPreviousStates(final OverallState state, final GenomeLoc pos, final PositionalState startingStates, boolean atCurrentPosition) {

        final byte refBase = startingStates.refBases[0];
        //if we're in BP resolution mode or a VC ends at the current position then the reference for the next output VC (refNextBase)
        // will be advanced one base
        final byte refNextBase = (atCurrentPosition) ? (startingStates.refBases.length > 1 ? startingStates.refBases[1] : (byte)'N' ): refBase;

        final List<VariantContext> stoppedVCs = new ArrayList<>(state.VCs.size());

        for ( int i = state.VCs.size() - 1; i >= 0; i-- ) {
            final VariantContext vc = state.VCs.get(i);
            //the VC for the previous state will be stopped if its position is previous to the current position or it we've moved to a new contig
            if ( vc.getStart() <= pos.getStart() || !vc.getChr().equals(pos.getContig())) {

                stoppedVCs.add(vc);

                // if it was ending anyways, then remove it from the future state
                if ( vc.getEnd() == pos.getStart()) {
                    state.samples.removeAll(vc.getSampleNames());
                    state.VCs.remove(i);
                    continue; //don't try to remove twice
                }

                //if ending vc is the same sample as a starting VC, then remove it from the future state
                if(startingStates.VCs.size() > 0 && !atCurrentPosition && startingStates.samples.containsAll(vc.getSampleNames())) {
                    state.samples.removeAll(vc.getSampleNames());
                    state.VCs.remove(i);
                }
            }
        }

        //output the stopped VCs if there is no previous output (state.prevPos == null) or our current position is past
        // the last write position (state.prevPos)
        //NOTE: BP resolution with have current position == state.prevPos because it gets output via a different control flow
        if ( !stoppedVCs.isEmpty() &&  (state.prevPos == null || pos.isPast(state.prevPos) )) {
            final GenomeLoc gLoc = genomeLocParser.createGenomeLoc(stoppedVCs.get(0).getChr(), pos.getStart());

            // we need the specialized merge if the site contains anything other than ref blocks
            final VariantContext mergedVC;
            if ( containsTrueAltAllele(stoppedVCs) )
                mergedVC = ReferenceConfidenceVariantContextMerger.merge(stoppedVCs, gLoc, refBase, false, false, annotationEngine);
            else
                mergedVC = referenceBlockMerge(stoppedVCs, state, pos.getStart());

            vcfWriter.add(mergedVC);
            state.prevPos = gLoc;
            state.refAfterPrevPos = refNextBase;
        }
    }

    /**
     * Combine a list of reference block VariantContexts.
     * We can't use GATKVariantContextUtils.simpleMerge() because it is just too slow for this sort of thing.
     *
     * @param VCs   the variant contexts to merge
     * @param state the state object
     * @param end   the end of this block (inclusive)
     * @return a new merged VariantContext
     */
    private VariantContext referenceBlockMerge(final List<VariantContext> VCs, final OverallState state, final int end) {

        final VariantContext first = VCs.get(0);

        // ref allele and start
        final Allele refAllele;
        final int start;
        if ( state.prevPos == null || !state.prevPos.getContig().equals(first.getChr()) || first.getStart() >= state.prevPos.getStart() + 1) {
            start = first.getStart();
            refAllele = first.getReference();
        } else {
            start = state.prevPos.getStart() + 1;
            refAllele = Allele.create(state.refAfterPrevPos, true);
        }

        // attributes
        final Map<String, Object> attrs = new HashMap<>(1);
        if ( !USE_BP_RESOLUTION && end != start )
            attrs.put(VCFConstants.END_KEY, Integer.toString(end));

        // genotypes
        final GenotypesContext genotypes = GenotypesContext.create();
        for ( final VariantContext vc : VCs ) {
            for ( final Genotype g : vc.getGenotypes() )
                genotypes.add(new GenotypeBuilder(g).alleles(GATKVariantContextUtils.noCallAlleles(g.getPloidy())).make());
        }

        return new VariantContextBuilder("", first.getChr(), start, end, Arrays.asList(refAllele, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE)).attributes(attrs).genotypes(genotypes).make();
    }

    /**
     * Does the given list of VariantContexts contain any with an alternate allele other than <NON_REF>?
     *
     * @param VCs  list of VariantContexts
     * @return true if there are one or more VCs that contain a true alternate allele, false otherwise
     */
    private boolean containsTrueAltAllele(final List<VariantContext> VCs) {
        if ( VCs == null ) throw new IllegalArgumentException("The list of VariantContexts cannot be null");

        for ( final VariantContext vc : VCs ) {
            if ( vc.getNAlleles() > 2 )
                return true;
        }
        return false;
    }

    @Override
    public void onTraversalDone(final OverallState state) {
        // there shouldn't be any state left unless the user cut in the middle of a gVCF block
        if ( !state.VCs.isEmpty() )
            logger.warn("You have asked for an interval that cuts in the middle of one or more gVCF blocks. Please note that this will cause you to lose records that don't end within your interval.");
    }
}
