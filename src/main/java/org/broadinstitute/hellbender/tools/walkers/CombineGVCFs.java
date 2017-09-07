package org.broadinstitute.hellbender.tools.walkers;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.MultiVariantWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
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

    private static final String GVCF_BLOCK = "GVCFBlock";

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
    @Advanced
    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", optional=true)
    private List<String> annotationGroupsToUse = new ArrayList<>(Arrays.asList(new String[]{StandardAnnotation.class.getSimpleName()}));

    @Advanced
    @Argument(fullName="annotationsToExclude", shortName="AX", doc="One or more specific annotations to exclude from recomputation", optional=true)
    private List<String> annotationsToExclude = new ArrayList<>();

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The combined GVCF output file", optional=false)
    private File outputFile;

    @Argument(fullName="convertToBasePairResolution", shortName="bpResolution", doc = "If specified, convert banded gVCFs to all-sites gVCFs", optional=true)
    protected boolean USE_BP_RESOLUTION = false;

    /**
     * To reduce file sizes our gVCFs group similar reference positions into bands.  However, there are cases when users will want to know that no bands
     * span across a given genomic position (e.g. when scatter-gathering jobs across a compute farm).  The option below enables users to break bands at
     * pre-defined positions.  For example, a value of 10,000 would mean that we would ensure that no bands span across chr1:10000, chr1:20000, etc.
     *
     * Note that the --convertToBasePairResolution argument is just a special case of this argument with a value of 1.
     */
    @Argument(fullName="breakBandsAtMultiplesOf", shortName="breakBandsAtMultiplesOf", doc = "If > 0, reference bands will be broken up at genomic positions that are multiples of this number", optional=true)
    protected int multipleAtWhichToBreakBands = 0;

    public static final String IGNORE_VARIANTS_THAT_START_OUTSIDE_INTERVAL = "ignore_variants_starting_outside_interval";
    /**
     * This option can only be activated if intervals are specified.
     *
     * This exists to mimic GATK3 interval traversal patterns
     */
    @Advanced
    @Argument(fullName= IGNORE_VARIANTS_THAT_START_OUTSIDE_INTERVAL,
            doc="Restrict variant output to sites that start within provided intervals",
            optional=true)
    private boolean ignoreIntervalsOutsideStart = false;

    /**
     * The rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate. Note that dbSNP is not used in any way for the calculations themselves.
     */

//    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();
//    public FeatureInput<VariantContext> getDbsnpRodBinding() { return dbsnp.dbsnp; }
//    public List<RodBinding<VariantContext>> getCompRodBindings() { return Collections.emptyList(); }
//    public RodBinding<VariantContext> getSnpEffRodBinding() { return null; }
//    public List<RodBinding<VariantContext>> getResourceRodBindings() { return Collections.emptyList(); }
//    public boolean alwaysAppendDbsnpId() { return false; }

    // the annotation engine
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();
    private VariantAnnotatorEngine annotationEngine;
    List<VariantContext> currentVariants = new ArrayList<>();
    PositionalState currentPositionalState;
    VariantContextWriter vcfWriter;
    ReferenceConfidenceVariantContextMerger referenceConfidenceVariantContextMerger = new ReferenceConfidenceVariantContextMerger();
    SAMSequenceDictionary sequenceDictionary;


    // STATE WHICH IS ACCUMULATED BETWEEN RUNS
    final LinkedList<VariantContext> VCs = new LinkedList<>();
    final Set<String> samples = new HashSet<>();
    SimpleInterval prevPos = null;
    byte refAfterPrevPos;
    byte[] storedReference;
    private OverlapDetector overlapDetector;
    private boolean hasReduced = false;
    private SimpleInterval storedReferenceLoc;


    /**
     * This method keeps track of all the variants it is passed and will feed all the variants that start at the same
     * site to the reduce method.
     *
     * @param variant Current variant being processed.
     * @param readsContext Reads overlapping the current variant. Will be an empty, but non-null, context object
     *                     if there is no backing source of reads data (in which case all queries on it will return
     *                     an empty array/iterator)
     * @param referenceContext Reference bases spanning the current variant. Will be an empty, but non-null, context object
     *                         if there is no backing source of reference data (in which case all queries on it will return
     *                         an empty array/iterator). Can request extra bases of context around the current variant's interval
     *                         by invoking {@link ReferenceContext#setWindow}
     *                         on this object before calling {@link ReferenceContext#getBases}
     * @param featureContext Features spanning the current variant. Will be an empty, but non-null, context object
     *                       if there is no backing source of Feature data (in which case all queries on it will return an
     */
    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        // Filtering out reads which start outside of the specified intervals
        if (ignoreIntervalsOutsideStart && !overlapDetector.overlapsAny(new SimpleInterval(variant.getContig(),variant.getStart(),variant.getStart()))) {
            return;
        }

        // Collecting all the reads that start at a particular base into one.
        if (currentVariants.isEmpty()) {
            currentVariants.add(variant);
        } else if (!currentVariants.get(0).getContig().equals(variant.getContig())
                || currentVariants.get(0).getStart()<variant.getStart()) {
            // Emptying any sites which should emit a new VC since the last one
            reduceQueuedState();
            currentVariants.clear();
            currentVariants.add(variant);
            // TODO BE VERY CLEAR ABOUT THIS
        } else {
            currentVariants.add(variant);
        }
        referenceContext.setWindow(1,1);
        updatePositionalState(currentVariants, referenceContext);
    }

    private void reduceQueuedState() {
        if (hasReduced == true) {
            // TODO eliminate the amount of string comparison
            createIntermediateVariants(prevPos == null ? new SimpleInterval(VCs.get(0).getContig(), VCs.get(0).getStart(),
                    (VCs.get(0).getContig().equals(currentPositionalState.loc.getContig())
                            ? currentPositionalState.loc.getStart() - 1
                            : VCs.get(0).getStart() + storedReference.length - 1 ))

                    : new SimpleInterval(prevPos.getContig(), prevPos.getEnd(),
                    prevPos.getContig().equals(currentPositionalState.loc.getContig())
                            ? currentPositionalState.loc.getStart() - 1
                            : prevPos.getStart() + storedReference.length - 1));// TODO this value is bad,
        }

        //TODO ^^^ YUCK!
//==============================================//TODO, these values are bad becaus of problems
        reduce(currentPositionalState);

        // Update the stored reference if it has a later stop position than the current stored reference
        if ( (storedReferenceLoc == null) ||
                (!currentPositionalState.loc.getContig().equals(storedReferenceLoc.getContig()) ) ||
                (storedReferenceLoc.getStart() + storedReference.length < currentPositionalState.loc.getStart() + currentPositionalState.refBases.length)) {
            storedReference = currentPositionalState.refBases;
            storedReferenceLoc = currentPositionalState.refBasesStart;
        }
    }


    /**
     * Method which calculates at which sites between the last created variant context and the added ones that a
     * break should be created and calls the appropriate method.
     *
     */
    @VisibleForTesting
    void createIntermediateVariants(SimpleInterval intervalToClose) {
        Set<Integer> sitesToStop = new HashSet<>();
        // Perform any band breaking that needs to be done since the last one
        if ( multipleAtWhichToBreakBands > 0) {
            // TODO figure out +1 from previous code
            for (int i = intervalToClose.getStart()/multipleAtWhichToBreakBands; i < (intervalToClose.getEnd()-1)/multipleAtWhichToBreakBands; i++) {
                sitesToStop.add(multipleAtWhichToBreakBands*(i+1));
            }
        }

        for (VariantContext vc : VCs) {
            if (vc.getNAlleles() > 2) {
                //TODO BE VERY SURE ABOUT ENDING AFTER WHEN WE WANT!!!!!!!!!!!!!!!
                for (int i = vc.getStart(); i <= vc.getEnd(); i++ ) {
                    sitesToStop.add(i);
                }
            } else if (vc.getEnd() <= intervalToClose.getEnd()) {
                sitesToStop.add(vc.getEnd());
            }
        }

        List<Integer> stoppedPlaces = new ArrayList<>(sitesToStop);
        stoppedPlaces.sort(Comparator.naturalOrder());//TODO comparitor magic

        for (int i : stoppedPlaces) {
            //TODO be sure about these reference bases and how to get them
            SimpleInterval loc = new SimpleInterval(intervalToClose.getContig(),i,i);
            if (( i <= intervalToClose.getEnd() && i>= intervalToClose.getStart()) && (overlapDetector==null || overlapDetector.overlapsAny(loc))) {
                byte[] refBases = Arrays.copyOfRange(storedReference, i - storedReferenceLoc.getStart(), i - storedReferenceLoc.getStart() + 1);
                PositionalState tmp = new PositionalState(Collections.emptyList(), refBases, new SimpleInterval(intervalToClose.getContig(), i, i));
                endPreviousStates(tmp.loc, refBases, tmp, true);
            }
        }

    }

    /**
     * Method which ensures that the currentPositionalState object to hold the longest set of reference bases
     *
     * @param currentVariants
     */
    private void updatePositionalState(List<VariantContext> currentVariants, ReferenceContext referenceContext) {
        if (currentVariants.size()==1 ) {
            currentPositionalState = new PositionalState(new ArrayList<>(currentVariants), referenceContext.getBases(), referenceContext.getInterval());
            currentPositionalState.refBasesStart = referenceContext.getWindow();
        } else {
            currentPositionalState.VCs.clear();
            currentPositionalState.VCs.addAll(currentVariants);
            currentPositionalState.refBases = (referenceContext.getBases().length > currentPositionalState.refBases.length?
                    referenceContext.getBases() : currentPositionalState.refBases);
            currentPositionalState.refBasesStart = referenceContext.getWindow();
        }
    }

    protected final class PositionalState {
        final List<VariantContext> VCs;
        final Set<String> samples = new HashSet<>();
        byte[] refBases;
        SimpleInterval loc;
        public SimpleInterval refBasesStart;

        public PositionalState(final List<VariantContext> VCs, final byte[] refBases, final SimpleInterval loc) {
            this.VCs = VCs;
            for(final VariantContext vc : VCs){
                samples.addAll(vc.getSampleNames());
            }
            this.refBases = refBases;
            this.loc = loc;
        }
    }


    @Override
    public void onTraversalStart() {
        final SortedSet<String> samples = getSamplesForVariants();

        final VCFHeader vcfHeader = new VCFHeader(getHeaderForVariants().getMetaDataInInputOrder(), samples);//TODO figure if the samples actually get imported

        // create the annotation engine
        annotationEngine = VariantAnnotatorEngine.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbsnp.dbsnp, Collections.EMPTY_LIST);

        setupVCFWriter(vcfHeader, new IndexedSampleList(samples));

        referenceConfidenceVariantContextMerger = new ReferenceConfidenceVariantContextMerger(annotationEngine);

        //now that we have all the VCF headers, initialize the annotations (this is particularly important to turn off RankSumTest dithering in integration tests)'
        sequenceDictionary = getBestAvailableSequenceDictionary();

        overlapDetector = hasIntervals() ? OverlapDetector.create(intervalArgumentCollection.getIntervals(sequenceDictionary)) :
                null;

        // optimization to prevent mods when we always just want to break bands
        if ( multipleAtWhichToBreakBands == 1 )
            USE_BP_RESOLUTION = true;
    }

    /**
     * Method which calls endPreviousStates at the appropriate places on the given a new startingStates object
     * and an OverallState object corresponding to the currently accumulated reads.
     *
     * @param startingStates
     * @return
     */
    public void reduce(final PositionalState startingStates) {
        hasReduced = true;
        if ( !startingStates.VCs.isEmpty() ) {
            if ( ! okayToSkipThisSite(startingStates) ) {
                SimpleInterval loc = startingStates.loc;
                endPreviousStates( new SimpleInterval(loc.getContig(),loc.getStart()-1,loc.getStart()-1), Arrays.copyOfRange(startingStates.refBases, 1,startingStates.refBases.length), startingStates, false);
            }
            VCs.addAll(startingStates.VCs);
            for(final VariantContext vc : VCs){
                samples.addAll(vc.getSampleNames());
            }
        }
    }


    private void setupVCFWriter(VCFHeader inputVCFHeader, SampleList samples) {
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());
        headerLines.addAll(getDefaultToolVCFHeaderLines());

        // Remove GCVFBlocks
        headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().startsWith(GVCF_BLOCK));

        headerLines.addAll(annotationEngine.getVCFAnnotationDescriptions());

        // add headers for annotations added by this tool
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
    private boolean breakBand(final SimpleInterval loc) {
        return USE_BP_RESOLUTION ||
                (loc != null && multipleAtWhichToBreakBands > 0 && (loc.getStart()+1) % multipleAtWhichToBreakBands == 0);  // add +1 to the loc because we want to break BEFORE this base
    }

    /**
     * Is it okay to skip the given position?
     *
     * @param startingStates  state information for this position
     * @return true if it is okay to skip this position, false otherwise
     */
    private boolean okayToSkipThisSite(final PositionalState startingStates) {
        final int thisPos = startingStates.loc.getStart();
        final SimpleInterval lastPosRun = prevPos;
        Set<String> intersection = new HashSet<String>(startingStates.samples);
        intersection.retainAll(samples);

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
     * @param pos   the position for the starting VCs
     * @param startingStates the state for the starting VCs
     * @param atCurrentPosition  indicates whether we output a variant at the current position, independent of VCF start/end, i.e. in BP resolution mode
     */
    private void endPreviousStates( final SimpleInterval pos, final byte[] refBases, final PositionalState startingStates, boolean atCurrentPosition) {

        final byte refBase = refBases[0];
        //if we're in BP resolution mode or a VC ends at the current position then the reference for the next output VC (refNextBase)
        // will be advanced one base
        final byte refNextBase = (atCurrentPosition) ? (refBases.length > 1 ? refBases[1] : (byte)'N' ): refBase;

        final List<VariantContext> stoppedVCs = new ArrayList<>(VCs.size());

        for ( int i = VCs.size() - 1; i >= 0; i-- ) {
            final VariantContext vc = VCs.get(i);
            //the VC for the previous state will be stopped if its position is previous to the current position or it we've moved to a new contig
            if ( vc.getStart() <= pos.getStart() || !vc.getContig().equals(pos.getContig())) {

                stoppedVCs.add(vc);

                // if it was ending anyways, then remove it from the future state
                if ( vc.getEnd() == pos.getStart()) {
                    samples.removeAll(vc.getSampleNames());
                    VCs.remove(i);
                    continue; //don't try to remove twice
                }

                //if ending vc is the same sample as a starting VC, then remove it from the future state
                if(startingStates.VCs.size() > 0 && !atCurrentPosition && startingStates.samples.containsAll(vc.getSampleNames())) {
                    samples.removeAll(vc.getSampleNames());
                    VCs.remove(i);
                }
            }
        }

        //output the stopped VCs if there is no previous output (state.prevPos == null) or our current position is past
        // the last write position (state.prevPos)
        //NOTE: BP resolution with have current position == state.prevPos because it gets output via a different control flow
        if ( !stoppedVCs.isEmpty() &&  (prevPos == null || IntervalUtils.isAfter(pos,prevPos,sequenceDictionary) )) {
            final SimpleInterval gLoc = new SimpleInterval(stoppedVCs.get(0).getContig(), pos.getStart(), pos.getStart());

            // we need the specialized merge if the site contains anything other than ref blocks
            final VariantContext mergedVC;
            if ( containsTrueAltAllele(stoppedVCs) )
                mergedVC = referenceConfidenceVariantContextMerger.merge(stoppedVCs, gLoc, refBase, false, false);
            else
                mergedVC = referenceBlockMerge(stoppedVCs, pos.getStart());

            vcfWriter.add(mergedVC);
            prevPos = gLoc;
            refAfterPrevPos = refNextBase;
        }
    }

    /**
     * Combine a list of reference block VariantContexts.
     * We can't use GATKVariantContextUtils.simpleMerge() because it is just too slow for this sort of thing.
     *
     * @param VCs   the variant contexts to merge
     * @param end   the end of this block (inclusive)
     * @return a new merged VariantContext
     */
    private VariantContext referenceBlockMerge(final List<VariantContext> VCs, final int end) {

        final VariantContext first = VCs.get(0);

        // ref allele and start
        final Allele refAllele;
        final int start;
        if ( prevPos == null || !prevPos.getContig().equals(first.getContig()) || first.getStart() >= prevPos.getStart() + 1) {
            start = first.getStart();
            refAllele = first.getReference();
        } else {
            start = prevPos.getStart() + 1;
            refAllele = Allele.create(refAfterPrevPos, true);
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

        return new VariantContextBuilder("", first.getContig(), start, end, Arrays.asList(refAllele, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE)).attributes(attrs).genotypes(genotypes).make();
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
    public Object onTraversalSuccess() {
        // Clearing the accumulator
        if(currentVariants.isEmpty()){
            logger.warn("You have asked for an interval does not contain any data in source GVCFs");

        } else{
            reduceQueuedState();

            // Create variants with whatever remains
            SimpleInterval interval = prevPos != null ? new SimpleInterval(prevPos.getContig(), prevPos.getStart(), prevPos.getStart() + storedReference.length + 1) :
                    currentPositionalState.loc;
            createIntermediateVariants(interval);
        }

        //TODO the reference bases are almost certainly wrong here
        // there shouldn't be any state left unless the user cut in the middle of a gVCF block
        if ( !VCs.isEmpty() )
            logger.warn("You have asked for an interval that cuts in the middle of one or more gVCF blocks. Please note that this will cause you to lose records that don't end within your interval.");
        vcfWriter.close();
        return null;
    }
}
