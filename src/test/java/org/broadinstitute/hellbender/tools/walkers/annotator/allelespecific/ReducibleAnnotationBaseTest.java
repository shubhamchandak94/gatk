package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.hadoop.yarn.webapp.hamlet.Hamlet;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

/**
 * A class for storing end-end integration tests over artificial data intended to test allele specific annotation implementations,
 * As of 8/29/17, test output was matched against GATK3 implementations of combineGVCFs().
 */
public abstract class ReducibleAnnotationBaseTest extends BaseTest {
    private FeatureDataSource<VariantContext> VCFA = new FeatureDataSource<>(getTestFile("NA12878.AS.chr20snippet.g.vcf"));
    private FeatureDataSource<VariantContext> VCFB = new FeatureDataSource<>(getTestFile("NA12892.AS.chr20snippet.g.vcf"));
    private FeatureDataSource<VariantContext> CombineVCFOutput = new FeatureDataSource<>(getTestFile("CombineGVCFs.output.vcf"));

    @Override
    public String getToolTestDataDir(){
        return "src/test/resources/" +this.getClass().getPackage().getName().replace(".","/") + "/";
    }

    @DataProvider (name = "interestingSitesCombineResults")
    public Object[][] makeAnnotateRsIDData() {
        List<Object[]> tests = new ArrayList<>();

        // these are hand picked sites from the allele specific unit tests for combinegvcfs that triggered a combine in GATK3
        Integer[] interestingLocs = {10087820,10433312, 10433322, 10433324, 10433326, 10433328, 10433382, 10433391, 10433468, 10433560, 10433575, 10433594, 10433955,
                10434177, 10434384, 10435067, 10434258, 10436227, 10684106};
        List<SimpleInterval> intervals = Arrays.stream(interestingLocs).map(m -> new SimpleInterval("20", m, m)).collect(Collectors.toList());
        for(SimpleInterval loc: intervals) {
            VariantContext a = VCFA.query(loc).next();
            VariantContext b = VCFB.query(loc).next();
            VariantContext result = CombineVCFOutput.query(loc).next();
            tests.add(new Object[]{Arrays.asList(a, b),result});
        }
        return tests.toArray(new Object[][]{});
    }

    // Method that determines whether two variant contexts have equivalent allele specific annotations regardless of allele ordering
    public static Boolean alleleSpecificAnnotationEquals(VariantContext a, VariantContext b, String annotation) {
        List<Allele> Aalleles = a.getAlleles();
        String[] Aannot = String.join(",",a.getAttributeAsStringList(annotation, "")).split("\\|",-1);
        String[] Bannot = String.join(",",b.getAttributeAsStringList(annotation, "")).split("\\|",-1);
        if (Arrays.equals(Aannot, Bannot)) {
            return true;
        } if (Aannot.length!=Bannot.length) {
            return false;
        }
        for (int i = 0; i < Aalleles.size(); i++) {
            Allele al = Aalleles.get(i);

            int k = b.getAlleleIndex(al);

            if (!Aannot[i].equals(Bannot[k])) {
                return false;
            }
        }
        return true;
     }

    /*
     * Methods that must be overriden in order for the automatic GATK3 combineGVCFs tests to be run.
     */
    protected abstract List<String> getAnnotationsToUse();
    protected abstract String getRawKey();



    // NOTE: this code is mimicking the behavior of GATK3 combineGVCFS insofar as it is important for the annotations
    @Test (dataProvider = "interestingSitesCombineResults")
    public void testAnnotationCombineGVCF(List<VariantContext> VCs, VariantContext result) throws Exception {
        ReferenceConfidenceVariantContextMerger merger = new MiniMerger(getAnnotationsToUse());
        VariantContext merged = merger.merge(VCs, new SimpleInterval(result.getContig(), result.getStart(), result.getStart()), result.getReference().getBases()[0],false, false );
        Assert.assertTrue(alleleSpecificAnnotationEquals(merged,result,getRawKey()));
    }


    /**
     * Subclass of ReferenceConfidenceVariantContextMerger which contains the logic necessary for handling of
     * allele specific annotaions to test the combine logic for client tools.
     */
    private class MiniMerger extends ReferenceConfidenceVariantContextMerger {
        VariantAnnotatorEngine annotatorEngine;

        public MiniMerger(List<String> annotationsToUse) {
            final List<String> annotationGroupsToUse = Collections.emptyList();
            final List<String> annotationsToExclude = Collections.emptyList();
            final FeatureInput<VariantContext> dbSNPBinding = null;
            final List<FeatureInput<VariantContext>> features = Collections.emptyList();
            annotatorEngine = VariantAnnotatorEngine.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding, features);
        }

        /**
         * Merges VariantContexts from gVCFs into a single hybrid.
         * Assumes that none of the input records are filtered.
         *
         * @param vcs     collection of unsorted genomic vcs
         * @param loc     the current location
         * @param refBase the reference allele to use if all contexts in the VC are spanning (i.e. don't start at the location in loc); if null, we'll return null in this case
         * @param removeNonRefSymbolicAllele if true, remove the <NON_REF> allele from the merged VC
         * @param samplesAreUniquified  if true, sample names have been uniquified
         * @return new VariantContext representing the merge of all vcs or null if it not relevant
         */
        public VariantContext merge(final List<VariantContext> vcs, final Locatable loc, final Byte refBase,
                                    final boolean removeNonRefSymbolicAllele, final boolean samplesAreUniquified) {
            Utils.nonEmpty(vcs);

            // establish the baseline info (sometimes from the first VC)
            final String name = vcs.get(0).getSource();

            // ref allele
            final Allele refAllele = determineReferenceAlleleGivenReferenceBase(vcs, loc, refBase);
            if ( refAllele == null ) {
                return null;
            }

            // In this list we hold the mapping of each variant context alleles.
            final List<VCWithNewAlleles> vcAndNewAllelePairs = new ArrayList<>(vcs.size());

            // cycle through and add info from the other vcs
            for ( final VariantContext vc : vcs ) {
                // if this context doesn't start at the current location then it must be a spanning event (deletion or ref block)
                final boolean isSpanningEvent = loc.getStart() != vc.getStart();
                vcAndNewAllelePairs.add(new VCWithNewAlleles(vc, isSpanningEvent ? replaceWithNoCallsAndDels(vc) : remapAlleles(vc, refAllele), isSpanningEvent));
            }

            final List<Allele> allelesList = collectTargetAlleles(vcAndNewAllelePairs, refAllele, removeNonRefSymbolicAllele);

            final Set<String> rsIDs = new LinkedHashSet<>(1); // most of the time there's one id
            int depth = 0;
            final Map<String, List<ReducibleAnnotationData<Object>>> annotationMap = new LinkedHashMap<>();

            final GenotypesContext genotypes = GenotypesContext.create();

            for ( final VCWithNewAlleles vcWithNewAlleles : vcAndNewAllelePairs ) {
                final VariantContext vc = vcWithNewAlleles.getVc();
                final List<Allele> remappedAlleles = vcWithNewAlleles.getNewAlleles();

                genotypes.addAll(mergeRefConfidenceGenotypes(vc, remappedAlleles, allelesList, samplesAreUniquified));
                depth += calculateVCDepth(vc);

                if ( loc.getStart() != vc.getStart() ) {
                    continue;
                }

                // special case ID (just preserve it)
                if ( vc.hasID() ) {
                    rsIDs.add(vc.getID());
                }

                // add attributes
                addReferenceConfidenceAttributes(vcWithNewAlleles, annotationMap);
            }

            final Map<String, Object> attributes = mergeAttributes(depth, allelesList, annotationMap);

            final String ID = rsIDs.isEmpty() ? VCFConstants.EMPTY_ID_FIELD : String.join(",", rsIDs);

            // note that in order to calculate the end position, we need a list of alleles that doesn't include anything symbolic
            final VariantContextBuilder builder = new VariantContextBuilder()
                    .source(name)
                    .id(ID)
                    .alleles(allelesList)
                    .chr(loc.getContig())
                    .start(loc.getStart())
                    .computeEndFromAlleles(nonSymbolicAlleles(allelesList), loc.getStart(), loc.getStart())
                    .genotypes(genotypes).unfiltered()
                    .attributes(new TreeMap<>(attributes)).log10PError(CommonInfo.NO_LOG10_PERROR);  // we will need to re-genotype later

            return builder.make();
        }

        protected <T extends Comparable<? super T>> void addReferenceConfidenceAttributes(final VCWithNewAlleles vcPair,
                                                                                                 final Map<String, List<ReducibleAnnotationData<Object>>> annotationMap) {
            for (final Map.Entry<String, Object> p : vcPair.getVc().getAttributes().entrySet()) {
                final String key = p.getKey();
                final List<Object> valueList = vcPair.getVc().getAttributeAsList(key);

                // add the annotation values to a list for combining later
                List<ReducibleAnnotationData<Object>> values = annotationMap.get(key);
                if (values == null) {
                    values = new ArrayList<>();
                    annotationMap.put(key, values);
                }
                String combinedString = "";
                for(int i=0; i < valueList.size(); i++) {
                    if (i > 0)
                        combinedString += ",";
                    combinedString += valueList.get(i);
                }

                ReducibleAnnotationData<Object> pairData = new AlleleSpecificAnnotationData<>(vcPair.getNewAlleles(), combinedString);
                values.add(pairData);
                annotationMap.put(key, values);
            }

        }


        public Map<String, Object> mergeAttributes(int depth, List<Allele> alleleList, Map<String, List<ReducibleAnnotationData<Object>>> annotationMap) {
            final Map<String, Object> attributes = new LinkedHashMap<>();

            attributes.putAll(annotatorEngine.combineAnnotations(alleleList, annotationMap));
            annotationMap.entrySet().stream()
                    .filter(p -> !p.getValue().isEmpty())
                    .forEachOrdered(p -> attributes.put(p.getKey(), (p.getValue()))
                    );

            if ( depth > 0 ) {
                attributes.put(VCFConstants.DEPTH_KEY, String.valueOf(depth));
            }

            // remove stale AC and AF based attributes
            removeStaleAttributesAfterMerge(attributes);
            return attributes;
        }
    }

}