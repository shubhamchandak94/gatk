package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by tsato on 8/7/17.
 */
public class PerContextData {
    private final String referenceContext;
    private final List<Allele> alleles;
    private final List<String> positions; // for debugging

    // the number of data points (i.e. loci) with this 3-mer in the reference context
    private int numLoci = 0;

    private List<Integer> depths;
    private List<Short> altDepths;
    private List<Short> altF1R2Depths;

    List<Double> responsibilities;

    // k-dimensional vectors of hyperparameters, where k is the number of available states of the latent variable z
    double[][] mixtureWeights; // pi. rows are the alleles, columns the mixture components. each row must add up to 1.0

    double[] alleleFractions; // f

    double[] altF1R2Fractions; // theta

    public PerContextData(final String referenceContext){
        this.referenceContext = referenceContext;
        alleles = new ArrayList<>(ContextDependentArtifactFilter.DEFAULT_INITIAL_LIST_SIZE);
        positions = new ArrayList<>(ContextDependentArtifactFilter.DEFAULT_INITIAL_LIST_SIZE);

        depths = new ArrayList<>(ContextDependentArtifactFilter.DEFAULT_INITIAL_LIST_SIZE);
        altDepths = new ArrayList<>(ContextDependentArtifactFilter.DEFAULT_INITIAL_LIST_SIZE);
        altF1R2Depths = new ArrayList<>(ContextDependentArtifactFilter.DEFAULT_INITIAL_LIST_SIZE);

        mixtureWeights = new double[ContextDependentArtifactFilterEngine.NUM_ALLELES][ContextDependentArtifactFilterEngine.NUM_STATUSES];
        alleleFractions = new double[ContextDependentArtifactFilterEngine.NUM_STATUSES];
        altF1R2Fractions = new double[ContextDependentArtifactFilterEngine.NUM_STATUSES];

        responsibilities = new ArrayList<>(ContextDependentArtifactFilter.DEFAULT_INITIAL_LIST_SIZE);
    }

    public void addNewExample(final int depth, final short altDepth, final short altF1R2Depth, final Allele allele, final AlignmentContext alignmentContext){
        depths.add(depth);
        altDepths.add(altDepth);
        altF1R2Depths.add(altF1R2Depth);
        alleles.add(allele);
        positions.add(String.format("%s:%s-%s", alignmentContext.getContig(), alignmentContext.getStart(), alignmentContext.getEnd()));

        numLoci++;
    }

    // debug method
    // TODO: choose a better name for this method (Philip Guo called it something else...)
    private void validateInternalStructures(){
        final double EPSILON = 1e-3;
        for (int i = 0; i < ContextDependentArtifactFilterEngine.NUM_ALLELES; i++){
            Utils.validate(Math.abs(MathUtils.sum(mixtureWeights[i]) - 1.0) < EPSILON, "mixture weights must add up to 1.0");
        }
    }

    // getters
    public String getReferenceContext(){ return referenceContext; }

    public List<Allele> getAlleles(){ return alleles; }

    public int getNumLoci(){ return numLoci; }

    public List<Integer> getDepths(){ return depths; }

    public List<Short> getAltDepths(){ return altDepths; }

    public List<Short> getAltF1R2Depths(){ return altF1R2Depths; }

}
