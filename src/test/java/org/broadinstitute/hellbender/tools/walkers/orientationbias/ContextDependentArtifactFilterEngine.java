package org.broadinstitute.hellbender.tools.walkers.orientationbias;

import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * Created by tsato on 7/26/17.
 */
public class ContextDependentArtifactFilterEngine {
    // z \in { F1R2, F2R1, Balanced Hom Ref, Balanced Het, Balanced Hom Var }. Thus |z| = K = 5.
    static final int NUM_STATUSES = 5;

    static final String[] ALL_ALLELES = new String[] { "A", "C", "G", "T" };

    // A, C, G, or T, since we only look at SNP sites
    static final int NUM_POSSIBLE_ALLELES = ALL_ALLELES.length;

    // When the increase in likelihood falls under this value, we call the algorithm converged
    static final double CONVERGENCE_THRESHOLD = 1e-4;

    String referenceContext;

    // Observed data
    final PerContextData data;

    // N by K matrix of per-sample posterior probabilities of latent variable z
    // with the current estimates of hyperparameters pi, f, and theta
    double[][] responsibilities;

    final int numExamples;

    final Allele[] alleles = new Allele[NUM_POSSIBLE_ALLELES];



    /*** Hyperparameters of the model ***/

    // pi is the A by K matrix of weights (prior?) for the categorical variable z
    // For each allele a, the K-dimensional vector pi_a must be a valid probability distribution over z
    // In other words, each row must add up to 1
    double[][] pi = new double[NUM_POSSIBLE_ALLELES][NUM_STATUSES];

    // K-dimensional vector of probabilities for the binomial m given z, which represents the number of alt reads at site n
    double[] f;

    // K-dimensional vector of probabilities for the binomial x given z, which represents the number of F1R2 alt reads at site n
    double[] theta;


    public ContextDependentArtifactFilterEngine(final PerContextData data){
        this.data = data;
        numExamples = data.getNumLoci();
        responsibilities = new double[numExamples][NUM_STATUSES];

        // initialize responsibilities
        final double initialProbability = 1.0 / NUM_STATUSES;
        for (int i = 0; i < data.getNumLoci(); i++ ) {
            Arrays.fill(responsibilities[i], initialProbability);
        }

        referenceContext = data.getReferenceContext();

        // populate alleles
        final String referenceBase = referenceContext.substring(1, 2);

        for (int i = 0; i < ALL_ALLELES.length ; i++){
            alleles[i] = Allele.create(ALL_ALLELES[i], ALL_ALLELES[i].equals(referenceBase));
        }
    }

    public Hyperparameters runEMAlgorithm(){
        boolean converged = false;
        double oldLikelihood = 0;
        double newLikelihood = 0;

        while (!converged){
            takeMstep();
            assert newLikelihood > oldLikelihood : "M step must increase the likelihood";

            takeEstep();
            converged = checkLikelihoodHasConverged(oldLikelihood, newLikelihood);
        }

        return new Hyperparameters(pi, f, theta);
    }

    private void takeEstep(){

    }

    // given the current posterior distributions over z, compute the maximum likelihood estimate for
    // the categorical weights (pi), allele fractions (f), and alt F1R2 fraction (theta)
    private void takeMstep(){
        // compute responsibility-based statistics based on the current responsibilities

        // K-dimensional vector of effective sample counts for each class of z, weighted by the responsibilities
        // N_k in the docs
        final double[] effectiveCounts = GATKProtectedMathUtils.sumArrayFunction(0, numExamples,
                i -> responsibilities[i]);
        assert effectiveCounts.length == NUM_STATUSES : "effectiveCount must be a k-dimensional vector";

        // A by K matrix of effective counts. Each column must add up to N_k
        final double[][] effectiveCountsGivenAllele = new double[NUM_POSSIBLE_ALLELES][NUM_STATUSES];

        // TODO: optimize
        for (int i = 0; i < NUM_STATUSES; i++){
            for (int n = 0; n < numExamples; n++) {
                if (! alleles[i].equals(data.getAlleles().get(n))){
                    continue;
                }

                effectiveCountsGivenAllele[i] = MathArrays.ebeAdd(effectiveCountsGivenAllele[i], responsibilities[i]);
            }
        }

        // TODO: somehow, code below doens't work. Investigate
//        IntStream.range(0, NUM_POSSIBLE_ALLELES).forEach(i ->
//            effectiveCountGivenAllele[i] = IntStream.range(0, numExamples)
//                    .filter(n -> data.getAlleles().get(i).equals(alleles[n]))
//                    .boxed()
//                    .reduce(new double[NUM_STATUSES], (acc, n_a) -> MathArrays.ebeAdd(acc, responsibilities[n_a]))
//        );

        assert MathArrays.equals(GATKProtectedMathUtils.sumArrayFunction(0, NUM_POSSIBLE_ALLELES, a -> effectiveCountsGivenAllele[a]),
                effectiveCounts) : "N_{ak} summed over alleles must add up to N_k";

        // A-dimensional vector of effective count of samples for each allele, where A = 4 = |{A, C, G, T}|
        // N_a in docs
        final double[] effectiveCountsOfAllele = new double[NUM_POSSIBLE_ALLELES];

        // TODO: we don't have a good way of adding up columns of a 2-dimensional array
        for (int a = 0; a < NUM_POSSIBLE_ALLELES; a++){
            effectiveCountsOfAllele[a] = MathUtils.sum(effectiveCountsGivenAllele[a]);
        }

        final int totalCount = (int) MathUtils.sum(effectiveCounts);
        assert totalCount == numExamples : "effectiveCount should add up to numExamples";

        // K-dimensional vector of sample means weighted by the responsibilities, \bar{m} in the docs
        final double[] weightedSampleMeanM = GATKProtectedMathUtils.sumArrayFunction(0, numExamples,
                i -> MathArrays.scale((double) data.getAltDepths().get(i), responsibilities[i]));

        // K-dimensional vector of sample means weighted by the responsibilities, \bar{x} in the docs
        final double[] weightedSampleMeanX = GATKProtectedMathUtils.sumArrayFunction(0, numExamples,
                i -> MathArrays.scale((double) data.getAltF1R2Depths().get(i), responsibilities[i]));

        // K-dimensional vector of mean read depths weighted by the responsibilities, \bar{R} in the docs
        final double[] weightedDepthR = GATKProtectedMathUtils.sumArrayFunction(0, numExamples,
                i -> MathArrays.scale((double) data.getDepths().get(i), responsibilities[i]));

        theta = MathArrays.ebeDivide(weightedSampleMeanX, weightedSampleMeanM);
        f = MathArrays.ebeDivide(weightedSampleMeanM, weightedDepthR);

        for (int a = 0; a < NUM_POSSIBLE_ALLELES; a++){
            pi[a] = MathArrays.scale(1/effectiveCountsOfAllele[a], effectiveCountsGivenAllele[a]);
        }
    }

    private boolean checkLikelihoodHasConverged(final double oldLikelihood, final double newLikelihood){
        return Math.abs(newLikelihood - oldLikelihood) < CONVERGENCE_THRESHOLD;
    }

    private static boolean checkLikelihoodHasConverged(final int numIterations){
        return numIterations > 10;
    }

    class Hyperparameters {
        double[][] pi;
        double[] f;
        double[] theta;

        public Hyperparameters(double[][] pi, double[] f, double[] theta){
            this.pi = pi;
            this.f = f;
            this.theta = theta;
        }

        double[][] getPi(){ return pi; }
        double[] getF(){ return f; }
        double[] getTheta(){ return theta; }
    }

}
