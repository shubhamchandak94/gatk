package org.broadinstitute.hellbender.tools.walkers.orientationbias;

import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.spark_project.guava.annotations.VisibleForTesting;

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
    static final int NUM_ALLELES = ALL_ALLELES.length;

    // When the increase in likelihood falls under this value, we call the algorithm converged
    static final double CONVERGENCE_THRESHOLD = 1e-3;

    // Regularizer (?) TODO: think this through
    static final double EPSILON = 1e-4;

    String referenceContext;

    // Observed data
    final PerContextData data;

    // N by K matrix of (natural) log posterior probabilities of latent variable z
    // with the current estimates of hyperparameters pi, f, and theta
    final double[][] log10Responsibilities;

    // experimental measure - eventually pick a log or linear version
    final double[][] responsibilities;

    final int numExamples;

    // A by K matrix of effective counts. Each column must add up to N_k
    // TODO: rename to CountsPerAllele?
    final double[][] effectiveCountsGivenAllele = new double[NUM_ALLELES][NUM_STATUSES];

    // K-dimensional vector of effective sample counts for each class of z, weighted by the the responsibilities. For a fixed k,
    // we sum up the counts over all alleles. N_k in the docs.
    // TODO: should be final
    @VisibleForTesting
    double[] effectiveCounts = new double[NUM_STATUSES];

    // A-dimensional vector of effective count of samples for each allele, where A = 4 = |{A, C, G, T}|, N_a in docs
    @VisibleForTesting
    double[] effectiveCountsOfAlleles = new double[NUM_ALLELES];

    /*** Hyperparameters of the model ***/

    // pi is the A by K matrix of weights (prior?) for the categorical variable z
    // For each allele a, the K-dimensional vector pi_a must be a valid probability distribution over z
    // In other words, each row must add up to 1
    final double[][] pi = new double[NUM_ALLELES][NUM_STATUSES];

    // K-dimensional vector of probabilities for the binomial m given z, which represents the number of alt reads at site n
    final double[] f = new double[NUM_STATUSES];

    // In case we observe no data - assume that the allele fraction given artifact is this value
    static final double DEFAULT_ARTIFACT_ALLELE_FRACTION = 0.3;

    // K-dimensional vector of probabilities for the binomial x given z, which represents the number of F1R2 alt reads at site n
    final double[] theta = new double[NUM_STATUSES];

    private int numIterations = 0;

    private static int MAX_ITERATIONS = 200;

    // one may plot the changes in L2 distance of parameters to make sure that EM is steadily moving towards the (local? global?) maximum
    // TODO: this being public is questionable
    public double[] l2distancesOfParameters = new double[MAX_ITERATIONS];


    public ContextDependentArtifactFilterEngine(final PerContextData data){
        this.data = data;
        numExamples = data.getNumLoci();
        log10Responsibilities = new double[numExamples][NUM_STATUSES];
        responsibilities = new double[numExamples][NUM_STATUSES];

        // initialize responsibilities in log 10 space
        final double initialLog10Probability = - Math.log10(NUM_STATUSES);
        for (int n = 0; n < data.getNumLoci(); n++ ) {
            Arrays.fill(log10Responsibilities[n], initialLog10Probability);
        }

        final double initialProbability = - 1.0/ NUM_STATUSES;
        for (int n = 0; n < data.getNumLoci(); n++ ) {
            Arrays.fill(responsibilities[n], initialProbability);
        }

        referenceContext = data.getReferenceContext();

        // populate alleles
        final String referenceBase = referenceContext.substring(1, 2);

        // we fix some of the parameters to entice the model to assign particular states to the indices into z
        // for instance, we fix the allele fraction parameter f for z = Balanced Het to be 0.5.
        f[States.BALANCED_HET.ordinal()] = 0.5;
        f[States.BALANCED_HOM_REF.ordinal()] = EPSILON;
        f[States.BALANCED_HOM_VAR.ordinal()] = 1 - EPSILON;

        // similarly, we may fix some of theta_z
        // TODO: or should we learn them?
        theta[States.F1R2.ordinal()] = 1 - EPSILON;
        theta[States.F2R1.ordinal()] = EPSILON;
        theta[States.BALANCED_HOM_REF.ordinal()] = 0.5;
        theta[States.BALANCED_HET.ordinal()] = 0.5;
        theta[States.BALANCED_HOM_VAR.ordinal()] = 0.5;
    }

    public Hyperparameters runEMAlgorithm(){
        boolean converged = false;
        double[][] oldPi = new double[NUM_ALLELES][NUM_STATUSES];

        for (int a = 0; a < NUM_ALLELES; a++){
             oldPi[a] = Arrays.copyOf(pi[a], NUM_STATUSES);
        }


        while (!converged && numIterations < MAX_ITERATIONS){
            // TODO: stylistic problems here, there's too much side-effect stuff
            takeMstep();

            // assert newLikelihood >= oldLikelihood : "M step must increase the likelihood";
            final double l2Distance = IntStream.range(0, NUM_ALLELES).mapToDouble(a -> MathArrays.distance(oldPi[a], pi[a])).sum();
            converged = l2Distance < CONVERGENCE_THRESHOLD;

            l2distancesOfParameters[numIterations] = l2Distance;

            for (int a = 0; a < NUM_ALLELES; a++){
                oldPi[a] = Arrays.copyOf(pi[a], NUM_STATUSES);
            }

            takeEstep();

            numIterations++;
        }

        return new Hyperparameters(pi, f, theta);
    }

    // Given the current estimates of the parameters pi, f, and theta, compute the log10Responsibilities
    // gamma_nk = p(z_nk|a)
    private void takeEstep(){
        for (int example = 0; example < numExamples; example++){
            final int n = example;

            Allele allele = data.getAlleles().get(n);
            final int depth = data.getDepths().get(n);
            final short altDepth = data.getAltDepths().get(n);
            final short altF1R2Depth = data.getAltF1R2Depths().get(n);

            assert allele.getBases().length == 1 : "the length of allele must be 1";
            final int a = BaseUtils.simpleBaseToBaseIndex(allele.getBases()[0]); // hack to work around the fact that java stream doesn't let you use a non-final variable

            final double[] log10AlleleFractionTerms = IntStream.range(0, NUM_STATUSES)
                    .mapToDouble(k -> new BinomialDistribution(depth, f[k]).logProbability(altDepth) * MathUtils.LOG10_OF_E)
                    .toArray();

            // When we have no alt reads (i.e. m = 0), the binomial over the number of F1R2 is a deterministic;
            // namely, it's always 1, and log(1) = 0
            final double[] logAltF1R2FractionTerms = altDepth == 0 ? new double[NUM_STATUSES] :
                    IntStream.range(0, NUM_STATUSES)
                            .mapToDouble(k -> new BinomialDistribution(altDepth, theta[k]).logProbability(altF1R2Depth) * MathUtils.LOG10_OF_E)
                            .toArray();

            assert log10AlleleFractionTerms.length == NUM_STATUSES : "alleleFractionFactors must have length K";
            assert logAltF1R2FractionTerms.length == NUM_STATUSES : "altF1R2FractionFactors must have length K";

            // TODO: might want to do this in log space, watch out for underflow.
            double[] log10UnnormalizedResponsibilities = IntStream.range(0, NUM_STATUSES)
                    .mapToDouble(k -> Math.log10(pi[a][k]) + log10AlleleFractionTerms[k] + logAltF1R2FractionTerms[k])
                    .toArray();

            log10Responsibilities[n] = MathUtils.normalizeLog10(log10UnnormalizedResponsibilities, true, false);
            responsibilities[n] = MathUtils.normalizeFromLog10ToLinearSpace(log10UnnormalizedResponsibilities);

            assert Math.abs(MathUtils.sumLog10(log10Responsibilities[n]) - 1.0) < EPSILON :
                    String.format("log responsibility for %dth example added up to %f", n,  MathUtils.sumLog10(log10Responsibilities[n]));
            assert Math.abs(MathUtils.sum(responsibilities[n]) - 1.0) < EPSILON :
                    String.format("responsibility for %dth example added up to %f", n,  MathUtils.sumLog10(responsibilities[n]));
        }
    }

    // given the current posterior distributions over z, compute the maximum likelihood estimate for
    // the categorical weights (pi), allele fractions (f), and alt F1R2 fraction (theta)
    private void takeMstep(){
        /*** compute responsibility-based statistics based on the current log10Responsibilities ***/

        // reset the effectiveCountsGivenAllele array
        for (int a = 0; a < NUM_ALLELES; a++){
            Arrays.fill(effectiveCountsGivenAllele[a], 0.0);
        }

        // TODO: optimize
        for (int n = 0; n < numExamples; n++) {
            Allele allele = data.getAlleles().get(n);
            final int a = BaseUtils.simpleBaseToBaseIndex(allele.getBases()[0]);

            // Warning; raising the log responsibilities by 10 may be naive REWORD
            effectiveCountsGivenAllele[a] = MathArrays.ebeAdd(effectiveCountsGivenAllele[a],
                    Arrays.stream(log10Responsibilities[n]).map(logp -> Math.pow(10, logp)).toArray());
        }

        // For a fixed k, sum up the effective counts over all alleles. N_k in the docs.
        effectiveCounts = GATKProtectedMathUtils.sumArrayFunction(0, NUM_ALLELES, a -> effectiveCountsGivenAllele[a]);

        assert effectiveCounts.length == NUM_STATUSES : "effectiveCount must be a k-dimensional vector";
        assert Math.abs(MathUtils.sum(effectiveCounts) - numExamples) < EPSILON : String.format("effective counts must add up to number of examples but got %f", MathUtils.sum(effectiveCounts));

        // TODO: we don't have a good way of adding up columns of a 2-dimensional array
        for (int a = 0; a < NUM_ALLELES; a++){
            effectiveCountsOfAlleles[a] = MathUtils.sum(effectiveCountsGivenAllele[a]);
        }

        assert Math.abs(MathUtils.sum(effectiveCountsOfAlleles) - numExamples) < EPSILON : "effectiveCountOfAlleles should add up to numExamples";

        // K-dimensional vector of sample means weighted by the log10Responsibilities, \bar{m} in the docs
        // TOOD: I should probabily not enter log space - where do we need it? Where do we multiply a lot of small numbers?
        final double[] weightedSampleMeanM = GATKProtectedMathUtils.sumArrayFunction(0, numExamples,
                n -> MathArrays.scale((double) data.getAltDepths().get(n), responsibilities[n]));

        // K-dimensional vector of sample means weighted by the log10Responsibilities, \bar{x} in the docs
        final double[] weightedSampleMeanX = GATKProtectedMathUtils.sumArrayFunction(0, numExamples,
                n -> MathArrays.scale((double) data.getAltF1R2Depths().get(n), responsibilities[n]));

        // K-dimensional vector of mean read depths weighted by the log10Responsibilities, \bar{R} in the docs
        final double[] weightedDepthR = GATKProtectedMathUtils.sumArrayFunction(0, numExamples,
                n -> MathArrays.scale((double) data.getDepths().get(n), responsibilities[n]));

        // TODO: should we learn theta?
        // theta = MathArrays.ebeDivide(weightedSampleMeanX, weightedSampleMeanM);

        // We update some allele fractions according to data and log10Responsibilities and keep others fixed
        final double[] updatedAlleleFractions = MathArrays.ebeDivide(weightedSampleMeanM, weightedDepthR);
        f[States.F1R2.ordinal()] = updatedAlleleFractions[States.F1R2.ordinal()] == 0 ? DEFAULT_ARTIFACT_ALLELE_FRACTION :
                updatedAlleleFractions[States.F1R2.ordinal()];
        f[States.F2R1.ordinal()] = updatedAlleleFractions[States.F2R1.ordinal()]  == 0 ? DEFAULT_ARTIFACT_ALLELE_FRACTION :
                updatedAlleleFractions[States.F2R1.ordinal()];

        for (int a = 0; a < NUM_ALLELES; a++){
            final double numExamplesWithThisAllele = effectiveCountsOfAlleles[a];
            pi[a] = numExamplesWithThisAllele < EPSILON ? new double[NUM_STATUSES] :
                    MathArrays.scale(1/numExamplesWithThisAllele, effectiveCountsGivenAllele[a]);
        }

        for (int a = 0; a < NUM_ALLELES; a++){
            final double sumProbabilities = Math.abs(MathUtils.sum(pi[a]));
            assert sumProbabilities - 1.0 < EPSILON || sumProbabilities == 0.0 : "pi[a] must add up to 1";
        }

        return;
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

    enum States {
        F1R2, // Orientation bias at the site, evidence for alt is predominantly F1R2 reads
        F2R1, // Orientation Bias at the site, evidence for alt is predominantly F1R1 reads
        BALANCED_HOM_REF, // No orientation bias, and the site is hom ref
        BALANCED_HET, // No orientation bias, and the site is het
        BALANCED_HOM_VAR // No orientation bias, and the site is hom var
    }

}
