package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.rsvd;

import org.apache.commons.math3.linear.RealMatrix;

/**
 * Structure to hold the result of the PoN reduction steps.
 */
final class ReductionResult {
    private final RealMatrix pseudoInverse;
    private final RealMatrix reducedCounts;
    private final RealMatrix reducedPseudoInverse;

    ReductionResult(final RealMatrix pseudoInverse, final RealMatrix reducedCounts, final RealMatrix reducedPseudoInverse) {
        this.pseudoInverse = pseudoInverse;
        this.reducedCounts = reducedCounts;
        this.reducedPseudoInverse = reducedPseudoInverse;
    }

    RealMatrix getPseudoInverse() {
        return pseudoInverse;
    }

    RealMatrix getReducedCounts() {
        return reducedCounts;
    }

    RealMatrix getReducedPseudoInverse() {
        return reducedPseudoInverse;
    }
}
