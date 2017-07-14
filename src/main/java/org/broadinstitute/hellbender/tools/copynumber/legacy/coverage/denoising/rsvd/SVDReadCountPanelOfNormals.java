package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.rsvd;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.DenoisedCopyRatioResult;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;

import java.util.List;

/**
 * Interface for the panel of normals (PoN) for SVD-based coverage denoising.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public interface SVDReadCountPanelOfNormals {
    /**
     * Returns the PoN version.
     */
    double getVersion();

    /**
     * Returns the number of eigensamples.
     */
    int getNumEigensamples();

    /**
     * Returns a modifiable copy of the original matrix of integer read-counts used to build the PoN
     * (no filtering will have been applied).
     */
    RealMatrix getOriginalReadCounts();

    /**
     * Returns a modifiable copy of the list of the original intervals that were used to build this PoN
     * (no filtering will have been applied).
     *
     * TODO replace this with a list of SimpleIntervals.  See https://github.com/broadinstitute/gatk/issues/3246
     */
    List<Target> getOriginalIntervals();

    /**
     * Returns a modifiable copy of the list of the intervals contained in this PoN after all filtering has been applied.
     *
     * TODO replace this with a list of SimpleIntervals.  See https://github.com/broadinstitute/gatk/issues/3246
     */
    List<Target> getPanelIntervals();

    /**
     * Returns an array containing the median (across all samples, before filtering)
     * of the fractional coverage at each panel interval (in the same order as in {@link #getPanelIntervals()}).
     * This is used to standardize samples.
     */
    double[] getPanelIntervalFractionalMedians();

    /**
     * Returns the singular values of the eigensamples in decreasing order.  This array has length {@code K}.
     */
    double[] getSingularValues();

    /**
     * Returns the matrix of right-singular vectors.
     * This matrix has has dimensions {@code MxK},
     * where {@code M} is the number of panel intervals (after filtering)
     * and {@code K} is the number of eigensamples.
     * Columns are sorted by singular value in decreasing order.
     */
    double[][] getRightSingular();

    /**
     * Returns the pseudoinverse of the matrix of right-singular vectors returned by {@link #getRightSingular()}.
     * This matrix has dimensions {@code KxM},
     * where {@code K} is the number of eigensamples
     * and {@code M} is the number of panel intervals (after filtering).
     */
    double[][] getRightSingularPseudoinverse();

    default DenoisedCopyRatioResult denoise(final ReadCountCollection readCounts,
                                            final int numEigensamples,
                                            final JavaSparkContext ctx) {
        return SVDDenoisingUtils.denoise(this, readCounts, numEigensamples, ctx);
    }
}