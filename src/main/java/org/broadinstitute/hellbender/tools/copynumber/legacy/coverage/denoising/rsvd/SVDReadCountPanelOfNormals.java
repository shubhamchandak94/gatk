package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.rsvd;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.DenoisedCopyRatioResult;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.ReadCountPanelOfNormals;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;

/**
 * Interface for the panel of normals for SVD-based coverage denoising.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public interface SVDReadCountPanelOfNormals extends ReadCountPanelOfNormals {
    /**
     * Returns the PoN version.
     */
    double getVersion();

    /**
     * Returns a modifiable copy of the list of intervals that were used to build this PoN
     * (no filtering will have been applied).
     */
    List<SimpleInterval> getAllIntervals();

    /**
     * Returns a modifiable copy of the list of intervals contained in this PoN after all filtering has been applied.
     */
    List<SimpleInterval> getPanelIntervals();

    /**
     * Returns a modifiable copy of an array containing the median (across all samples, before filtering)
     * of the fractional coverage at each interval (in the same order as in {@link #getAllIntervals()}).
     */
    double[] getAllIntervalFractionalMedians();

    /**
     * Returns a matrix with dimensions {@code TxE}, where {@code T} is the number of panel targets (after filtering)
     * and {@code E} is the number of eigensamples, to be used for denoising.
     */
    RealMatrix getRightSingularVectors();

    /**
     * Returns a pseudoinverse matrix with dimensions {@code ExT}, where {@code E} is the number of eigensamples
     * and {@code T} is the number of panel targets (after filtering), to be used for denoising.
     */
    RealMatrix getPanelPseudoinverse();

    /**
     * Returns a modifiable copy of the list of read-count file names for samples that were used to build this PoN
     * (no filtering will have been applied).
     */
    List<String> getAllSampleFileNames();

    /**
     * Returns a modifiable copy of the list of read-count file names for samples contained in this PoN after
     * all filtering has been applied.
     */
    List<String> getPanelSampleFileNames();

    @Override
    default DenoisedCopyRatioResult denoise(final ReadCountCollection readCounts, final JavaSparkContext ctx) {
        return SVDDenoisingUtils.tangentNormalize(this, readCounts, ctx);
    }
}