package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.pca;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.CoveragePanelOfNormals;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;

/**
 * Interface for the panel of normals for PCA coverage denoising.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public interface PCACoveragePoN extends CoveragePanelOfNormals<PCATangentNormalizationResult> {
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
     * Returns a modifiable copy of an array containing the median of the proportional coverage
     * (calculated across all samples) at each interval (in the same order as in {@link #getAllIntervals()}).
     */
    double[] getAllIntervalProportionalMedians();

    /**
     * Returns a matrix with dimensions {@code TxE}, where {@code T} is the number of panel targets (after filtering)
     * and {@code E} is the number of eigensamples, to be used for denoising.
     */
    RealMatrix getPanel();

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
    default PCATangentNormalizationResult denoise(final ReadCountCollection readCounts, final JavaSparkContext ctx) {
        return PCATangentNormalizationUtils.tangentNormalize(this, readCounts, ctx);
    }
}