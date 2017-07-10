package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;

/**
 * Interface for the panel of normals for coverage denoising.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public interface CoveragePanelOfNormals<T extends DenoisedCopyRatioProfile> {
    T denoise(final ReadCountCollection readCounts, final JavaSparkContext ctx);
}