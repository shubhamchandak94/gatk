package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;

/**
 * Interface for the panel of normals for coverage denoising.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public interface ReadCountPanelOfNormals {
    int getNumEigensamples();

    DenoisedCopyRatioResult denoise(ReadCountCollection readCounts, int numberOfEigensamples, JavaSparkContext ctx);
}