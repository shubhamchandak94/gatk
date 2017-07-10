package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising;

import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;

/**
 * Interface for a copy-ratio profile that has been denoised by a {@link CoveragePanelOfNormals}.
 */
@FunctionalInterface
public interface DenoisedCopyRatioProfile {
    ReadCountCollection getProfile();
}
