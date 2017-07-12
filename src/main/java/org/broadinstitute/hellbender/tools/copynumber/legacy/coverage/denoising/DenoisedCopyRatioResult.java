package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;

/**
 * Represents a copy-ratio profile that has been denoised by a {@link ReadCountPanelOfNormals}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class DenoisedCopyRatioResult {
    private final ReadCountCollection standardizedProfile;
    private final ReadCountCollection denoisedProfile;

    public DenoisedCopyRatioResult(final ReadCountCollection standardizedProfile,
                                       final ReadCountCollection denoisedProfile) {
        this.standardizedProfile = standardizedProfile;
        this.denoisedProfile = denoisedProfile;
    }

    public ReadCountCollection getStandardizedProfile() {
        return standardizedProfile;
    }

    public ReadCountCollection getDenoisedProfile() {
        return denoisedProfile;
    }

    /**
     * Write results along with the command line that produced them to specified files.
     */
    public void write(final File standardizedProfileFile,
                      final File denoisedProfileFile) {
        Utils.nonNull(standardizedProfileFile);
        Utils.nonNull(denoisedProfileFile);
        writeProfile(standardizedProfileFile, standardizedProfile, "Standardized copy-ratio profile");
        writeProfile(denoisedProfileFile, denoisedProfile, "Denoised copy-ratio profile");
    }

    private void writeProfile(final File file, final ReadCountCollection profile, final String title) {
        try {
            ReadCountCollectionUtils.write(file, profile,"title = " + title);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(file, e.getMessage());
        }
    }
}
