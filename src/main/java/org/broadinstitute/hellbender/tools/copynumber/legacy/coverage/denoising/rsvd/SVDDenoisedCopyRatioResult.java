package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.rsvd;

import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Represents a copy-ratio profile that has been standardized and denoised by an {@link SVDReadCountPanelOfNormals}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SVDDenoisedCopyRatioResult {
    //TODO clean this up once Targets are removed and new ReadCountCollection is in place
    private final List<Target> intervals;
    private final List<String> columnNames;
    private final RealMatrix standardizedProfile;
    private final RealMatrix denoisedProfile;

    public SVDDenoisedCopyRatioResult(final List<Target> intervals,
                                      final List<String> sampleNames,
                                      final RealMatrix standardizedProfile,
                                      final RealMatrix denoisedProfile) {
        Utils.nonEmpty(intervals);
        Utils.nonEmpty(sampleNames);
        Utils.nonNull(standardizedProfile);
        Utils.nonNull(denoisedProfile);
        Utils.validateArg(intervals.size() == standardizedProfile.getColumnDimension(),
                "Number of intervals and columns in standardized profile must match.");
        Utils.validateArg(intervals.size() == denoisedProfile.getColumnDimension(),
                "Number of intervals and columns in denoised profile must match.");
        Utils.validateArg(sampleNames.size() == standardizedProfile.getRowDimension(),
                "Number of sample names and rows in standardized profile must match.");
        Utils.validateArg(sampleNames.size() == denoisedProfile.getRowDimension(),
                "Number of sample names and rows in denoised profile must match.");
        this.intervals = intervals;
        this.columnNames = sampleNames;
        this.standardizedProfile = standardizedProfile;
        this.denoisedProfile = denoisedProfile;
    }

    public void write(final File standardizedProfileFile,
                      final File denoisedProfileFile) {
        Utils.nonNull(standardizedProfileFile);
        Utils.nonNull(denoisedProfileFile);
        writeProfile(standardizedProfileFile, standardizedProfile, "Standardized copy-ratio profile");
        writeProfile(denoisedProfileFile, denoisedProfile, "Denoised copy-ratio profile");
    }

    private void writeProfile(final File file, final RealMatrix profile, final String title) {
        try {
            final ReadCountCollection rcc = new ReadCountCollection(intervals, columnNames, profile.transpose());
            ReadCountCollectionUtils.write(file, rcc,"title = " + title);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(file, e.getMessage());
        }
    }
}
