package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.pca;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.DenoisedCopyRatioProfile;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;

/**
 * Stores the results of a tangent normalization.
 */
public final class PCATangentNormalizationResult implements DenoisedCopyRatioProfile {
    private final ReadCountCollection tangentNormalized;
    private final ReadCountCollection preTangentNormalized;

    public PCATangentNormalizationResult(final ReadCountCollection tangentNormalized, final ReadCountCollection preTangentNormalized) {
        this.tangentNormalized = tangentNormalized;
        this.preTangentNormalized = preTangentNormalized;
    }

    public ReadCountCollection getTangentNormalized() {
        return tangentNormalized;
    }

    public ReadCountCollection getPreTangentNormalized() {
        return preTangentNormalized;
    }

    @Override
    public ReadCountCollection getProfile() {
        return getTangentNormalized();
    }

    /**
     * Write results along with the command line that produced them to specified files.
     */
    public void write(final String commandLine,
                      final File tangentNormalizedFile,
                      final File preTangentNormalizedFile) {
        Utils.nonNull(commandLine);
        writeTangentNormalizedOutput(tangentNormalizedFile, commandLine);
        writePreTangentNormalizationOutput(preTangentNormalizedFile, commandLine);
    }

    /**
     * Writes the tangent-normalized counts in the format described in {@link ReadCountCollectionUtils}.
     */
    private void writeTangentNormalizedOutput(final File file, final String commandLine) {
        writeOutput(file, tangentNormalized, commandLine, "Tangent-normalized coverage profile");
    }

    /**
     * Writes the pre-tangent-normalized counts in the format described in {@link ReadCountCollectionUtils}.
     * If the output file is null, we do not write a file.
     */
    private void writePreTangentNormalizationOutput(final File file, final String commandLine) {
        writeOutput(file, preTangentNormalized, commandLine, "Pre-tangent-normalized coverage profile");
    }

    private void writeOutput(final File file, final ReadCountCollection counts, final String commandLine, final String title) {
        if (file == null) {
            return;
        }
        try {
            ReadCountCollectionUtils.write(file, counts,
                    "fileFormat = tsv",
                    "commandLine = " + commandLine,
                    "title = " + title);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(file, ex.getMessage());
        }
    }
}
