package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.svd.SVDReadCountPanelOfNormals;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Represents a copy-ratio profile that has been standardized and denoised by an {@link SVDReadCountPanelOfNormals}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyRatioSegmentationResult {
    private final List<SimpleInterval> segments;
    private final List<Double> columnNames;
    private final RealMatrix standardizedProfile;
    private final RealMatrix denoisedProfile;

    class CopyRatioSegment {
        private final SimpleInterval interval;
        private final double meanLog2CopyRatio;

        CopyRatioSegment(final SimpleInterval interval,
                         final List<Double> log2DenoisedCopyRatios) {
            Utils.nonNull(interval);
            Utils.nonEmpty(log2DenoisedCopyRatios);
            this.interval = interval;
            this.meanLog2CopyRatio = log2DenoisedCopyRatios.stream().mapToDouble(Double::doubleValue).average().getAsDouble();
        }
    }

    public CopyRatioSegmentationResult(final List<CopyRatioSegment> segments) {
        Utils.nonEmpty(segments);
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
