package org.broadinstitute.hellbender.tools.copynumber.legacy;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.LegacyCopyNumberArgument;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * TODO
 *
 * @author Samuel Lee &lt;davidben@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Segment a denoised copy-ratio profile into regions of constant copy ratio.",
        oneLineSummary = "Segment a denoised copy-ratio profile into regions of constant copy ratio.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class PerformCopyRatioSegmentation extends CommandLineProgram {
    @Argument(
            doc = "Input file containing denoised copy-ratio profile (output of DenoiseReadCounts).",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME
    )
    protected String inputDenoisedProfileFile;

    @Argument(
            doc = "Output file for copy-ratio segments.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    protected File outputSegmentsFile;

    @Override
    public Object doWork() {
        final String sampleName = ReadCountCollectionUtils.getSampleNameForCLIsFromReadCountsFile(new File(inputDenoisedProfileFile));
        final ReadCountCollection rcc;
        try {
            rcc = ReadCountCollectionUtils.parse(new File(inputDenoisedProfileFile));
        } catch (final IOException ex) {
            throw new UserException.BadInput("could not read input file");
        }

        return "SUCCESS";
    }
}
