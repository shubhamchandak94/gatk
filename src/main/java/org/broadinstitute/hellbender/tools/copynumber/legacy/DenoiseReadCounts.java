package org.broadinstitute.hellbender.tools.copynumber.legacy;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5Library;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.DenoisedCopyRatioResult;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.ReadCountPanelOfNormals;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.rsvd.HDF5RandomizedSVDCoveragePoN;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.IOException;
import java.util.stream.IntStream;

/**
 * Denoises read counts given the panel of normals (PoN) created by {@link CreateReadCountPanelOfNormals} to produce
 * a copy-ratio profile.
 *
 * <h3>Examples</h3>
 * <p>
 *     The following command is for either whole exome sequencing (WES) or whole genome sequencing (WGS) data.
 * </p>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" DenoiseReadCounts \
 *   --input tumor.coverage.tsv \
 *   --panelOfNormals panel_of_normals.pon \
 *   --tangentNormalized tumor.tn.tsv \
 *   --preTangentNormalized tumor.preTN.tsv
 * </pre>
 *
 * The resulting copy-ratio profile is log2 transformed.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Denoise read counts using a panel of normals",
        oneLineSummary = "Denoise read counts using a panel of normals",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class DenoiseReadCounts extends SparkCommandLineProgram {
    @Argument(
            doc = "Input read-count file containing integer read counts in genomic intervals for a single case sample.",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME
    )
    protected File inputReadCountsFile;

    @Argument(
            doc = "Input HDF5 file containing the panel of normals (output of CreateReadCountPanelOfNormals).",
            fullName = ExomeStandardArgumentDefinitions.PON_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME
    )
    protected File inputPoNFile;

    @Argument(
            doc = "Output file for standardized (pre-tangent-normalized) copy-ratio profile.",
            fullName = ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME
    )
    protected File standardizedProfileFile;

    @Argument(
            doc = "Output file for fully denoised (tangent-normalized) copy-ratio profile.",
            fullName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME
    )
    protected File denoisedProfileFile;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        if (! new HDF5Library().load(null)){ //Note: passing null means using the default temp dir.
            throw new UserException.HardwareFeatureException("Cannot load the required HDF5 library. " +
                    "HDF5 is currently supported on x86-64 architecture and Linux or OSX systems.");
        }

        validateInputFiles();

        try (final HDF5File hdf5PoNFile = new HDF5File(inputPoNFile)) {  //HDF5File implements AutoCloseable
            //load input files
            final ReadCountCollection rcc = ReadCountCollectionUtils.parse(inputReadCountsFile);
            final ReadCountPanelOfNormals pon = new HDF5RandomizedSVDCoveragePoN(hdf5PoNFile, logger);

            //check that read-count collection contains single sample and integer counts
            validateReadCounts(rcc);

            //perform denoising and write result
            final DenoisedCopyRatioResult denoisedCopyRatioResult = pon.denoise(rcc, ctx);
            denoisedCopyRatioResult.write(standardizedProfileFile, denoisedProfileFile);

            logger.info("Read counts successfully denoised.");
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(inputReadCountsFile);
        }
    }

    private void validateInputFiles() {
        IOUtils.canReadFile(inputPoNFile);
        IOUtils.canReadFile(inputReadCountsFile);
    }

    //TODO remove these checks once ReadCountCollection is refactored to only store single sample, non-negative integer counts
    private void validateReadCounts(final ReadCountCollection readCountCollection) {
        if (readCountCollection.columnNames().size() != 1) {
            throw new UserException.BadInput("Read-count file must contain counts for only a single sample.");
        }
        if (readCountCollection.targets().isEmpty()) {
            throw new UserException.BadInput("Read-count file must contain counts for at least one genomic interval.");
        }
        final double[] readCounts = readCountCollection.counts().getColumn(0);
        if (!IntStream.range(0, readCounts.length).allMatch(i -> (readCounts[i] >= 0) &&((int) readCounts[i] == readCounts[i]))) {
            throw new UserException.BadInput("Read-count file must contain non-negative integer counts.");
        }
    }
}
