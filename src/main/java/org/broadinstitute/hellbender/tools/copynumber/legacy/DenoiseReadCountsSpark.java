package org.broadinstitute.hellbender.tools.copynumber.legacy;

import org.apache.logging.log4j.Level;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5Library;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.pca.HDF5PCACoveragePoN;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.pca.PCACoveragePoN;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.pca.PCATangentNormalizationResult;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Normalizes read counts given the PanelOfNormals (PoN).
 *
 * <p> A note to developers:  If this is extended to use Spark, please be wary that the parallelization in tangent normalization is
 * by case sample, which may not yield benefits for most use cases (which are one sample)  </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 *
 * <h3>Examples</h3>
 * <p>
 *     The following command is for either whole exome sequencing (WES) or whole genome sequencing (WGS) data.
 * </p>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" NormalizeSomaticReadCounts \
 *   --input tumor.coverage.tsv \
 *   --panelOfNormals panel_of_normals.pon \
 *   --tangentNormalized tumor.tn.tsv \
 *   --preTangentNormalized tumor.preTN.tsv
 * </pre>
 *
 * The resulting data are log2 transformed. Currently, the tool can only produce log2-transformed counts.
 */
@CommandLineProgramProperties(
        summary = "Normalize PCOV read counts using a panel of normals",
        oneLineSummary = "Normalize proportional coverage (PCOV) read counts using a panel of normals",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class DenoiseReadCountsSpark extends SparkCommandLineProgram {
    @Argument(
            doc = "read counts input file.  This can only contain one sample.",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME
    )
    protected File readCountsFile;

    @Argument(
            doc = "panel of normals HDF5 file",
            shortName = ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.PON_FILE_LONG_NAME,
            optional = false
    )
    protected File ponFile;

    @Argument(
            doc = "Tangent normalized counts output",
            shortName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            optional = false
    )
    protected File tangentNormalizationOutFile;

    @Argument(
            doc = "Pre-tangent normalization counts",
            shortName = ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            optional = true
    )
    protected File preTangentNormalizationOutFile;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        if (! new HDF5Library().load(null)){ //Note: passing null means using the default temp dir.
            throw new UserException.HardwareFeatureException("Cannot load the required HDF5 library. " +
                    "HDF5 is currently supported on x86-64 architecture and Linux or OSX systems.");
        }
        IOUtils.canReadFile(ponFile);
        try (final HDF5File hdf5PoNFile = new HDF5File(ponFile)) {  //HDF5File implements AutoCloseable
            final PCACoveragePoN pon = new HDF5PCACoveragePoN(hdf5PoNFile, logger);
            final ReadCountCollection readCounts = ReadCountCollectionUtils.parse(readCountsFile);
            final PCATangentNormalizationResult tangentNormalizationResult = pon.denoise(readCounts, ctx);
            tangentNormalizationResult.write(getCommandLine(), tangentNormalizationOutFile, preTangentNormalizationOutFile);
            logger.info("Read counts successfully denoised.");
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(readCountsFile);
        }
    }
}
