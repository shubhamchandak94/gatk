package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.rsvd;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.MatrixSummaryUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.svd.SVD;
import org.broadinstitute.hellbender.utils.svd.SVDFactory;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

/**
 * This class contains utility methods for creating and writing the PCA coverage panel of normals to an {@link HDF5RandomizedSVDReadCountPanelOfNormals}.
 */
public final class HDF5PCACoveragePoNCreationUtils {
    public static final double EPSILON = 1E-9;
    private static final double INV_LN_2 = 1.0 / Math.log(2);

    private static final Logger logger = LogManager.getLogger(HDF5PCACoveragePoNCreationUtils.class);

    private HDF5PCACoveragePoNCreationUtils() {}

    /**
     * Create an HDF5 coverage PoN file with the output from {@link CombineReadCounts}.
     *
     * When using this method in testing, see {@link CreatePanelOfNormals} for good default values.
     * @param ctx If no spark context is available, specify {@code null}
     * @param outputHDF5Filename  final filename for the output PoN HDF5 file
     * @param openMode            desired {@link HDF5File.OpenMode} (if {@code HDF5File.OpenMode.READ_ONLY}, an exception will be thrown if not a dry run)
     * @param inputPCovFile  pcov file from {@link CombineReadCounts} with samples as columns.  Never {@code null}
     * @param initialTargets  the set of targets that was used to generate the {@code inputPCovFile}.
     * @param sampleNameBlacklist  sample names in {@code inputPCovFile} that should be ignored.  Names must match exactly.  Never {@code null}.  Entires in this parameter that are not in {@code inputPCovFile} are ignored.
     * @param targetFactorPercentileThreshold  the percentile of extreme target factor values to cut.
     * @param maximumPercentageZeroColumns  the maximum percentage of zero values in a sample (across targets) before the sample is filtered.
     * @param maximumPercentageZeroTargets  the maximum percentage of zero values in a target (across samples) before the target is filtered.
     * @param extremeColumnMedianCountPercentileThreshold Percentile to cut columns after they are ranked by median value
     * @param countTruncatePercentile percentile (on either end) to truncate extreme values (i.e. extreme high or low values will be set to the high or low count truncation value, which is based on this percentlie)
     * @param numberOfEigensamples  desired number of eigensamples in the final PoN.  If missing, will apply a rule to determine.
     * @param isDryRun  whether this is a dry run.  If you wish to do diagnostics and skip the actual writing of a PoN file, set this to true.
     */
    public static void create(final JavaSparkContext ctx,
                              final File outputHDF5Filename,
                              final HDF5File.OpenMode openMode,
                              final File inputPCovFile,
                              final TargetCollection<Target> initialTargets,
                              final List<String> sampleNameBlacklist,
                              final double targetFactorPercentileThreshold,
                              final double maximumPercentageZeroColumns,
                              final double maximumPercentageZeroTargets,
                              final double extremeColumnMedianCountPercentileThreshold,
                              final double countTruncatePercentile,
                              final OptionalInt numberOfEigensamples,
                              final boolean isDryRun) {
        Utils.nonNull(outputHDF5Filename);
        IOUtils.canReadFile(inputPCovFile);
        Utils.nonNull(initialTargets, "Target collection cannot be null.");
        Utils.nonNull(sampleNameBlacklist, "Blacklist sample list cannot be null.  Use empty list if no blacklisting is desired.");
        ParamUtils.inRange(targetFactorPercentileThreshold, 0, 100, "Target factor percentile threshold must be in range [0, 100].");
        ParamUtils.inRange(maximumPercentageZeroColumns, 0, 100, "Maximum percentage of zero-columns must be in range [0, 100].");
        ParamUtils.inRange(maximumPercentageZeroTargets, 0, 100, "Maximum percentage of zero-targets must be in range [0, 100].");
        ParamUtils.inRange(extremeColumnMedianCountPercentileThreshold, 0, 50, "Extreme column median percentile threshold must be in range [0, 50].");
        ParamUtils.inRange(countTruncatePercentile, 0, 50, "Count truncation threshold percentile threshold must be in range [0, 50].");
        Utils.nonNull(numberOfEigensamples, "Number of eigensamples cannot be null.");

        final ReadCountCollection inputPCov = readReadCountsFromFile(inputPCovFile, initialTargets);    //ignore missing targets

        final Pair<ReadCountCollection, double[]> inputSubsetByUsableTargets = subsetReadCountsToUsableTargets(inputPCov, targetFactorPercentileThreshold, logger);

        // Remove samples listed on the blacklist and normalize read-counts by removing the target factor component
        ReadCountCollection normalizedCounts = inputSubsetByUsableTargets.getLeft();
        normalizedCounts = normalizedCounts.subsetColumns(Sets.difference(new HashSet<>(normalizedCounts.columnNames()), new HashSet<>(sampleNameBlacklist)) );
        final double[] targetFactors = inputSubsetByUsableTargets.getLeft();
        SVDDenoisingUtils.factorNormalize(normalizedCounts.counts(), targetFactors);

        // Impute zeros as median values for columns and targets with too many zeros, remove targets with extreme medians, and truncate extreme counts
        final ReadCountCollection logNormalizedCounts = cleanNormalizedCounts(normalizedCounts, logger, maximumPercentageZeroColumns, maximumPercentageZeroTargets, extremeColumnMedianCountPercentileThreshold, countTruncatePercentile);

        // Normalize by the median and log_2 scale the read counts.
        normalizeAndLogReadCounts(logNormalizedCounts, logger);

        // Subtract the median of medians
        subtractMedianOfMedians(logNormalizedCounts, logger);

        // Perform the SVD and calculate the pseudoinverse
        final ReductionResult reduction = calculateReducedPanelAndPInverses(logNormalizedCounts, numberOfEigensamples, logger, ctx);

        // Calculate the target variances
        final List<String> panelTargetNames = logNormalizedCounts.targets().stream().map(Target::getName).collect(Collectors.toList());
        final double[] targetVariances = calculateTargetVariances(normalizedCounts, panelTargetNames, reduction, ctx);

        // Write the PoN to HDF5 file
        if (!isDryRun) {
            HDF5RandomizedSVDReadCountPanelOfNormals.write(outputHDF5Filename, openMode, initialTargets.targets(), normalizedCounts, logNormalizedCounts, targetFactors, targetVariances, reduction);
        }
    }



    /**
     * Performs
     * {@link ReadCountCollectionUtils#removeColumnsWithTooManyZeros},
     * {@link ReadCountCollectionUtils#removeTargetsWithTooManyZeros},
     * {@link ReadCountCollectionUtils#removeColumnsWithExtremeMedianCounts},
     * {@link ReadCountCollectionUtils#imputeZeroCountsAsTargetMedians}, and
     * {@link ReadCountCollectionUtils#truncateExtremeCounts}.
     *
     * The input {@link ReadCountCollection} should contain proportional coverage.
     */
    @VisibleForTesting
    static ReadCountCollection cleanNormalizedCounts(final ReadCountCollection readCounts, final Logger logger,
                                                     double maximumPercentageZeroColumns,
                                                     double maximumPercentageZeroTargets,
                                                     double extremeColumnMedianCountPercentileThreshold,
                                                     double countTruncatePercentile) {
        // Remove column and targets with too many zeros and targets with extreme median coverage:
        final int maximumColumnZerosCount = calculateMaximumZerosCount(readCounts.targets().size(), maximumPercentageZeroColumns);
        final int maximumTargetZerosCount = calculateMaximumZerosCount(readCounts.columnNames().size(), maximumPercentageZeroTargets);
        ReadCountCollection cleanedCounts = ReadCountCollectionUtils.removeColumnsWithTooManyZeros(readCounts, maximumColumnZerosCount, false, logger);
        cleanedCounts = ReadCountCollectionUtils.removeTargetsWithTooManyZeros(cleanedCounts, maximumTargetZerosCount, false, logger);
        cleanedCounts = ReadCountCollectionUtils.removeColumnsWithExtremeMedianCounts(cleanedCounts, extremeColumnMedianCountPercentileThreshold, logger);

        // Impute zero counts to be same as the median for the same target.
        // This happens in-place.
        ReadCountCollectionUtils.imputeZeroCountsAsTargetMedians(cleanedCounts, logger);

        // Truncate extreme count values.
        // This happens in-place.
        ReadCountCollectionUtils.truncateExtremeCounts(cleanedCounts, countTruncatePercentile, logger);
        return cleanedCounts;
    }

//    /**
//     * Determine the variance for each target in the PoN (panel targets).
//     *
//     * @return      array of doubles where each double corresponds to a target in the PoN (panel targets)
//     */
//    private static double[] calculateTargetVariances(final ReadCountCollection normalizedCounts,
//                                                     final List<String> panelTargetNames,
//                                                     final ReductionResult reduction,
//                                                     final JavaSparkContext ctx) {
//        Utils.nonNull(panelTargetNames);
//        Utils.nonNull(normalizedCounts);
//        Utils.nonNull(reduction);
//        final SVDDenoisedCopyRatioProfile allNormals =
//                SVDDenoisingUtils.tangentNormalizeNormalsInPoN(normalizedCounts, panelTargetNames, reduction.getReducedCounts(), reduction.getReducedPseudoinverse(), ctx);
//        final RealMatrix allSampleProjectedTargets = allNormals.getTangentNormalized().counts();
//
//        return MatrixSummaryUtils.getRowVariances(allSampleProjectedTargets);
//    }

    /**
     * Calculates the absolute number of zeroes tolerated given the user requested percentage
     * and the number of counts involved.
     *
     * @param totalCounts   a non-negative number of counts to check
     * @param percentage    a value in [0, 100]
     * @return              a value between 0 and {@code totalCounts}
     */
    private static int calculateMaximumZerosCount(final int totalCounts, final double percentage) {
        ParamUtils.isPositiveOrZero(totalCounts, "Total counts must be non-negative.");
        ParamUtils.inRange(percentage, 0, 100, "Percentage must be in [0, 100].");
        return (int) Math.ceil(totalCounts * percentage / 100.0);
    }
}
