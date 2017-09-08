package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.tools.exome.HashedListTargetCollection;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.tools.walkers.mutect.FilterMutectCalls;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Given pileup data from {@link GetPileupSummaries}, calculates the fraction of reads coming from cross-sample contamination.
 *
 * <p>
 *     The resulting contamination table is used with {@link FilterMutectCalls}.
 * </p>
 *
 * <p>This tool and GetPileupSummaries together replace GATK3's ContEst.</p>
 *
 * <p>
 *     The resulting table provides the fraction contamination, one line per sample, e.g. SampleID--TAB--Contamination.
 *     The file has no header.
 * </p>
 *
 * <h3>Example</h3>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" CalculateContamination \
 *   -I pileups.table \
 *   -O contamination.table
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Calculate contamination",
        oneLineSummary = "Calculate contamination",
        programGroup = VariantProgramGroup.class
)
@DocumentedFeature
public class CalculateContamination extends CommandLineProgram {

    private static final Logger logger = LogManager.getLogger(CalculateContamination.class);

    // we consider allele fractions in a small range around 0.5 to be heterozygous.  Beyond that range is LoH.
    public static final double MIN_HET_AF = 0.4;
    public static final double MAX_HET_AF = 0.6;

    private static final double INITIAL_CONTAMINATION_GUESS = 0.1;

    public static final double LOH_RATIO_DIFFERENCE_THRESHOLD = 0.2;
    private static final double LOH_Z_SCORE_THRESHOLD = 3.0;
    private static final double HET_CONTAMINATION_CONVERGENCE_THRESHOLD = 0.01;

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc="The input table", optional = false)
    private File inputPileupSummariesTable;

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output table", optional = false)
    private final File outputTable = null;

    private static final int CNV_SCALE = 10_000_000;

    private static final Median MEDIAN = new Median();

    public enum BiallelicGenotypes {
        HOM_REF, HET, HOM_ALT
    }

    @Override
    public Object doWork() {
        final List<PileupSummary> pileupSummaries = PileupSummary.readPileupSummaries(inputPileupSummariesTable);
        List<List<PileupSummary>> neighborhoods = splitSites(pileupSummaries);
        double hetContamination = INITIAL_CONTAMINATION_GUESS;
        double homAltContamination = INITIAL_CONTAMINATION_GUESS;
        int iteration = 0;

        while (++iteration < 10) {  // loop over contamination convergence
            final double contaminationForLambda = hetContamination;
            final List<ContaminationStats> neighborhoodStats = neighborhoods.stream()
                    .map(nbhd -> ContaminationStats.getStats(nbhd, contaminationForLambda))
                    .collect(Collectors.toList());

            final double medianHetCountRatio = MEDIAN.evaluate(neighborhoodStats.stream().mapToDouble(ContaminationStats::ratioOfActualToExpectedHets).toArray());
            final double medianHomAltCountRatio = MEDIAN.evaluate(neighborhoodStats.stream().mapToDouble(ContaminationStats::ratioOfActualToExpectedHomAlts).toArray());

            // total stats over all neighborhoods that we don't flag for loss of heterozygosity
            final ContaminationStats genomeStats = new ContaminationStats();

            for (final ContaminationStats stats : neighborhoodStats) {
                final double hetRatioDifference = stats.ratioOfActualToExpectedHets() - medianHetCountRatio;
                final double hetRatioZ = hetRatioDifference / stats.getStdOfHetCount();

                final double homRatioDifference = stats.ratioOfActualToExpectedHomAlts() - medianHomAltCountRatio;
                final double homRatioZ = homRatioDifference / stats.getStdOfHomAltCount();

                // too few hets or too many hom alts indicates LoH
                if (hetRatioDifference < -LOH_RATIO_DIFFERENCE_THRESHOLD && hetRatioZ < -LOH_Z_SCORE_THRESHOLD
                        || homRatioDifference > LOH_RATIO_DIFFERENCE_THRESHOLD && homRatioZ > LOH_Z_SCORE_THRESHOLD) {
                    logger.info(String.format("Discarding region with %d hets %d hom alts versus %.2f expected hets and " +
                            "%.2f expected hom alts due to possible loss of heterozygosity", stats.getHetCount(), stats.getHomAltCount(), stats.getExpectedHetCount(), stats.getExpectedHomAltCount()));
                } else {
                    genomeStats.increment(stats);
                }
            }

            final double newHetContamination = genomeStats.contaminationFromHets();
            final boolean converged = Math.abs(hetContamination - newHetContamination) < HET_CONTAMINATION_CONVERGENCE_THRESHOLD;
            hetContamination = newHetContamination;
            homAltContamination = genomeStats.contaminationFromHomAlts();
            logger.info(String.format("In iteration %d we estimate a contamination of %.4f based on hets and %.4f based on hom alts.", iteration, newHetContamination, homAltContamination));
            if (converged) {
                ContaminationRecord.writeContaminationTable(Arrays.asList(new ContaminationRecord(ContaminationRecord.Level.WHOLE_BAM.toString(), homAltContamination, genomeStats.standardErrorOfContaminationFromHomAlts())), outputTable);
                break;
            }
        }






        return "SUCCESS";
    }


    // split list of sites into CNV-scale-sized sublists in order to flag individual sublists for loss of heterozygosity
    private static List<List<PileupSummary>> splitSites(final List<PileupSummary> sites) {
        final List<List<PileupSummary>> result = new ArrayList<>();

        final TargetCollection<PileupSummary> tc = new HashedListTargetCollection(sites);

        int currentIndex = 0;
        while (currentIndex < tc.targetCount()) {
            final PileupSummary currentSite = tc.target(currentIndex);
            final SimpleInterval nearbyRegion = new SimpleInterval(currentSite.getContig(), currentSite.getStart(), currentSite.getEnd() + CNV_SCALE);
            final IndexRange nearbyIndices = tc.indexRange(nearbyRegion);
            final List<PileupSummary> nearbySites = IntStream.range(nearbyIndices.from, nearbyIndices.to).mapToObj(tc::target).collect(Collectors.toList());
            result.add(nearbySites);
            currentIndex = nearbyIndices.to;
        }

        return result;
    }
}
