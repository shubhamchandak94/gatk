package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.distribution.BinomialDistribution;
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
import org.broadinstitute.hellbender.utils.MathUtils;
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
    public static final double P_VALUE_THRESHOLD_FOR_HETS = 0.4;

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc="The input table", optional = false)
    private File inputPileupSummariesTable;

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output table", optional = false)
    private final File outputTable = null;

    private static final int CNV_SCALE = 1000000;

    @Override
    public Object doWork() {
        final List<PileupSummary> pileupSummaries = PileupSummary.readPileupSummaries(inputPileupSummariesTable);
        final List<PileupSummary> homAltSites = findConfidentHomAltSites(pileupSummaries);
        final Pair<Double, Double> contaminationAndError = homAltSites.isEmpty() ? Pair.of(0.0, 0.0) : calculateContamination(homAltSites);
        ContaminationRecord.writeContaminationTable(Arrays.asList(new ContaminationRecord(ContaminationRecord.Level.WHOLE_BAM.toString(), contaminationAndError.getLeft(), contaminationAndError.getRight())), outputTable);

        return "SUCCESS";
    }

    private Pair<Double, Double> calculateContamination(List<PileupSummary> homAltSites) {

        final long totalReadCount = homAltSites.stream().mapToLong(PileupSummary::getTotalCount).sum();
        final long totalRefCount = homAltSites.stream().mapToLong(PileupSummary::getRefCount).sum();

        // if eg ref is A, alt is C, then # of ref reads due to error is roughly (# of G read + # of T reads)/2
        final long errorRefCount = homAltSites.stream().mapToLong(PileupSummary::getOtherAltCount).sum() / 2;
        final long contaminationRefCount = Math.max(totalRefCount - errorRefCount, 0);
        final double totalDepthWeightedByRefFrequency = homAltSites.stream()
                .mapToDouble(ps -> ps.getTotalCount() * (1 - ps.getAlleleFrequency()))
                .sum();
        final double contamination = contaminationRefCount / totalDepthWeightedByRefFrequency;
        final double standardError = Math.sqrt(contamination / totalDepthWeightedByRefFrequency);

        logger.info(String.format("In %d homozygous variant sites we find %d reference reads due to contamination and %d" +
                        " due to to sequencing error out of a total %d reads.", homAltSites.size(), contaminationRefCount, errorRefCount, totalReadCount));
        logger.info(String.format("Based on population data, we would expect %d reference reads in a contaminant with equal depths at these sites.", (long) totalDepthWeightedByRefFrequency));
        logger.info(String.format("Therefore, we estimate a contamination of %.3f.", contamination));
        logger.info(String.format("The error bars on this estimate are %.5f.", standardError));
        return Pair.of(contamination, standardError);
    }

    private static List<PileupSummary> findConfidentHomAltSites(List<PileupSummary> sites) {
        if (sites.isEmpty()) {
            return new ArrayList<>();
        }

        // the probability of a hom alt is f^2
        final double expectedNumberOfHomAlts = sites.stream()
                .mapToDouble(PileupSummary::getAlleleFrequency).map(MathUtils::square).sum();

        // the variance in the Bernoulli count with hom alt probability p = f^2 is p(1-p)
        final double stdOfNumberOfHomAlts = Math.sqrt(sites.stream()
                .mapToDouble(PileupSummary::getAlleleFrequency).map(MathUtils::square).map(x -> x*(1-x)).sum());

        logger.info(String.format("We expect %.3f +/- %.3f hom alts", expectedNumberOfHomAlts, stdOfNumberOfHomAlts));
        final TargetCollection<PileupSummary> tc = new HashedListTargetCollection<>(sites);
        final double averageCoverage = sites.stream().mapToInt(PileupSummary::getTotalCount).average().getAsDouble();

        final List<PileupSummary> potentialHomAltSites = sites.stream()
                .filter(s -> s.getAltFraction() > 0.8)
                .collect(Collectors.toList());

        logger.info(String.format("We find %d potential hom alt sites", potentialHomAltSites.size()));

        final List<PileupSummary> filteredHomAltSites = new ArrayList<>();
        for (final PileupSummary site : potentialHomAltSites) {
            logger.info(String.format("Considering hom alt site %s:%d.", site.getContig(), site.getStart()));
            final SimpleInterval nearbySpan = new SimpleInterval(site.getContig(), Math.max(1, site.getStart() - CNV_SCALE), site.getEnd() + CNV_SCALE);
            final List<PileupSummary> nearbySites = tc.targets(nearbySpan);

            final double averageNearbyCopyRatio = nearbySites.stream().mapToDouble(s -> s.getTotalCount()/averageCoverage).average().orElseGet(() -> 0);
            logger.info(String.format("The average copy ratio in the vicinity of this site is %.2f", averageNearbyCopyRatio));
            final double expectedNumberOfNearbyConfidentHets = (1 - P_VALUE_THRESHOLD_FOR_HETS) * nearbySites.stream().mapToDouble(PileupSummary::getAlleleFrequency).map(x -> 2*x*(1-x)).sum();
            final long numberOfNearbyConfidentHets = nearbySites.stream().filter(ps -> isConfidentHet(ps, P_VALUE_THRESHOLD_FOR_HETS)).count();
            logger.info(String.format("We expect %.1f hets near here and found %d", expectedNumberOfNearbyConfidentHets, numberOfNearbyConfidentHets));
            if (numberOfNearbyConfidentHets > 0.5 * expectedNumberOfHomAlts) {
                if (averageNearbyCopyRatio > 0.6 && averageNearbyCopyRatio < 3.0) {
                    filteredHomAltSites.add(site);
                } else {
                    logger.info("We reject this site due to anomalous copy ratio");
                }
            } else {
                logger.info("We reject this site due to potential loss of heterozygosity.");
            }

            //TODO: as extra security, filter out sites that are near too many hom alts

        }

        logger.info(String.format("We excluded %d candidate hom alt sites.", potentialHomAltSites.size() - filteredHomAltSites.size()));

        return filteredHomAltSites;
    }

    // Can we reject the null hypothesis that a site is het?
    // the pValue threshold is (almost) the rejection rate of true hets -- almost because the binomial is a discrete distribution
    private static boolean isConfidentHet(final PileupSummary ps, final double pValueThreshold) {
        final int altCount = ps.getAltCount();
        final int totalCount = ps.getTotalCount();
        final BinomialDistribution binomialDistribution = new BinomialDistribution(null, totalCount, 0.5);

        // this pValue is the probability under the null hypothesis of an allele fraction at least as far from 1/2 as this
        // the factor of 2 comes from the symmetry of the binomial distribution -- the two-sided p-value is twice the one-side p-value
        final double pValue = 2 * binomialDistribution.cumulativeProbability(Math.min(altCount, totalCount - altCount));

        return pValue > pValueThreshold;
    }
}
