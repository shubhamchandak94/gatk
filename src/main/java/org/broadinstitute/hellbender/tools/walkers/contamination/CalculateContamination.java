package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.apache.commons.lang.mutable.MutableDouble;
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
import org.broadinstitute.hellbender.utils.IndexRange;
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

    // we consider allele fractions in a small range around 0.5 to be heterozygous.  Beyond that range is LoH.
    private static final double MIN_HET_AF = 0.4;
    private static final double MAX_HET_AF = 0.6;

    private static final double INITIAL_CONTAMINATION_GUESS = 0.1;

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

    private enum BiallelicGenotypes {
        HOM_REF, HET, HOM_ALT
    }

    // statistics relevant to computing contamination and blacklisting possible LoH regions
    private static class Stats {
        private double homAltCount = 0;
        private double hetCount = 0;

        private double readCountInHomAltSites;
        private double refCountInHomAltSites;
        private double otherAltCountInHomAltSites;

        double homAltDepthWeightedByRefFrequency = 0; //TODO: ?????

        public void increment(final EnumMap<BiallelicGenotypes, Double> posteriors, final PileupSummary ps) {
            final double homAltResponsibility = posteriors.get(BiallelicGenotypes.HOM_ALT);
            final double hetResponsibility = posteriors.get(BiallelicGenotypes.HET);

            homAltCount += homAltResponsibility;
            this.hetCount += hetResponsibility;

            readCountInHomAltSites += homAltResponsibility * ps.getTotalCount();
            refCountInHomAltSites += homAltResponsibility * ps.getRefCount();
            otherAltCountInHomAltSites += homAltResponsibility * ps.getOtherAltCount();
        }

        public static Stats getStats(final Collection<PileupSummary> pileupSummaries, final double contamination) {
            final Stats result = new Stats();
            pileupSummaries.forEach(ps -> result.increment(genotypePosteriors(ps, contamination), ps));
            return result;
        }

        public void increment(final Stats other) {
            this.homAltCount += other.homAltCount;
            this.hetCount += other.hetCount;

            this.readCountInHomAltSites += other.readCountInHomAltSites;
            this.refCountInHomAltSites += other.refCountInHomAltSites;
            this.otherAltCountInHomAltSites += other.otherAltCountInHomAltSites;

            //TODO: as Stats expands, need to increment new members
        }

        public boolean isLossOfHeterozygosity() {
            return false;
            //TODO: this is a stub

    }

    @Override
    public Object doWork() {
        final List<PileupSummary> pileupSummaries = PileupSummary.readPileupSummaries(inputPileupSummariesTable);

        List<List<PileupSummary>> neighborhoods = splitSites(pileupSummaries);
        double contamination = INITIAL_CONTAMINATION_GUESS;

        while (true) {  // loop over contamination convergence
            final Stats genomeStats = new Stats();

            neighborhoods.stream()
                    .map( nbhd -> Stats.getStats(nbhd, contamination))
                    .filter(stats -> !stats.isLossOfHeterozygosity())
                    .forEach(genomeStats::increment);

            }




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

        final double expectedNumberOfHomAlts = expectedHomAltCount(sites);
        final double stdOfNumberOfHomAlts = standardDeviationOfHomAltCount(sites);

        final TargetCollection<PileupSummary> allSites = new HashedListTargetCollection<>(sites);
        final double medianCoverage = MEDIAN.evaluate(sites.stream().mapToDouble(PileupSummary::getTotalCount).toArray());

        final TargetCollection<PileupSummary> potentialHomAltSites = new HashedListTargetCollection<>(sites.stream()
                .filter(s -> s.getAltFraction() > 0.8)
                .collect(Collectors.toList()));

        logger.info(String.format("We expect %.1f +/- %.1f and find %d hom alts.", expectedNumberOfHomAlts, stdOfNumberOfHomAlts, potentialHomAltSites.targetCount()));

        final List<PileupSummary> filteredHomAltSites = new ArrayList<>();
        for (final PileupSummary site : potentialHomAltSites.targets()) {

            final SimpleInterval nearbySpan = new SimpleInterval(site.getContig(), Math.max(1, site.getStart() - CNV_SCALE), site.getEnd() + CNV_SCALE);
            final List<PileupSummary> nearbySites = allSites.targets(nearbySpan);
            final int localHomAltCount = potentialHomAltSites.targets(nearbySpan).size();
            final double expectedLocalHomAltCount = expectedHomAltCount(nearbySites);
            final long localHetCount =
            final double expectedLocalHetCount = expectedHetCount(nearbySites);
            final double localCopyRatio = MEDIAN.evaluate(nearbySites.stream().mapToDouble(s -> s.getTotalCount()).toArray()) / medianCoverage;

            final boolean tooFewHets = localHetCount < 0.5 * expectedLocalHetCount;
            final boolean tooManyHomAlts = localHomAltCount > 3 * expectedLocalHomAltCount && localHomAltCount > expectedLocalHetCount / 4;
            if (tooFewHets && tooManyHomAlts) {
                logger.info(String.format("Rejecting candidate hom alt site %s:%d due to suspected loss of heterozygosity.  " +
                                "We expect %.1f and observe %d hets nearby, while we expect %.1f and observe %d hom alts nearby.", site.getContig(), site.getStart(),
                        expectedLocalHetCount, localHetCount, expectedLocalHomAltCount, localHomAltCount));
            } else if (localCopyRatio < 0.6 || localCopyRatio > 3.0) {
                logger.info(String.format("We reject this site due to anomalous copy ratio %.3f", localCopyRatio));
            } else {
                filteredHomAltSites.add(site);
            }
        }

        logger.info(String.format("We excluded %d candidate hom alt sites.", potentialHomAltSites.targetCount() - filteredHomAltSites.size()));

        return filteredHomAltSites;
    }

    // the probability of a hom alt is f^2
    private static double expectedHetCount(List<PileupSummary> sites) {
        return sites.stream().mapToDouble(PileupSummary::getAlleleFrequency).map(x -> 2 * x * (1 - x)).sum();
    }

    // the probability of a hom alt is f^2
    private static double expectedHomAltCount(List<PileupSummary> sites) {
        return sites.stream().mapToDouble(PileupSummary::getAlleleFrequency).map(MathUtils::square).sum();
    }

    // the variance in the Bernoulli count with hom alt probability p = f^2 is p(1-p)
    private static double standardDeviationOfHomAltCount(List<PileupSummary> sites) {
        return Math.sqrt(sites.stream().mapToDouble(PileupSummary::getAlleleFrequency).map(MathUtils::square).map(x -> x*(1-x)).sum());
    }



    // contamination is a current rough estimate of contamination
    private static EnumMap<BiallelicGenotypes, Double> genotypePosteriors(final PileupSummary ps, final double contamination) {
        final double alleleFrequency = ps.getAlleleFrequency();
        final double homRefPrior = MathUtils.square(1 - alleleFrequency);
        final double hetPrior = 2 * alleleFrequency * (1 - alleleFrequency);
        final double homAltPrior = MathUtils.square(alleleFrequency);

        final int totalCount = ps.getTotalCount();
        final int altCount = ps.getAltCount();

        final double maxHomRefFraction = contamination;
        final double minHomAltFraction = 1 - contamination;
        final double minHetFraction = MIN_HET_AF - contamination / 2;
        final double maxHetFraction = MAX_HET_AF + contamination / 2;

        final double homRefLikelihood = MathUtils.uniformBinomialProbability(totalCount, altCount, 0, maxHomRefFraction);
        final double homAltLikelihood = MathUtils.uniformBinomialProbability(totalCount, altCount, minHomAltFraction, 1);
        final double hetLikelihood = MathUtils.uniformBinomialProbability(totalCount, altCount, minHetFraction, maxHetFraction);

        final double[] unnormalized = new double[] {homRefLikelihood * homRefPrior, hetLikelihood * hetPrior, homAltLikelihood * homAltPrior};
        final double[] normalized = MathUtils.normalizeFromRealSpace(unnormalized, true);

        final EnumMap<BiallelicGenotypes, Double> result = new EnumMap<BiallelicGenotypes, Double>(BiallelicGenotypes.class);
        result.put(BiallelicGenotypes.HOM_REF, normalized[0]);
        result.put(BiallelicGenotypes.HET, normalized[1]);
        result.put(BiallelicGenotypes.HOM_ALT, normalized[2]);

        return result;
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
