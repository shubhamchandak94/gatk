package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.Collection;
import java.util.EnumMap;

// statistics relevant to computing contamination and blacklisting possible LoH regions
class ContaminationStats {
    public double getHomAltCount() {
        return homAltCount;
    }

    public double getHetCount() {
        return hetCount;
    }

    private double homAltCount;
    private double hetCount;

    public double getExpectedHomAltCount() {
        return expectedHomAltCount;
    }

    public double getExpectedHetCount() {
        return expectedHetCount;
    }

    private double expectedHomAltCount;
    private double expectedHetCount;

    public double getStdOfHomAltCount() {
        return Math.sqrt(varianceOfHomAltCount);
    }

    public double getStdOfHetCount() {
        return Math.sqrt(varianceOfHetCount);
    }

    private double varianceOfHomAltCount;
    private double varianceOfHetCount;

    private double readCountInHomAltSites;
    private double refCountInHomAltSites;
    private double otherAltCountInHomAltSites;

    private double refCountInHetSites;
    private double altCountInHetSites;

    private double expectedRefInHomAltPerUnitContamination;
    private double expectedRefExcessInHetPerUnitContamination;

    public void increment(final EnumMap<CalculateContamination.BiallelicGenotypes, Double> posteriors, final PileupSummary ps) {
        final double homAltResponsibility = posteriors.get(CalculateContamination.BiallelicGenotypes.HOM_ALT);
        final double hetResponsibility = posteriors.get(CalculateContamination.BiallelicGenotypes.HET);

        homAltCount += homAltResponsibility;
        hetCount += hetResponsibility;

        final double homAltPrior = MathUtils.square(ps.getAlleleFrequency());
        final double hetPrior = 2 * ps.getAlleleFrequency() * ( 1 - ps.getAlleleFrequency());
        expectedHomAltCount += homAltPrior;
        expectedHetCount += hetPrior;
        varianceOfHomAltCount += homAltPrior * (1 - homAltPrior);
        varianceOfHetCount += hetPrior * (1 - hetPrior);

        readCountInHomAltSites += homAltResponsibility * ps.getTotalCount();
        refCountInHomAltSites += homAltResponsibility * ps.getRefCount();
        otherAltCountInHomAltSites += homAltResponsibility * ps.getOtherAltCount();

        refCountInHetSites += hetResponsibility * ps.getRefCount();
        altCountInHetSites += hetResponsibility * ps.getAltCount();

        expectedRefInHomAltPerUnitContamination += homAltResponsibility * ps.getTotalCount() * ps.getRefFrequency();
        expectedRefExcessInHetPerUnitContamination += hetResponsibility * ps.getTotalCount() * ( ps.getRefFrequency() - ps.getAlleleFrequency());
    }

    public void increment(final ContaminationStats other) {
        this.homAltCount += other.homAltCount;
        this.hetCount += other.hetCount;

        this.expectedHomAltCount = other.expectedHomAltCount;
        this.expectedHetCount += other.expectedHetCount;
        this.varianceOfHomAltCount += other.varianceOfHomAltCount;
        this.varianceOfHetCount += other.varianceOfHetCount;

        this.readCountInHomAltSites += other.readCountInHomAltSites;
        this.refCountInHomAltSites += other.refCountInHomAltSites;
        this.otherAltCountInHomAltSites += other.otherAltCountInHomAltSites;

        this.refCountInHetSites += other.refCountInHetSites;
        this.altCountInHetSites += other.altCountInHetSites;

        this.expectedRefInHomAltPerUnitContamination += other.expectedRefInHomAltPerUnitContamination;
        this.expectedRefExcessInHetPerUnitContamination += other.expectedRefExcessInHetPerUnitContamination;
    }

    public double contaminationFromHomAlts() {
        // if ref is A, alt is C, then # of ref reads due to error is roughly (# of G read + # of T reads)/2
        final double refInHomAltDueToError = otherAltCountInHomAltSites / 2;
        final double refCountInHomAltDueToContamination = Math.max(refCountInHomAltSites - refInHomAltDueToError, 0);
        return refCountInHomAltDueToContamination / expectedRefInHomAltPerUnitContamination;
    }

    public double standardErrorOfContaminationFromHomAlts() {
        return Math.sqrt(contaminationFromHomAlts() / expectedRefInHomAltPerUnitContamination);
    }

    public double contaminationFromHets() {
        final double refExcessInHetSites = refCountInHetSites - altCountInHetSites;
        return Math.max(refExcessInHetSites / expectedRefExcessInHetPerUnitContamination, 0);
    }

    public boolean isLossOfHeterozygosity() {
        return false;
        //TODO: this is a stub
    }

    public double ratioOfActualToExpectedHets() {
        return hetCount / expectedHetCount;
    }

    public double ratioOfActualToExpectedHomAlts() {
        return homAltCount / expectedHomAltCount;
    }

    public static ContaminationStats getStats(final Collection<PileupSummary> pileupSummaries, final double contamination) {
        final ContaminationStats result = new ContaminationStats();
        pileupSummaries.forEach(ps -> result.increment(genotypePosteriors(ps, contamination), ps));
        return result;
    }

    // contamination is a current rough estimate of contamination
    private static EnumMap<CalculateContamination.BiallelicGenotypes, Double> genotypePosteriors(final PileupSummary ps, final double contamination) {
        final double alleleFrequency = ps.getAlleleFrequency();
        final double homRefPrior = MathUtils.square(1 - alleleFrequency);
        final double hetPrior = 2 * alleleFrequency * (1 - alleleFrequency);
        final double homAltPrior = MathUtils.square(alleleFrequency);

        final int totalCount = ps.getTotalCount();
        final int altCount = ps.getAltCount();

        final double homRefLikelihood = MathUtils.uniformBinomialProbability(totalCount, altCount, 0, contamination);
        final double homAltLikelihood = MathUtils.uniformBinomialProbability(totalCount, altCount, 1 - contamination, 1);
        final double hetLikelihood = MathUtils.uniformBinomialProbability(totalCount, altCount, CalculateContamination.MIN_HET_AF - contamination / 2, CalculateContamination.MAX_HET_AF + contamination / 2);

        final double[] unnormalized = new double[CalculateContamination.BiallelicGenotypes.values().length];
        unnormalized[CalculateContamination.BiallelicGenotypes.HOM_REF.ordinal()] = homRefLikelihood * homRefPrior;
        unnormalized[CalculateContamination.BiallelicGenotypes.HET.ordinal()] = hetLikelihood * hetPrior;
        unnormalized[CalculateContamination.BiallelicGenotypes.HOM_ALT.ordinal()] = homAltLikelihood * homAltPrior;
        if (Double.isNaN(unnormalized[0])) {
            int j = 2;
            double test = MathUtils.uniformBinomialProbability(totalCount, altCount, 0, contamination);
        }
        final double[] normalized = MathUtils.normalizeFromRealSpace(unnormalized, true);

        final EnumMap<CalculateContamination.BiallelicGenotypes, Double> result = new EnumMap<CalculateContamination.BiallelicGenotypes, Double>(CalculateContamination.BiallelicGenotypes.class);
        result.put(CalculateContamination.BiallelicGenotypes.HOM_REF, normalized[0]);
        result.put(CalculateContamination.BiallelicGenotypes.HET, normalized[1]);
        result.put(CalculateContamination.BiallelicGenotypes.HOM_ALT, normalized[2]);

        return result;
    }
}
