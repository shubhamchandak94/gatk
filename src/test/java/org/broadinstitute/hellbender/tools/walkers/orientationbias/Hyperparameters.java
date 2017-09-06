package org.broadinstitute.hellbender.tools.walkers.orientationbias;

import org.apache.commons.math3.util.Pair;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import htsjdk.variant.variantcontext.Allele;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;

import static org.broadinstitute.hellbender.tools.walkers.orientationbias.ContextDependentArtifactFilterEngine.ALL_ALLELES;
import static org.broadinstitute.hellbender.tools.walkers.orientationbias.ContextDependentArtifactFilterEngine.States;

/**
 * Created by tsato on 9/5/17.
 */
class Hyperparameters {
    String referenceContext;
    double[][] pi;
    double[] f;
    double[] theta;

    public Hyperparameters(final String referenceContext, final double[][] pi, final double[] f, final double[] theta) {
        this.referenceContext = referenceContext;
        this.pi = pi;
        this.f = f;
        this.theta = theta;
    }

    double[] getPiForAllele(final Allele allele){
        Utils.validateArg(allele.getBases().length == 1, "allele has to be single base");
        return pi[BaseUtils.simpleBaseToBaseIndex(allele.getBases()[0])];
    }

    double[][] getPi() {
        return pi;
    }

    double[] getF() {
        return f;
    }

    double[] getTheta() {
        return theta;
    }

    String getReferenceContext() { return referenceContext; }

    /** Reading and writing the learned hyperparameters of the model **/

    /** Writing **/
    public static void writeHyperparameters(final Collection<Hyperparameters> hyperparameters, final File outputTable) {
        try ( HyperparameterTableWriter writer = new HyperparameterTableWriter(outputTable)) {
            for (Hyperparameters hps : hyperparameters){
                // each hyperparameter object contains four rows worth of information (one row per allele)
                // because each call to writeRecord() prints out one row, we cannot
                for (String allele : ALL_ALLELES) {
                    writer.writeRecord(new Pair(hyperparameters, Allele.create(allele)));
                }
            }
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", outputTable), e);
        }

    }

    // Perhaps I can abstract this, pretty much copied the code from PileupSummary
    private static class HyperparameterTableWriter extends TableWriter<Pair<Hyperparameters, Allele>> {
        private HyperparameterTableWriter(final File output) throws IOException {
            super(output, HyperparameterTableColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final Pair<Hyperparameters, Allele> pair, final DataLine dataLine) {
            final Hyperparameters hps = pair.getFirst();
            final Allele allele = pair.getSecond();
            final double[] pi = hps.getPiForAllele(allele);
            final double[] f = hps.getF();

            // it'd be nice to set() less manually...
            // Note that allele fraction f is not allele-specific, thus the same f array will be printed
            // four times for each context
            dataLine.set(HyperparameterTableColumn.CONTEXT.toString(), hps.getReferenceContext())
                    .set(HyperparameterTableColumn.ALLELE.toString(), allele.toString())
                    .set(HyperparameterTableColumn.PI_F1R2.toString(), pi[States.F1R2.ordinal()])
                    .set(HyperparameterTableColumn.PI_F2R1.toString(), pi[States.F2R1.ordinal()])
                    .set(HyperparameterTableColumn.PI_HOMREF.toString(), pi[States.BALANCED_HOM_REF.ordinal()])
                    .set(HyperparameterTableColumn.PI_HET.toString(), pi[States.BALANCED_HET.ordinal()])
                    .set(HyperparameterTableColumn.PI_HOMVAR.toString(), pi[States.BALANCED_HOM_VAR.ordinal()])
                    .set(HyperparameterTableColumn.F_F1R2.toString(), f[States.F1R2.ordinal()])
                    .set(HyperparameterTableColumn.F_F2R1.toString(), f[States.F2R1.ordinal()])
                    .set(HyperparameterTableColumn.F_HOMREF.toString(), f[States.BALANCED_HOM_REF.ordinal()])
                    .set(HyperparameterTableColumn.F_HET.toString(), f[States.BALANCED_HET.ordinal()])
                    .set(HyperparameterTableColumn.F_HOMVAR.toString(), f[States.BALANCED_HOM_VAR.ordinal()]);
            // TODO: potentially add thetas
        }
    }

    /** Reading **/
    public static List<Hyperparameters> readHyperparameters(final File tableFile) {
        try( HyperParametemrTableReader reader = new HyperParametemrTableReader(tableFile) ) {
            return reader.toList();
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", tableFile));
        }
    }


    private static class HyperParametemrTableReader extends TableReader<Pair<Hyperparameters, Allele>> {
        @Override
        protected Pair<Hyperparameters, Allele> createRecord(final DataLine dataLine) {
            final String referenceContext = dataLine.get(HyperparameterTableColumn.CONTEXT);
            final Allele allele = Allele.create(dataLine.get(HyperparameterTableColumn.ALLELE));
            final double[][] pi;
            final double[] f;

            new Hyperparameters(referenceContext, pi, f, theta);
        }
    }

    public PileupSummaryTableReader(final File file) throws IOException { super(file); }

    @Override
    protected PileupSummary createRecord(final DataLine dataLine) {
        final String contig = dataLine.get(PileupSummaryTableColumn.CONTIG);
        final int position = dataLine.getInt(PileupSummaryTableColumn.POSITION);
        final int refCount = dataLine.getInt(PileupSummaryTableColumn.REF_COUNT);
        final int altCount = dataLine.getInt(PileupSummaryTableColumn.ALT_COUNT);
        final int otherAltCount = dataLine.getInt(PileupSummaryTableColumn.OTHER_ALT_COUNT);
        final double alleleFrequency = dataLine.getDouble(PileupSummaryTableColumn.ALT_ALLELE_FREQUENCY);

        return new PileupSummary(contig, position, refCount, altCount, otherAltCount, alleleFrequency);
    }
}

    private enum HyperparameterTableColumn {
        CONTEXT("context"),
        ALLELE("allele"),
        PI_F1R2("pi_f1r2"), PI_F2R1("pi_f2r1"), PI_HOMREF("pi_homref"), PI_HET("pi_het"), PI_HOMVAR("pi_homvar"),
        F_F1R2("f_f1r2"), F_F2R1("f_f2r1"), F_HOMREF("f_homref"), F_HET("f_het"), F_HOMVAR("f_homvar");
        // add \theta as needed

        private String columnName;

        HyperparameterTableColumn(final String columnName){ this.columnName = columnName; }

        @Override
        public String toString() { return columnName; }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());


    }

}

