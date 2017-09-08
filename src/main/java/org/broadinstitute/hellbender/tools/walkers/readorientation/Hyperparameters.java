package org.broadinstitute.hellbender.tools.walkers.readorientation;

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
import java.util.Arrays;
import java.util.List;

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

    double[] getPiForAllele(final Allele allele) {
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

    String getReferenceContext() {
        return referenceContext;
    }

    /** Reading and writing the learned hyperparameters of the model **/

    /**
     * Writing
     **/
    // TODO: I can abstract this, pretty much copied the code from PileupSummary
    private static class HyperparameterTableWriter extends TableWriter<Hyperparameters> {
        private HyperparameterTableWriter(final File output) throws IOException {
            super(output, HyperparameterTableColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final Hyperparameters hps, final DataLine dataLine) {
            // it'd be nice to set() less manually...
            // Note that allele fraction f is not allele-specific, thus the same f array will be printed
            // four times for each context
            dataLine.set(HyperparameterTableColumn.CONTEXT.toString(), hps.getReferenceContext())
                    .set(HyperparameterTableColumn.PI_A.toString(), doubleArrayToString(hps.getPi()[BaseUtils.simpleBaseToBaseIndex("A".getBytes()[0])]))
                    .set(HyperparameterTableColumn.PI_C.toString(), doubleArrayToString(hps.getPi()[BaseUtils.simpleBaseToBaseIndex("C".getBytes()[0])]))
                    .set(HyperparameterTableColumn.PI_G.toString(), doubleArrayToString(hps.getPi()[BaseUtils.simpleBaseToBaseIndex("G".getBytes()[0])]))
                    .set(HyperparameterTableColumn.PI_T.toString(), doubleArrayToString(hps.getPi()[BaseUtils.simpleBaseToBaseIndex("T".getBytes()[0])]))
                    .set(HyperparameterTableColumn.F.toString(), doubleArrayToString(hps.getF()))
                    .set(HyperparameterTableColumn.THETA.toString(), doubleArrayToString(hps.getTheta()));
        }
    }

    /** Converts a double array to a comma-separated string without brackets.
     *  Contrasts Arrays.toString, which includes brackets
     */
    private static String doubleArrayToString(final double[] xs){
        Utils.validateArg(xs.length > 0, "xs must not be an empty (uninitialized?) array");
        StringBuilder sb = new StringBuilder(String.valueOf(xs[0]));
        for (int i = 1; i < xs.length ; i++){
            sb.append("," + String.valueOf(xs[i]));
        }
        return sb.toString();
    }

    public static void writeHyperparameters(final List<Hyperparameters> hyperparameters, final File outputTable) {
        try (HyperparameterTableWriter writer = new HyperparameterTableWriter(outputTable)) {
            writer.writeAllRecords(hyperparameters);
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", outputTable), e);
        }
    }

    /**
     * Reading
     **/
    public static List<Hyperparameters> readHyperparameters(final File table) {
        try (HyperParameterTableReader reader = new HyperParameterTableReader(table)) {
            return reader.toList();
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", table), e);
        }
    }


    private static class HyperParameterTableReader extends TableReader<Hyperparameters> {
        private HyperParameterTableReader(final File table) throws IOException {
            super(table);
        }

        @Override
        protected Hyperparameters createRecord(final DataLine dataLine) {
            final String referenceContext = dataLine.get(HyperparameterTableColumn.CONTEXT);
            final double[] piA = Arrays.stream(dataLine.get(HyperparameterTableColumn.PI_A).split(","))
                    .mapToDouble(Double::parseDouble).toArray();
            final double[] piC = Arrays.stream(dataLine.get(HyperparameterTableColumn.PI_C).split(","))
                    .mapToDouble(Double::parseDouble).toArray();
            final double[] piG = Arrays.stream(dataLine.get(HyperparameterTableColumn.PI_G).split(","))
                    .mapToDouble(Double::parseDouble).toArray();
            final double[] piT = Arrays.stream(dataLine.get(HyperparameterTableColumn.PI_T).split(","))
                    .mapToDouble(Double::parseDouble).toArray();

            final double[][] pi = new double[ContextDependentArtifactFilterEngine.NUM_ALLELES][ContextDependentArtifactFilterEngine.NUM_STATUSES];
            pi[BaseUtils.simpleBaseToBaseIndex("A".getBytes()[0])] = piA;
            pi[BaseUtils.simpleBaseToBaseIndex("C".getBytes()[0])] = piC;
            pi[BaseUtils.simpleBaseToBaseIndex("G".getBytes()[0])] = piG;
            pi[BaseUtils.simpleBaseToBaseIndex("T".getBytes()[0])] = piT;

            final double[] f = Arrays.stream(dataLine.get(HyperparameterTableColumn.F).split(","))
                    .mapToDouble(Double::parseDouble).toArray();
            ;
            final double[] theta = Arrays.stream(dataLine.get(HyperparameterTableColumn.THETA).split(","))
                    .mapToDouble(Double::parseDouble).toArray();
            ;

            return new Hyperparameters(referenceContext, pi, f, theta);
        }
    }

    private enum HyperparameterTableColumn {
        CONTEXT("context"),
        PI_A("pi_A"), PI_C("pi_C"), PI_G("pi_G"), PI_T("pi_T"),
        F("allele_fraction"),
        THETA("alt_f1r2_fraction");

        private String columnName;

        HyperparameterTableColumn(final String columnName) {
            this.columnName = columnName;
        }

        @Override
        public String toString() {
            return columnName;
        }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());


    }
}

