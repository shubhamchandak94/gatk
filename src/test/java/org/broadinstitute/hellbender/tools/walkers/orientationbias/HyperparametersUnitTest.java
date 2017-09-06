package org.broadinstitute.hellbender.tools.walkers.orientationbias;

import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import static org.broadinstitute.hellbender.tools.walkers.orientationbias.ContextDependentArtifactFilterEngine.NUM_ALLELES;
import static org.broadinstitute.hellbender.tools.walkers.orientationbias.ContextDependentArtifactFilterEngine.NUM_STATUSES;
import static org.testng.Assert.*;

/**
 * Created by tsato on 9/6/17.
 */
public class HyperparametersUnitTest {
    private static final double EPSILON = 1e-4;

    @Test
    public void test() throws IOException {
        final double[][] pi = new double[NUM_ALLELES][NUM_STATUSES];
        for (int a = 0; a < NUM_ALLELES; a++){
            Arrays.fill(pi[a], 1.0/NUM_STATUSES);
        }

        final double[] f1 = new double[NUM_STATUSES];
        Arrays.fill(f1, 0.1);

        final double[] f2 = new double[NUM_STATUSES];
        Arrays.fill(f2, 0.2);

        final double[] f3 = new double[NUM_STATUSES];
        Arrays.fill(f3, 0.3);

        final double[] theta1 = new double[NUM_STATUSES];
        Arrays.fill(theta1, 0.9);

        final double[] theta2 = new double[NUM_STATUSES];
        Arrays.fill(theta2, 0.8);

        final double[] theta3 = new double[NUM_STATUSES];
        Arrays.fill(theta3, 0.7);

        List<Hyperparameters> hyperParametersBefore = Arrays.asList(
                new Hyperparameters("ACT", pi, f1, theta1),
                new Hyperparameters("GTA", pi, f2, theta2),
                new Hyperparameters("CCC", pi, f3, theta3));
        final File table = File.createTempFile("hyperparameters", ".tsv");

        Hyperparameters.writeHyperparameters(hyperParametersBefore, table);
        List<Hyperparameters> hyperparametersAfter = Hyperparameters.readHyperparameters(table);

        Hyperparameters hyp1 = hyperparametersAfter.get(0);
        for (int a = 0; a < NUM_ALLELES; a++){
            ArrayAsserts.assertArrayEquals(hyp1.getPi()[a], pi[a], EPSILON);
        }
        ArrayAsserts.assertArrayEquals(hyp1.getF(), f1, EPSILON);

        Hyperparameters hyp2 = hyperparametersAfter.get(1);
        for (int a = 0; a < NUM_ALLELES; a++){
            ArrayAsserts.assertArrayEquals(hyp2.getPi()[a], pi[a], EPSILON);
        }
        ArrayAsserts.assertArrayEquals(hyp2.getF(), f2, EPSILON);

        Hyperparameters hyp3 = hyperparametersAfter.get(2);
        for (int a = 0; a < NUM_ALLELES; a++){
            ArrayAsserts.assertArrayEquals(hyp3.getPi()[a], pi[a], EPSILON);
        }
        ArrayAsserts.assertArrayEquals(hyp3.getF(), f3, EPSILON);
    }
}