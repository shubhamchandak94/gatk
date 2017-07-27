package org.broadinstitute.hellbender.tools.copynumber.temporary;

import htsjdk.samtools.util.Lazy;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * NOTE: There is code duplication in here, since some of the older code is going to be removed in the future.
 *
 * HDF5 paths are chosen to reflect that these files will ultimately only store data for a single sample,
 * but code supports multiple samples.
 */
public class HDF5ReadCountCollection {
    private static final String INTERVALS_GROUP_NAME = "/intervals";
    private static final String SAMPLE_NAME_PATH = "/sample_name/value";
    private static final String READ_COUNTS_PATH = "/read_counts/values";

    private final HDF5File file;
    private Lazy<List<Locatable>> intervals;
    private Lazy<List<String>> sampleNames;

    public HDF5ReadCountCollection(final HDF5File file) {
        Utils.nonNull(file, "The input file cannot be null.");
        this.file = file;
        intervals = new Lazy<>(() -> HDF5Utils.readIntervals(file, INTERVALS_GROUP_NAME));
        sampleNames = new Lazy<>(() -> readNames(file, SAMPLE_NAME_PATH));
    }

    public List<Locatable> getIntervals() {
        return intervals.get();
    }

    public List<String> getSampleNames() {
        return sampleNames.get();
    }

    public RealMatrix getReadCounts() {
        final double[][] values = file.readDoubleMatrix(READ_COUNTS_PATH);
        if (values.length != sampleNames.get().size()) {
            throw new GATKException(String.format("Wrong number of elements in the read counts recovered " +
                            "from file '%s': number of read counts found in file (%d) != number of samples (%d).",
                    file.getFile(), values.length, sampleNames.get().size()));
        }
        if (values[0].length != intervals.get().size()) {
            throw new GATKException(String.format("Wrong number of elements in the read counts recovered " +
                            "from file '%s': number of read counts found in file (%d) != number of intervals (%d).",
                    file.getFile(), values[0].length, intervals.get().size()));
        }
        return new Array2DRowRealMatrix(values);
    }

    public List<Target> getTargets() {
        return IntStream.range(0, intervals.get().size()).boxed()
                .map(i -> new Target((SimpleInterval) intervals.get().get(i))).collect(Collectors.toList());
    }

    /**
     *  Create (or modify) HDF5 file.
     *
     * @param outFile Path to the final HDF5 file.  Not {@code null}
     * @param targets the intervals and names for each target.  Not {@code null}
     * @param values a T x S matrix where T is the number of targets and S is the number of samples.  Not {@code null}
     * @param sampleNames the column names for each of the value columns.  Not {@code null}
     */
    public static void write(final File outFile,
                             final List<Target> targets, final double[][] values, final List<String> sampleNames) {
        Utils.nonNull(outFile);
        Utils.nonNull(targets);
        Utils.nonNull(values);
        Utils.nonNull(sampleNames);

        if (values.length != sampleNames.size()) {
            throw new GATKException("The shape of the values array (" + values.length + " x " + values[0].length + ") does not match the number of samples (" + sampleNames.size() + ").");
        }
        if (values[0].length != targets.size()) {
            throw new GATKException("The shape of the values array (" + values.length + " x " + values[0].length + ") does not match the number of targets (" + targets.size() + ").");
        }

        try (final HDF5File file = new HDF5File(outFile, HDF5File.OpenMode.CREATE)) {
            final HDF5ReadCountCollection hdf5ReadCountCollection = new HDF5ReadCountCollection(file);
            hdf5ReadCountCollection.writeIntervals(targets);
            hdf5ReadCountCollection.writeNames(SAMPLE_NAME_PATH, sampleNames);
            hdf5ReadCountCollection.writeReadCounts(values);
        }
    }

    private static List<String> readNames(final HDF5File reader, final String namesPath) {
        final String[] values = reader.readStringArray(namesPath);
        return Collections.unmodifiableList(Arrays.asList(values));
    }

    private <T extends Locatable> void writeIntervals(final List<T> intervals) {
        HDF5Utils.writeIntervals(file, INTERVALS_GROUP_NAME, intervals);
    }

    private void writeNames(final String path, final List<String> names) {
        file.makeStringArray(path, names.toArray(new String[names.size()]));
    }

    private void writeReadCounts(final double[][] readCounts) {
        file.makeDoubleMatrix(READ_COUNTS_PATH, readCounts);
    }
}
