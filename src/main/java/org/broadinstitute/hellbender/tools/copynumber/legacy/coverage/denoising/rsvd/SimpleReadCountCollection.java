package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.rsvd;

import htsjdk.samtools.util.Locatable;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.supercsv.cellprocessor.ParseInt;
import org.supercsv.cellprocessor.constraint.NotNull;
import org.supercsv.cellprocessor.ift.CellProcessor;
import org.supercsv.comment.CommentStartsWith;
import org.supercsv.io.CsvListReader;
import org.supercsv.io.ICsvListReader;
import org.supercsv.prefs.CsvPreference;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * //TODO replace this class with updated ReadCountCollection
 *
 * Simple class to enable input and output of a TSV containing a list of {@link Locatable} and integer read counts.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SimpleReadCountCollection {
    private static final Logger logger = LogManager.getLogger(SimpleReadCountCollection.class);

    private static final CsvPreference TAB_SKIP_COMMENTS_PREFERENCE = new CsvPreference.Builder(CsvPreference.TAB_PREFERENCE)
            .skipComments(new CommentStartsWith("#")).build();

    private final List<Locatable> intervals;
    private final RealMatrix readCounts;

    private SimpleReadCountCollection(final List<Locatable> intervals,
                                      final RealMatrix readCounts) {
        Utils.nonEmpty(intervals);
        Utils.nonNull(readCounts);
        Utils.validateArg(readCounts.getColumnDimension() == 1, "Read-count matrix must contain only a single column.");
        Utils.validateArg(intervals.size() == readCounts.getRowDimension(), "Number of intervals and read counts must match.");
        Utils.validateArg(Arrays.stream(readCounts.getColumn(0)).noneMatch(c -> c < 0), "Read counts must all be non-negative integers.");
        Utils.validateArg(intervals.stream().distinct().count() == intervals.size(), "Intervals must all be unique.");

        this.intervals = intervals;
        this.readCounts = readCounts;
    }

    public List<Locatable> getIntervals() {
        return intervals;
    }

    public RealMatrix getReadCounts() {
        return readCounts;
    }

    public static SimpleReadCountCollection read(final File file) {
        IOUtils.canReadFile(file);
        final int numLines = countLines(file);  //this includes comment lines and column headers
        final List<Locatable> intervals = new ArrayList<>(numLines);
        final List<Integer> readCounts = new ArrayList<>(numLines);
        try (final FileReader fileReader = new FileReader(file);
             final ICsvListReader listReader = new CsvListReader(fileReader, TAB_SKIP_COMMENTS_PREFERENCE)) {
            listReader.getHeader(true);
            final CellProcessor[] processors = new CellProcessor[]{
                    new NotNull(),  //contig
                    new ParseInt(), //start
                    new ParseInt(), //stop
                    null,           //ignore name
                    new ParseInt()  //count
            };

            List<Object> row;
            while ((row = listReader.read(processors)) != null) {
                final Locatable interval = new SimpleInterval(row.get(0).toString(), (int) row.get(1), (int) row.get(2));
                intervals.add(interval);
                readCounts.add((int) row.get(4));
            }
            final RealMatrix readCountsMatrix = new Array2DRowRealMatrix(readCounts.size(), 1);
            readCountsMatrix.setColumn(0, readCounts.stream().mapToDouble(Integer::doubleValue).toArray());
            return new SimpleReadCountCollection(intervals, readCountsMatrix);
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(file);
        }
    }

    private static int countLines(final File file) {
        try {
            return (int) Files.lines(file.toPath()).count();
        } catch (final IOException e) {
            throw new UserException.BadInput("Could not determine number of lines in TSV file.");
        }
    }
}
