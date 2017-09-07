package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.SegmentTableColumn;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

/**
 * Represents a legacy copy-ratio segmentation.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyRatioSegmentationResult {
    //TODO update column headers; we keep the old ones for now
    private static final TableColumnCollection COPY_RATIO_SEGMENT_FILE_TABLE_COLUMNS = SegmentTableColumn.MEAN_AND_NO_CALL_COLUMNS;

    private final List<CopyRatioSegment> segments;

    static class CopyRatioSegment {
        private final SimpleInterval interval;
        private final int numDataPoints;
        private final double meanLog2CopyRatio;

        /**
         * @param denoisedCopyRatios    in log2 space
         */
        CopyRatioSegment(final SimpleInterval interval,
                         final List<Double> denoisedCopyRatios) {
            Utils.nonNull(interval);
            Utils.nonEmpty(denoisedCopyRatios);
            this.interval = interval;
            numDataPoints = denoisedCopyRatios.size();
            meanLog2CopyRatio = denoisedCopyRatios.stream().mapToDouble(Double::doubleValue).average().getAsDouble();
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final CopyRatioSegment that = (CopyRatioSegment) o;
            if (numDataPoints != that.numDataPoints) {
                return false;
            }
            if (Double.compare(that.meanLog2CopyRatio, meanLog2CopyRatio) != 0) {
                return false;
            }
            return interval.equals(that.interval);
        }

        @Override
        public int hashCode() {
            int result;
            long temp;
            result = interval.hashCode();
            result = 31 * result + numDataPoints;
            temp = Double.doubleToLongBits(meanLog2CopyRatio);
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            return result;
        }
    }

    public CopyRatioSegmentationResult(final List<CopyRatioSegment> segments) {
        Utils.nonEmpty(segments);
        this.segments = Collections.unmodifiableList(segments);
    }

    public void write(final File file,
                      final String sampleName) {
        Utils.nonNull(file);
        Utils.nonNull(sampleName);
        try (final TableWriter<CopyRatioSegment> writer =
                     TableUtils.writer(file, COPY_RATIO_SEGMENT_FILE_TABLE_COLUMNS,
                             (segment, dataLine) ->
                                     dataLine.append(sampleName)
                                             .append(segment.interval.getContig())
                                             .append(segment.interval.getStart())
                                             .append(segment.interval.getEnd())
                                             .append(segment.numDataPoints)
                                             .append(Math.pow(2., segment.meanLog2CopyRatio)))) {
            writer.writeAllRecords(segments);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(file, e);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final CopyRatioSegmentationResult that = (CopyRatioSegmentationResult) o;
        return segments.equals(that.segments);
    }

    @Override
    public int hashCode() {
        return segments.hashCode();
    }
}