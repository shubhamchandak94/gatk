package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.*;
import java.util.*;

public final class SVFileUtils {

    private static final String REFERENCE_GAP_INTERVAL_FILE_COMMENT_LINE_PROMPT = "#";

    /**
     * Creates a directory, in local FS, HDFS, or Google buckets to write to.
     */
    public static boolean createDirToWriteTo(final String pathString) {
        try {
            Utils.nonNull(pathString);
            if ( java.nio.file.Files.exists(java.nio.file.Paths.get(pathString)) )
                throw new IOException("Directory to be created already exists: " + pathString);

            final boolean isSuccessful;
            if (BucketUtils.isHadoopUrl(pathString)) {
                isSuccessful = org.apache.hadoop.fs.FileSystem.get(new org.apache.hadoop.conf.Configuration()).mkdirs(new org.apache.hadoop.fs.Path(pathString));
            } else {
                final java.nio.file.Path dir = java.nio.file.Files.createDirectory(IOUtils.getPath(pathString));
                isSuccessful = java.nio.file.Files.isDirectory(dir) && java.nio.file.Files.isWritable(dir);
            }
            return isSuccessful;
        } catch (final IOException x) {
            throw new GATKException("Could not create directory: " + pathString, x);
        }
    }

    public static void writeLinesToSingleFile(final Iterator<String> linesToWrite, final String fileName) {
        try ( final OutputStream writer =
                      new BufferedOutputStream(BucketUtils.createFile(fileName)) ) {
            while (linesToWrite.hasNext()) {
                writer.write(linesToWrite.next().getBytes()); writer.write('\n');
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't write "+fileName, ioe);
        }
    }

    public static void writeSAMFile(final Iterator<SAMRecord> alignments, final SAMFileHeader header, final String outputName,
                                    final boolean preOrdered) {
        try (final SAMFileWriter writer = ReadUtils.createCommonSAMWriter(new File(outputName), null, header,
                preOrdered, true, false) ) {
            alignments.forEachRemaining(writer::addAlignment);
        } catch ( final UncheckedIOException ie) {
            throw new GATKException("Can't write SAM file to the specified location: " + outputName, ie);
        }
    }

    /** Read intervals from file. */
    public static List<SVInterval> readIntervalsFile(final String intervalsFile,
                                                     final Map<String, Integer> contigNameMap ) {
        final List<SVInterval> intervals;
        try ( final BufferedReader rdr =
                      new BufferedReader(new InputStreamReader(BucketUtils.openFile(intervalsFile))) ) {
            final int INTERVAL_FILE_LINE_LENGTH_GUESS = 25;
            final long sizeGuess = BucketUtils.fileSize(intervalsFile)/INTERVAL_FILE_LINE_LENGTH_GUESS;
            intervals = new ArrayList<>((int)sizeGuess);
            String line;
            int lineNo = 0;
            while ( (line = rdr.readLine()) != null ) {
                ++lineNo;
                if (line.startsWith(REFERENCE_GAP_INTERVAL_FILE_COMMENT_LINE_PROMPT)) {
                    continue;
                }
                final String[] tokens = line.split("\t");
                if ( tokens.length != 3 ) {
                    throw new GATKException("Interval file "+intervalsFile+" line "+
                                            lineNo+" did not contain 3 columns: "+line);
                }
                try {
                    final Integer contigId = contigNameMap.get(tokens[0]);
                    if ( contigId == null ) throw new GATKException("contig name "+tokens[0]+" not in dictionary");
                    final int start = Integer.valueOf(tokens[1]);
                    final int end = Integer.valueOf(tokens[2]);
                    intervals.add(new SVInterval(contigId, start, end));
                }
                catch ( final Exception e ) {
                    throw new GATKException("Unable to parse interval file "+intervalsFile+" line "+lineNo+": "+line, e);
                }
            }
        }
        catch ( final IOException ioe ) {
            throw new GATKException("Unable to read intervals from "+intervalsFile, ioe);
        }
        return intervals;
    }

    /** Write intervals to a file. */
    public static void writeIntervalsFile(final String intervalsFile,
                                          final Collection<SVInterval> intervals, final List<String> contigNames ) {

        writeLinesToSingleFile(intervals.stream().map(interval ->
                contigNames.get(interval.getContig()) + "\t" + interval.getStart() + "\t" + interval.getEnd() + "\n").iterator(),
                intervalsFile);
    }
}
