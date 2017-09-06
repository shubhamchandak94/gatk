package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.*;
import java.util.Collection;
import java.util.Set;

public class SVKmerUtils {

    /**
     * Read a file of kmers.
     * Each line must be exactly
     * {@link org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection#KMER_SIZE}
     * characters long, and must match [ACGT]*.
     */
    public static Set<SVKmer> readKmersFile(final int kSize, final String kmersFile,
                                            final SVKmer kmer ) {
        final Set<SVKmer> kmers;

        try ( final BufferedReader rdr =
                      new BufferedReader(new InputStreamReader(BucketUtils.openFile(kmersFile))) ) {
            final long fileLength = BucketUtils.fileSize(kmersFile);
            kmers = new HopscotchSet<>((int)(fileLength/(kSize+1)));
            String line;
            while ( (line = rdr.readLine()) != null ) {
                if ( line.length() != kSize ) {
                    throw new GATKException("SVKmer kill set contains a line of length " + line.length() +
                            " but we were expecting K=" + kSize);
                }

                final SVKmerizer kmerizer = new SVKmerizer(line, kSize, 1, new SVKmerLong(kSize));
                if ( !kmerizer.hasNext() ) {
                    throw new GATKException("Unable to kmerize the kmer kill set string '" + line + "'.");
                }

                kmers.add(kmerizer.next());
            }
        }
        catch ( final IOException ioe ) {
            throw new GATKException("Unable to read kmers from "+kmersFile, ioe);
        }

        return kmers;
    }

    /** Write kmers to file. */
    public static <KType extends SVKmer> void writeKmersFile(final int kSize, final String kmersFile,
                                                             final Collection<KType> kmers ) {
        try ( final Writer writer =
                      new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(kmersFile))) ) {
            for ( final KType kmer : kmers ) {
                writer.write(kmer.toString(kSize));
                writer.write('\n');
            }
        }
        catch ( final IOException ioe ) {
            throw new GATKException("Unable to write kmers to "+kmersFile, ioe);
        }
    }
}
