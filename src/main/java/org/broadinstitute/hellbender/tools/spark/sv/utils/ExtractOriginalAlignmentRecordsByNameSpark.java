package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashSet;

/**
 * Created by shuang on 3/3/17.
 */
@CommandLineProgramProperties(summary="Find original alignment records that have the requested read names and outputs them.",
        oneLineSummary="Dump original alignment records that have the requested read names.",
        usageExample = "ExtractOriginalAlignmentRecordsByNameSpark \\" +
                "-I /path/to/my/dir/allAlignments.sam \\" +
                "-O /path/to/my/dir/extractedOriginalAlignments.sam \\" +
                "--readNameFile /path/to/my/dir/readNames.txt",
        omitFromCommandLine = true,
        programGroup = StructuralVariationSparkProgramGroup.class)
@BetaFeature
public final class ExtractOriginalAlignmentRecordsByNameSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(ExtractOriginalAlignmentRecordsByNameSpark.class);

    @Argument(doc = "file containing list of read names", fullName = "readNameFile")
    private String readNameFile;

    @Argument(doc = "file to write reads to", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputSAM;

    @Argument(doc = "to require RG tag on reads or not [false]", shortName = "rg",
            fullName = "requireRG", optional = true)
    private boolean require = false;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    protected void runTool( final JavaSparkContext ctx ) {

        final Broadcast<HashSet<String>> namesToLookForBroadcast = ctx.broadcast(parseReadNames());

        final JavaRDD<GATKRead> reads =
                getUnfilteredReads().repartition(80)
                        .filter(read -> namesToLookForBroadcast.getValue().contains(read.getName())).cache();
        localLogger.info("Found these many alignments: " + reads.count());

        writeReads(ctx, outputSAM, reads);
    }

    private HashSet<String> parseReadNames() {

        try ( final BufferedReader rdr =
                      new BufferedReader(new InputStreamReader(BucketUtils.openFile(readNameFile))) ) {
            final HashSet<String> namesToLookFor = new HashSet<>();
            String line;
            while ( (line = rdr.readLine()) != null ) {
                namesToLookFor.add(line.replace("@", "")
                                       .replace("/1", "")
                                       .replace("/2", ""));
            }
            localLogger.info("Number of read names: " + namesToLookFor.size());
            return namesToLookFor;
        }
        catch ( final IOException ioe ) {
            throw new GATKException("Unable to read names file from " + readNameFile, ioe);
        }
    }
}
