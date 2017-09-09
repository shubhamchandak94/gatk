package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.*;
import scala.Tuple2;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

final class SuspectedTransLocDetector implements VariantDetectorFromLocalAssemblyContigAlignments {

    @SuppressWarnings("unchecked")
    private static final List<String> EMPTY_INSERTION_MAPPINGS = Collections.EMPTY_LIST;

    @Override
    public void inferSvAndWriteVCF(final JavaRDD<AlignedContig> localAssemblyContigs, final String vcfOutputFileName,
                                   final Broadcast<ReferenceMultiSource> broadcastReference,
                                   final Broadcast<SAMSequenceDictionary> sequenceDictionaryBroadcast,
                                   final Logger toolLogger) {
        localAssemblyContigs.cache();
        toolLogger.info(localAssemblyContigs.count() + " chimeras indicating strand-switch-less breakpoints");

        final JavaPairRDD<ChimericAlignment, byte[]> chimeraAndSequence =
                localAssemblyContigs
                        .filter(tig ->
                                SimpleStrandSwitchVariantDetector.splitPairStrongEnoughEvidenceForCA(tig.alignmentIntervals.get(0),
                                        tig.alignmentIntervals.get(1), SimpleStrandSwitchVariantDetector.MORE_RELAXED_ALIGNMENT_MIN_MQ,
                                        0))
                        .mapToPair(SuspectedTransLocDetector::convertAlignmentIntervalsToChimericAlignment).cache();

        final JavaRDD<VariantContext> annotatedBNDs =
                chimeraAndSequence
                        .mapToPair(pair -> new Tuple2<>(new NovelAdjacencyReferenceLocations(pair._1, pair._2), pair._1))
                        .groupByKey()
                        .mapToPair(noveltyAndEvidence -> inferBNDType(noveltyAndEvidence, broadcastReference.getValue()))
                        .flatMap(noveltyTypeAndEvidence ->
                                AnnotatedVariantProducer
                                        .produceAnnotatedBNDmatesVcFromNovelAdjacency(noveltyTypeAndEvidence._1,
                                                noveltyTypeAndEvidence._2._1, noveltyTypeAndEvidence._2._2, broadcastReference).iterator());

//        SVVCFWriter.writeVCF(vcfOutputFileName.replace(".vcf", "_transBND.vcf"),
//                sequenceDictionaryBroadcast.getValue(), annotatedBNDs, toolLogger);
    }

    private static Tuple2<ChimericAlignment, byte[]> convertAlignmentIntervalsToChimericAlignment(final AlignedContig contig) {
        // TODO: 9/9/17 this default empty insertion mapping treatment is temporary and should be fixed later
        return new Tuple2<>(new ChimericAlignment(contig.alignmentIntervals.get(0), contig.alignmentIntervals.get(1),
                EMPTY_INSERTION_MAPPINGS, contig.contigName), contig.contigSequence);
    }

    private static Tuple2<NovelAdjacencyReferenceLocations, Tuple2<List<SvType>, Iterable<ChimericAlignment>>>
    inferBNDType(final Tuple2<NovelAdjacencyReferenceLocations, Iterable<ChimericAlignment>> noveltyAndEvidence,
                 final ReferenceMultiSource reference) {

        final NovelAdjacencyReferenceLocations novelAdjacency = noveltyAndEvidence._1;
        final Iterable<ChimericAlignment> chimericAlignments = noveltyAndEvidence._2;
        final BreakEndVariantType bkpt_1 = null, bkpt_2 = null;
//        if (novelAdjacency.strandSwitch == StrandSwitch.FORWARD_TO_REVERSE) {
//            bkpt_1 = new BreakEndVariantType.INV55BND(novelAdjacency, true, reference);
//            bkpt_2 = new BreakEndVariantType.INV55BND(novelAdjacency, false, reference);
//        } else if (novelAdjacency.strandSwitch == StrandSwitch.REVERSE_TO_FORWARD) {
//            bkpt_1 = new BreakEndVariantType.INV33BND(novelAdjacency, true, reference);
//            bkpt_2 = new BreakEndVariantType.INV33BND(novelAdjacency, false, reference);
//        } else {
//            throw new GATKException("Wrong type of novel adjacency sent to wrong analysis pathway: no strand-switch being sent to strand-switch path. \n" +
//                    Utils.stream(chimericAlignments).map(ChimericAlignment::onErrStringRep).collect(Collectors.toList()));
//        }

        return new Tuple2<>(novelAdjacency, new Tuple2<>(Arrays.asList(bkpt_1, bkpt_2), chimericAlignments));
    }
}
