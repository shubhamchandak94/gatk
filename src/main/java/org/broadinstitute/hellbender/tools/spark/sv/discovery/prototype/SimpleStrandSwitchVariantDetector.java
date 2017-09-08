package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.RDDUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;


final class SimpleStrandSwitchVariantDetector implements VariantDetectorFromLocalAssemblyContigAlignments {

    @SuppressWarnings("unchecked")
    private static final List<String> EMPTY_INSERTION_MAPPINGS = Collections.EMPTY_LIST;

    static final int MORE_RELAXED_ALIGNMENT_MIN_LENGTH = 30;
    static final int MORE_RELAXED_ALIGNMENT_MIN_MQ = 20;

    @Override
    public void inferSvAndWriteVCF(final JavaRDD<AlignedContig> contigs, final String vcfOutputFileName,
                                   final Broadcast<ReferenceMultiSource> broadcastReference,
                                   final Broadcast<SAMSequenceDictionary> sequenceDictionaryBroadcast,
                                   final Logger toolLogger) {

        contigs.cache();
        toolLogger.info(contigs.count() + " chimeras indicating either 1) simple strand-switch breakpoints, or 2) inverted duplication.");

        // split between suspected inv dup VS strand-switch breakpoints
        // logic flow: first modify alignments (heuristically) if there are overlaps on read between the alignments,
        //             then split the input reads into two classes--those judged by IsLikelyInvertedDuplication are likely invdup and those aren't
        //             finally send the two split reads down different path, one for invdup and one for BND records
        final Tuple2<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>> invDupAndStrandSwitchBreakpoints =
                RDDUtils.split(contigs.map(ContigAlignmentsModifier::removeOverlap),
                        contig -> BreakpointComplications.isLikelyInvertedDuplication(contig.alignmentIntervals.get(0),
                                contig.alignmentIntervals.get(1)), true);

        final JavaRDD<VariantContext> simpleStrandSwitchBkpts =
                dealWithSimpleStrandSwitchBkpts(invDupAndStrandSwitchBreakpoints._2, broadcastReference, toolLogger);
        SVVCFWriter.writeVCF(vcfOutputFileName.replace(".vcf", "_simpleSS.vcf"),
                sequenceDictionaryBroadcast.getValue(), simpleStrandSwitchBkpts, toolLogger);

        // TODO: 8/23/17 add inv dup code in the next pr
    }
    // =================================================================================================================

    /**
     * @throws IllegalArgumentException if the assumption that the input aligned assembly contig has 2 alignments
     *                                  mapped to the same chr with strand switch is invalid
     */
    static Tuple2<ChimericAlignment, byte[]> convertAlignmentIntervalsToChimericAlignment
    (final AlignedContig contigWith2AIMappedToSameChrAndStrandSwitch) {
        Utils.validateArg(SvDiscoverFromLocalAssemblyContigAlignmentsSpark.isLikelyInvBreakpointOrInsInv(contigWith2AIMappedToSameChrAndStrandSwitch),
                "assumption that input aligned assembly contig has 2 alignments mapped to the same chr with strand switch is invalid.\n" +
                        contigWith2AIMappedToSameChrAndStrandSwitch.toString());

        final AlignmentInterval intervalOne = contigWith2AIMappedToSameChrAndStrandSwitch.alignmentIntervals.get(0),
                                intervalTwo = contigWith2AIMappedToSameChrAndStrandSwitch.alignmentIntervals.get(1);

        // TODO: 8/28/17 this default empty insertion mapping treatment is temporary and should be fixed later
        return new Tuple2<>(new ChimericAlignment(intervalOne, intervalTwo, EMPTY_INSERTION_MAPPINGS,
                contigWith2AIMappedToSameChrAndStrandSwitch.contigName), contigWith2AIMappedToSameChrAndStrandSwitch.contigSequence);
    }

    /**
     * Roughly similar to {@link ChimericAlignment#nextAlignmentMayBeNovelInsertion(AlignmentInterval, AlignmentInterval, Integer)}:
     *  1) either alignment may have very low mapping quality (a more relaxed mapping quality threshold);
     *  2) either alignment may consume only a "short" part of the contig, or if assuming that the alignment consumes
     *     roughly the same amount of ref bases and read bases, has isAlignment that is too short
     */
    static boolean splitPairStrongEnoughEvidenceForCA(final AlignmentInterval intervalOne,
                                                      final AlignmentInterval intervalTwo,
                                                      final int mapQThresholdInclusive,
                                                      final int alignmentLengthThresholdInclusive) {

        if (intervalOne.mapQual < mapQThresholdInclusive || intervalTwo.mapQual < mapQThresholdInclusive)
            return false;

        final int overlap = AlignmentInterval.overlapOnContig(intervalOne, intervalTwo);

        final int x = intervalOne.endInAssembledContig - intervalOne.startInAssembledContig + 1,
                  y = intervalTwo.endInAssembledContig - intervalTwo.startInAssembledContig + 1;

        return Math.min(x - overlap, y - overlap) >= alignmentLengthThresholdInclusive;
    }

    // workflow manager for simple strand-switch alignment contigs
    private JavaRDD<VariantContext> dealWithSimpleStrandSwitchBkpts(final JavaRDD<AlignedContig> contigs,
                                                                    final Broadcast<ReferenceMultiSource> broadcastReference,
                                                                    final Logger toolLogger) {

        final JavaPairRDD<ChimericAlignment, byte[]> simpleStrandSwitchBkpts =
                contigs
                        .filter(tig ->
                                splitPairStrongEnoughEvidenceForCA(tig.alignmentIntervals.get(0), tig.alignmentIntervals.get(1),
                                        MORE_RELAXED_ALIGNMENT_MIN_MQ,  MORE_RELAXED_ALIGNMENT_MIN_LENGTH))
                        .mapToPair(SimpleStrandSwitchVariantDetector::convertAlignmentIntervalsToChimericAlignment).cache();

        toolLogger.info(simpleStrandSwitchBkpts.count() + " chimeras indicating simple strand-switch breakpoints.");

        return simpleStrandSwitchBkpts
                .mapToPair(pair -> new Tuple2<>(new NovelAdjacencyReferenceLocations(pair._1, pair._2), pair._1))
                .groupByKey()
                .mapToPair(noveltyAndEvidence -> inferBNDType(noveltyAndEvidence, broadcastReference.getValue()))
                .flatMap(noveltyTypeAndEvidence ->
                        AnnotatedVariantProducer
                                .produceAnnotatedBNDmatesVcFromNovelAdjacency(noveltyTypeAndEvidence._1,
                                        noveltyTypeAndEvidence._2._1, noveltyTypeAndEvidence._2._2, broadcastReference).iterator());
    }

    private static Tuple2<NovelAdjacencyReferenceLocations, Tuple2<List<SvType>, Iterable<ChimericAlignment>>>
    inferBNDType(final Tuple2<NovelAdjacencyReferenceLocations, Iterable<ChimericAlignment>> noveltyAndEvidence,
                 final ReferenceMultiSource reference) {

        final NovelAdjacencyReferenceLocations novelAdjacency = noveltyAndEvidence._1;
        final Iterable<ChimericAlignment> chimericAlignments = noveltyAndEvidence._2;
        final BreakEndVariantType bkpt_1, bkpt_2;
        if (novelAdjacency.strandSwitch == StrandSwitch.FORWARD_TO_REVERSE) {
            bkpt_1 = new BreakEndVariantType.INV55BND(novelAdjacency, true, reference);
            bkpt_2 = new BreakEndVariantType.INV55BND(novelAdjacency, false, reference);
        } else if (novelAdjacency.strandSwitch == StrandSwitch.REVERSE_TO_FORWARD) {
            bkpt_1 = new BreakEndVariantType.INV33BND(novelAdjacency, true, reference);
            bkpt_2 = new BreakEndVariantType.INV33BND(novelAdjacency, false, reference);
        } else {
            throw new GATKException("Wrong type of novel adjacency sent to wrong analysis pathway: no strand-switch being sent to strand-switch path. \n" +
                    Utils.stream(chimericAlignments).map(ChimericAlignment::onErrStringRep).collect(Collectors.toList()));
        }

        return new Tuple2<>(novelAdjacency, new Tuple2<>(Arrays.asList(bkpt_1, bkpt_2), chimericAlignments));
    }

    // TODO: 8/23/17 add inv dup code in the next PR
}
