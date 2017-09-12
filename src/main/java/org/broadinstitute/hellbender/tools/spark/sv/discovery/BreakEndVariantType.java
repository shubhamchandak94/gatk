package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.io.IOException;
import java.util.Collections;
import java.util.Map;

public abstract class BreakEndVariantType extends SvType {

    protected final boolean isTheUpstreamMate;

    public final boolean isTheUpstreamMate() {
        return isTheUpstreamMate;
    }

    protected static SimpleInterval getMateRefLoc(final NovelAdjacencyReferenceLocations narl, final boolean forUpstreamLoc) {
        return forUpstreamLoc ? narl.leftJustifiedRightRefLoc : narl.leftJustifiedLeftRefLoc;
    }

    @Override
    public final String toString() {
        return GATKSVVCFConstants.BREAKEND_STR;
    }

    // see VCF spec 4.2 or 4.3 for BND format ALT allele field for SV
    BreakEndVariantType(final String id, final Map<String, String> typeSpecificExtraAttributes,
                        final byte[] bases, final boolean bracketPointsLeft, final SimpleInterval novelAdjRefLoc,
                        final boolean basesFirst, final boolean isTheUpstreamMate) {
        super(id, constructAltAllele(bases, bracketPointsLeft, novelAdjRefLoc, basesFirst), INAPPLICABLE_LENGTH, typeSpecificExtraAttributes);
        this.isTheUpstreamMate = isTheUpstreamMate;
    }

    private static Allele constructAltAllele(final byte[] bases, final boolean bracketPointsLeft, final SimpleInterval novelAdjRefLoc,
                                             final boolean basesFirst) {
        final String s = bracketPointsLeft ? "]" + novelAdjRefLoc.getContig() + ":" + novelAdjRefLoc.getStart() + "]"
                                           : "[" + novelAdjRefLoc.getContig() + ":" + novelAdjRefLoc.getStart() + "[";
        return Allele.create( basesFirst ? new String(bases) + s  : s + new String(bases) );
    }

    static abstract class StrandSwitchBND extends BreakEndVariantType {

        StrandSwitchBND(final String id, final Map<String, String> typeSpecificExtraAttributes,
                        final byte[] bases, final boolean bracketPointsLeft, final SimpleInterval novelAdjRefLoc,
                        final boolean basesFirst, boolean isTheUpstreamMate) {
            super(id, typeSpecificExtraAttributes, bases, bracketPointsLeft, novelAdjRefLoc, basesFirst, isTheUpstreamMate);
        }


        static byte[] extractBasesForAltAllele(final NovelAdjacencyReferenceLocations narl, final boolean forUpstreamLoc,
                                               final ReferenceMultiSource reference) {
            try {
                final byte[] ref = reference
                        .getReferenceBases(null, forUpstreamLoc ? narl.leftJustifiedLeftRefLoc :
                                                                                narl.leftJustifiedRightRefLoc)
                        .getBases();
                final String ins = narl.complication.getInsertedSequenceForwardStrandRep();
                if (ins.isEmpty()) {
                    return ref;
                } else {
                    return forUpstreamLoc ? ArrayUtils.addAll(ref, ins.getBytes())
                            : ArrayUtils.addAll(ref, SequenceUtil.reverseComplement(ins).getBytes());
                }
            } catch (final IOException ioex) {
                throw new GATKException("Could not read reference for extracting reference bases.", ioex);
            }
        }

        static SimpleInterval getTheOtherRefLoc(final NovelAdjacencyReferenceLocations narl, final boolean forUpstreamLoc) {
            return forUpstreamLoc ? narl.leftJustifiedRightRefLoc : narl.leftJustifiedLeftRefLoc;
        }
    }

    public static final class INV55BND extends StrandSwitchBND {

        /**
         * Technically, a strand switch breakpoint should have two VCF records, hence we also have {@code forUpstreamLoc}.
         */
        private INV55BND(final NovelAdjacencyReferenceLocations narl, final boolean forUpstreamLoc,
                         final ReferenceMultiSource reference) {
            super(getIDString(narl, forUpstreamLoc),
                    Collections.singletonMap(GATKSVVCFConstants.INV55, ""),
                    extractBasesForAltAllele(narl, forUpstreamLoc, reference), true,
                    getMateRefLoc(narl, forUpstreamLoc), true, forUpstreamLoc);
        }

        public static Tuple2<BreakEndVariantType, BreakEndVariantType> getOrderedMates(final NovelAdjacencyReferenceLocations narl,
                                                                                       final ReferenceMultiSource reference) {
            final BreakEndVariantType bkpt_1 = new BreakEndVariantType.INV55BND(narl, true, reference);
            final BreakEndVariantType bkpt_2 = new BreakEndVariantType.INV55BND(narl, false, reference);
            return new Tuple2<>(bkpt_1, bkpt_2);
        }

        private static String getIDString(final NovelAdjacencyReferenceLocations narl, final boolean forUpstreamLoc) {

            return GATKSVVCFConstants.BREAKEND_STR + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    GATKSVVCFConstants.INV55 + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    narl.leftJustifiedLeftRefLoc.getContig() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    narl.leftJustifiedLeftRefLoc.getEnd() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    narl.leftJustifiedRightRefLoc.getStart() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    (forUpstreamLoc ? "1" : "2") ;
        }
    }

    public static final class INV33BND extends StrandSwitchBND {

        /**
         * Technically, a strand switch breakpoint should have two VCF records, hence we also have {@code forUpstreamLoc}.
         */
        private INV33BND(final NovelAdjacencyReferenceLocations narl, final boolean forUpstreamLoc,
                         final ReferenceMultiSource reference) {
            super(getIDString(narl, forUpstreamLoc),
                    Collections.singletonMap(GATKSVVCFConstants.INV33, ""),
                    extractBasesForAltAllele(narl, forUpstreamLoc, reference), false,
                    getMateRefLoc(narl, forUpstreamLoc), false, forUpstreamLoc);
        }

        public static Tuple2<BreakEndVariantType, BreakEndVariantType> getOrderedMates(final NovelAdjacencyReferenceLocations narl,
                                                                                       final ReferenceMultiSource reference) {
            final BreakEndVariantType bkpt_1 = new BreakEndVariantType.INV33BND(narl, true, reference);
            final BreakEndVariantType bkpt_2 = new BreakEndVariantType.INV33BND(narl, false, reference);
            return new Tuple2<>(bkpt_1, bkpt_2);
        }

        private static String getIDString(final NovelAdjacencyReferenceLocations narl, final boolean forUpstreamLoc) {

            return GATKSVVCFConstants.BREAKEND_STR + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    GATKSVVCFConstants.INV33 + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    narl.leftJustifiedLeftRefLoc.getContig() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    narl.leftJustifiedLeftRefLoc.getEnd() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    narl.leftJustifiedRightRefLoc.getStart() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    (forUpstreamLoc ? "1" : "2") ;
        }
    }

    /**
     * Generic variant type for suspected "translocations", including what could be a breakpoint for a duplication or
     * an inter-chromosome novel adjacency.
     */
    public static final class TransLocBND extends BreakEndVariantType {
        private static Map<String, String> emptyMap = Collections.emptyMap();

        private TransLocBND(final NovelAdjacencyReferenceLocations narl, final boolean forUpstreamLoc,
                            final ReferenceMultiSource reference, final SAMSequenceDictionary referenceDictionary,
                            final boolean basesFirst, final boolean bracketPointsLeft) {
            super(getIDString(narl, forUpstreamLoc), emptyMap,
                    extractBasesForAltAllele(narl, forUpstreamLoc, reference, referenceDictionary), bracketPointsLeft,
                    getMateRefLoc(narl, forUpstreamLoc), basesFirst, forUpstreamLoc);
        }

        public static Tuple2<BreakEndVariantType, BreakEndVariantType> getOrderedMates(final NovelAdjacencyReferenceLocations narl,
                                                                                       final ReferenceMultiSource reference,
                                                                                       final SAMSequenceDictionary referenceDictionary) {
            final boolean isSameChr = narl.leftJustifiedLeftRefLoc.getContig().equals(narl.leftJustifiedRightRefLoc.getContig());
            final BreakEndVariantType bkpt_1, bkpt_2;
            if (isSameChr) {
                bkpt_1 = new BreakEndVariantType.TransLocBND(narl, true, reference, referenceDictionary, false, true);
                bkpt_2 = new BreakEndVariantType.TransLocBND(narl, false, reference, referenceDictionary, true, false);
            } else {
                if (narl.strandSwitch == StrandSwitch.NO_SWITCH) {
                    final boolean isFirstOfPartner = IntervalUtils.compareContigs(narl.leftJustifiedLeftRefLoc,
                                                                                  narl.leftJustifiedRightRefLoc, referenceDictionary) < 0;
                    bkpt_1 = new BreakEndVariantType.TransLocBND(narl, true, reference, referenceDictionary,
                            isFirstOfPartner, !isFirstOfPartner);
                    bkpt_2 = new BreakEndVariantType.TransLocBND(narl, false, reference, referenceDictionary,
                            !isFirstOfPartner, isFirstOfPartner);
                } else {
                    bkpt_1 = new BreakEndVariantType.TransLocBND(narl, true, reference, referenceDictionary,
                            narl.strandSwitch == StrandSwitch.FORWARD_TO_REVERSE, narl.strandSwitch == StrandSwitch.FORWARD_TO_REVERSE);
                    bkpt_2 = new BreakEndVariantType.TransLocBND(narl, false, reference, referenceDictionary,
                            narl.strandSwitch != StrandSwitch.FORWARD_TO_REVERSE, narl.strandSwitch != StrandSwitch.FORWARD_TO_REVERSE);
                }
            }
            return new Tuple2<>(bkpt_1, bkpt_2);
        }

        private static String getIDString(final NovelAdjacencyReferenceLocations narl, final boolean forUpstreamLoc) {
            return GATKSVVCFConstants.BREAKEND_STR + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    narl.leftJustifiedLeftRefLoc.getContig() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    narl.leftJustifiedLeftRefLoc.getEnd() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    narl.leftJustifiedRightRefLoc.getContig() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    narl.leftJustifiedRightRefLoc.getStart() + GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    (forUpstreamLoc ? "1" : "2");
        }

        static byte[] extractBasesForAltAllele(final NovelAdjacencyReferenceLocations narl, final boolean forUpstreamLoc,
                                               final ReferenceMultiSource reference, final SAMSequenceDictionary referenceDictionary) {
            try {
                final SimpleInterval refLoc = forUpstreamLoc ? narl.leftJustifiedLeftRefLoc : narl.leftJustifiedRightRefLoc;
                final byte[] ref = reference.getReferenceBases(null, refLoc).getBases();
                final String ins = narl.complication.getInsertedSequenceForwardStrandRep();
                if (ins.isEmpty()) {
                    return ref;
                } else {
                    if (narl.leftJustifiedLeftRefLoc.getContig().equals(narl.leftJustifiedRightRefLoc.getContig())
                            || (narl.strandSwitch == StrandSwitch.NO_SWITCH && IntervalUtils.compareContigs(narl.leftJustifiedLeftRefLoc, narl.leftJustifiedRightRefLoc, referenceDictionary) > 0)
                            || narl.strandSwitch == StrandSwitch.REVERSE_TO_FORWARD) {
                        return forUpstreamLoc ? ArrayUtils.addAll(ins.getBytes(), ref) : ArrayUtils.addAll(ref, ins.getBytes());
                    } else {
                        return forUpstreamLoc ? ArrayUtils.addAll(ref, ins.getBytes()) : ArrayUtils.addAll(ins.getBytes(), ref);
                    }
                }
            } catch (final IOException ioex) {
                throw new GATKException("Could not read reference for extracting reference bases.", ioex);
            }
        }
    }
}
