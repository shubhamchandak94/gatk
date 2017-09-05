package org.broadinstitute.hellbender.utils.smithwaterman;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWAlignerArguments;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWAlignerNativeBinding;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWAlignmentResult;

public abstract class SWNativeAlignerWrapper implements SmithWatermanAligner {
    private final SWAlignerNativeBinding aligner;

    public SWNativeAlignerWrapper(final SWAlignerNativeBinding aligner, final SWAlignerArguments.Weights weights, final SWAlignerArguments.OverhangStrategy overhangStrategy) {
        this.aligner = aligner;
        this.aligner.initialize(new SWAlignerArguments(overhangStrategy, weights));
    }

    @Override
    public SmithWatermanAlignment align(final byte[] ref, final byte[] alt){
        final SWAlignmentResult alignment = aligner.align(ref, alt);
        return new SWNativeResultWrapper(alignment);
    }

    private static final class SWNativeResultWrapper implements SmithWatermanAlignment {
        private final Cigar cigar;
        private final int alignmentOffset;

        public SWNativeResultWrapper(final SWAlignmentResult nativeResult) {
            this.cigar = TextCigarCodec.decode(nativeResult.cigar);
            this.alignmentOffset = nativeResult.alignment_offset;
        }

        @Override
        public Cigar getCigar() {
            return cigar;
        }

        @Override
        public int getAlignmentOffset() {
            return alignmentOffset;
        }
    }
}
