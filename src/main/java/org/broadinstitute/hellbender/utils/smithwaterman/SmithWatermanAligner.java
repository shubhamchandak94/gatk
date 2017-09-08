package org.broadinstitute.hellbender.utils.smithwaterman;

import org.broadinstitute.gatk.nativebindings.smithwaterman.SWAlignerArguments;

import java.io.Closeable;
import java.util.function.BiFunction;

public interface SmithWatermanAligner extends Closeable{

    // match=1, mismatch = -1/3, gap=-(1+k/3)
    SWAlignerArguments.Parameters ORIGINAL_DEFAULT_PARAMETERS = new SWAlignerArguments.Parameters(3, -1, -4, -3);
    SWAlignerArguments.Parameters STANDARD_NGS_PARAMETERS = new SWAlignerArguments.Parameters(25, -50, -110, -6);

    SmithWatermanAlignment align(final byte[] ref, final byte[] alt);

    @Override
    default void close() {}

    enum Implementation {
        FASTEST_AVAILABLE(SWPairwiseAlignment::new),
        JAVA(SWPairwiseAlignment::new);

        private final BiFunction<SWAlignerArguments.Parameters, SWAlignerArguments.OverhangStrategy, SmithWatermanAligner> createAligner;

        Implementation(final BiFunction<SWAlignerArguments.Parameters, SWAlignerArguments.OverhangStrategy, SmithWatermanAligner> createAligner ){
                this.createAligner = createAligner;
        }

        private SmithWatermanAligner createAligner(final SWAlignerArguments.Parameters parameters, final SWAlignerArguments.OverhangStrategy strategy){
            return createAligner.apply(parameters, strategy);
        }
    }

    static SmithWatermanAligner getAligner(final SWAlignerArguments.Parameters parameters, final SWAlignerArguments.OverhangStrategy strategy, final Implementation type) {
        return type.createAligner(parameters, strategy);
    }
}
