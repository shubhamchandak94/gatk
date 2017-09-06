package org.broadinstitute.hellbender.utils.smithwaterman;

import java.io.Closeable;

public interface SmithWatermanAligner extends Closeable{
    SmithWatermanAlignment align(final byte[] ref, final byte[] alt);

    @Override
    default void close() {}
}
