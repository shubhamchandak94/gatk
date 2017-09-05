package org.broadinstitute.hellbender.utils.smithwaterman;

public interface SmithWatermanAligner {
    SmithWatermanAlignment align(final byte[] ref, final byte[] alt);

}
