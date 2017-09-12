package org.broadinstitute.hellbender.cmdline;

import org.broadinstitute.hellbender.exceptions.GATKException;

/**
 * The GATK and Picard CommandLineProgram classes have no common ancestor. This class acts as a shim for use
 * with Picard tools from within GATK.
 *
 * NOTE that this class does not have it's own CommandLineProgramProperties annotation.
 */
public class PicardCommandLineProgramExecutor extends CommandLineProgram {

    // Our wrapped Picard command line program, to which we forward subsequent calls.
    final private picard.cmdline.CommandLineProgram picardCommandLineProgram;

    public PicardCommandLineProgramExecutor(final picard.cmdline.CommandLineProgram picardCommandLineProgram) {
        this.picardCommandLineProgram = picardCommandLineProgram;
    }

    /**
     * This is the entry point for Picard tools that are called from GATK.
     */
    public Object instanceMain(final String[] argv) {
        return picardCommandLineProgram.instanceMain(argv);
    }

    @Override
    public Object doWork() {
        // This method should never be called directly. Call instanceMain instead.
        throw new GATKException.ShouldNeverReachHereException(
                String.format("Attempt to call the doWork method on the Picard tool \"%s\" directly.",
                        picardCommandLineProgram.getClass().getName()));
    }
}
