package org.broadinstitute.hellbender.utils.smithwaterman;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWAlignerArguments;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Pairwise discrete smith-waterman alignment
 *
 * ************************************************************************
 * ****                    IMPORTANT NOTE:                             ****
 * ****  This class assumes that all bytes come from UPPERCASED chars! ****
 * ************************************************************************
 */
public final class SWPairwiseAlignment implements SmithWatermanAligner {

    /**
     * The state of a trace step through the matrix
     */
    protected enum State {
        MATCH,
        INSERTION,
        DELETION,
        CLIP
    }

    private final SWAlignerArguments.Parameters parameters;
    private final SWAlignerArguments.OverhangStrategy overhangStrategy;

    private static final boolean cutoff = false;


    /**
     * The SW scoring matrix, stored for debugging purposes if keepScoringMatrix is true
     */
    private int[][] SW = null;

    /**
     * Only for testing purposes in the SWPairwiseAlignmentMain function
     * set to true to keep SW scoring matrix after align call
     */
    private static final boolean keepScoringMatrix = false;

    /**
     * Create a new SW pairwise aligner
     *
     * @param parameters the SW parameters to use
     * @param strategy   the overhang strategy to use
     */
    public SWPairwiseAlignment(final SWAlignerArguments.Parameters parameters, final SWAlignerArguments.OverhangStrategy strategy) {
        Utils.nonNull(parameters);
        Utils.nonNull(strategy);
        this.parameters = parameters;
        this.overhangStrategy = strategy;
    }

    /**
     * Aligns the alternate sequence to the reference sequence
     *
     * @param reference  ref sequence
     * @param alternate  alt sequence
     */
    @Override
    public SmithWatermanAlignment align(final byte[] reference, final byte[] alternate) {
        if ( reference == null || reference.length == 0 || alternate == null || alternate.length == 0 ) {
            throw new IllegalArgumentException("Non-null, non-empty sequences are required for the Smith-Waterman calculation");
        }

        // avoid running full Smith-Waterman if there is an exact match of alternate in reference
        int matchIndex = -1;
        if (overhangStrategy == SWAlignerArguments.OverhangStrategy.SOFTCLIP || overhangStrategy == SWAlignerArguments.OverhangStrategy.IGNORE) {
            // Use a substring search to find an exact match of the alternate in the reference
            // NOTE: This approach only works for SOFTCLIP and IGNORE overhang strategies
            matchIndex = Utils.lastIndexOf(reference, alternate);
        }

        if (matchIndex != -1) {
            // generate the alignment result when the substring search was successful
            final List<CigarElement> lce = new ArrayList<>(alternate.length);
            lce.add(makeElement(State.MATCH, alternate.length));
            return  new SWPairwiseAlignmentResult(AlignmentUtils.consolidateCigar(new Cigar(lce)), matchIndex);
        }
        else {
            // run full Smith-Waterman
            final int n = reference.length+1;
            final int m = alternate.length+1;
            final int[][] sw = new int[n][m];
            if ( keepScoringMatrix ) {
                SW = sw;
            }
            final int[][] btrack=new int[n][m];

            calculateMatrix(reference, alternate, sw, btrack, overhangStrategy, parameters);
            return calculateCigar(sw, btrack, overhangStrategy); // length of the segment (continuous matches, insertions or deletions)
        }
    }

    /**
     * Calculates the SW matrices for the given sequences
     *  @param reference  ref sequence
     * @param alternate  alt sequence
     * @param sw         the Smith-Waterman matrix to populate
     * @param btrack     the back track matrix to populate
     * @param overhangStrategy    the strategy to use for dealing with overhangs
     * @param parameters
     */
    private static void calculateMatrix(final byte[] reference, final byte[] alternate, final int[][] sw, final int[][] btrack, final SWAlignerArguments.OverhangStrategy overhangStrategy, SWAlignerArguments.Parameters parameters) {
        if ( reference.length == 0 || alternate.length == 0 ) {
            throw new IllegalArgumentException("Non-null, non-empty sequences are required for the Smith-Waterman calculation");
        }

        final int ncol = sw[0].length;//alternate.length+1; formerly m
        final int nrow = sw.length;// reference.length+1; formerly n

        final int MATRIX_MIN_CUTOFF;   // never let matrix elements drop below this cutoff
        if ( cutoff ) {
            MATRIX_MIN_CUTOFF = 0;
        } else {
            MATRIX_MIN_CUTOFF = (int) -1.0e8;
        }

        final int lowInitValue= Integer.MIN_VALUE/2;
        final int[] best_gap_v = new int[ncol+1];
        Arrays.fill(best_gap_v, lowInitValue);
        final int[] gap_size_v = new int[ncol+1];
        final int[] best_gap_h = new int[nrow+1];
        Arrays.fill(best_gap_h, lowInitValue);
        final int[] gap_size_h = new int[nrow+1];

        // we need to initialize the SW matrix with gap penalties if we want to keep track of indels at the edges of alignments
        if ( overhangStrategy == SWAlignerArguments.OverhangStrategy.INDEL || overhangStrategy == SWAlignerArguments.OverhangStrategy.LEADING_INDEL ) {
            // initialize the first row
            final int[] topRow=sw[0];
            topRow[1]= parameters.getGapOpenPenalty();
            int currentValue = parameters.getGapOpenPenalty();
            for ( int i = 2; i < topRow.length; i++ ) {
                currentValue += parameters.getGapExtendPenalty();
                topRow[i]=currentValue;
            }
            // initialize the first column
            sw[1][0]= parameters.getGapOpenPenalty();
            currentValue = parameters.getGapOpenPenalty();
            for ( int i = 2; i < sw.length; i++ ) {
                currentValue += parameters.getGapExtendPenalty();
                sw[i][0]=currentValue;
            }
        }
        // build smith-waterman matrix and keep backtrack info:
        int[] curRow=sw[0];

        //access is pricey if done enough times so we extract those out
        final int w_open = parameters.getGapOpenPenalty();
        final int w_extend = parameters.getGapExtendPenalty();
        final int w_match = parameters.getMatchValue();
        final int w_mismatch = parameters.getMismatchPenalty();

        //array length checks are expensive in tight loops so extract the length out
        for ( int i = 1, sw_length = sw.length; i < sw_length ; i++ ) {
            final byte a_base = reference[i-1]; // letter in a at the current pos
            final int[] lastRow=curRow;
            curRow=sw[i];
            final int[] curBackTrackRow=btrack[i];

            //array length checks are expensive in tight loops so extract the length out
            for ( int j = 1, curRow_length = curRow.length; j < curRow_length; j++) {
                final byte b_base = alternate[j-1]; // letter in b at the current pos
                // in other words, step_diag = sw[i-1][j-1] + wd(a_base,b_base);
                final int step_diag = lastRow[j-1] + (a_base == b_base ? w_match : w_mismatch);

                // optimized "traversal" of all the matrix cells above the current one (i.e. traversing
                // all 'step down' events that would end in the current cell. The optimized code
                // does exactly the same thing as the commented out loop below. IMPORTANT:
                // the optimization works ONLY for linear w(k)=wopen+(k-1)*wextend!!!!

                // if a gap (length 1) was just opened above, this is the cost of arriving to the current cell:
                int prev_gap = lastRow[j] + w_open;
                best_gap_v[j] += w_extend; // for the gaps that were already opened earlier, extending them by 1 costs w_extend
                if (  prev_gap > best_gap_v[j]  ) {
                    // opening a gap just before the current cell results in better score than extending by one
                    // the best previously opened gap. This will hold for ALL cells below: since any gap
                    // once opened always costs w_extend to extend by another base, we will always get a better score
                    // by arriving to any cell below from the gap we just opened (prev_gap) rather than from the previous best gap
                    best_gap_v[j] = prev_gap;
                    gap_size_v[j] = 1; // remember that the best step-down gap from above has length 1 (we just opened it)
                } else {
                    // previous best gap is still the best, even after extension by another base, so we just record that extension:
                    gap_size_v[j]++;
                }

                final int step_down = best_gap_v[j] ;
                final int kd = gap_size_v[j];

                // optimized "traversal" of all the matrix cells to the left of the current one (i.e. traversing
                // all 'step right' events that would end in the current cell. The optimized code
                // does exactly the same thing as the commented out loop below. IMPORTANT:
                // the optimization works ONLY for linear w(k)=wopen+(k-1)*wextend!!!!

                prev_gap =curRow[j-1]  + w_open; // what would it cost us to open length 1 gap just to the left from current cell
                best_gap_h[i] += w_extend; // previous best gap would cost us that much if extended by another base
                if ( prev_gap > best_gap_h[i] ) {
                    // newly opened gap is better (score-wise) than any previous gap with the same row index i; since
                    // gap penalty is linear with k, this new gap location is going to remain better than any previous ones
                    best_gap_h[i] = prev_gap;
                    gap_size_h[i] = 1;
                } else {
                    gap_size_h[i]++;
                }

                final int step_right = best_gap_h[i];
                final int ki = gap_size_h[i];

                //priority here will be step diagonal, step right, step down
                final boolean diagHighestOrEqual = (step_diag >= step_down)
                                                && (step_diag >= step_right);

                if ( diagHighestOrEqual ) {
                    curRow[j]= Math.max(MATRIX_MIN_CUTOFF, step_diag);
                    curBackTrackRow[j]=0;
                }
                else if(step_right>=step_down) { //moving right is the highest
                    curRow[j]= Math.max(MATRIX_MIN_CUTOFF, step_right);
                    curBackTrackRow[j]=-ki; // negative = horizontal
                }
                else  {
                    curRow[j]= Math.max(MATRIX_MIN_CUTOFF, step_down);
                    curBackTrackRow[j]= kd; // positive=vertical
                }
            }
        }
    }

    /*
     * Class to store the result of calculating the CIGAR from the back track matrix
     */
    private static final class SWPairwiseAlignmentResult implements SmithWatermanAlignment {
        private final Cigar cigar;
        private final int alignmentOffset;

        SWPairwiseAlignmentResult(final Cigar cigar, final int alignmentOffset) {
            this.cigar = cigar;
            this.alignmentOffset = alignmentOffset;
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

    /**
     * Calculates the CIGAR for the alignment from the back track matrix
     *
     * @param sw                   the Smith-Waterman matrix to use
     * @param btrack               the back track matrix to use
     * @param overhangStrategy    the strategy to use for dealing with overhangs
     * @return non-null SWPairwiseAlignmentResult object
     */
    private static SWPairwiseAlignmentResult calculateCigar(final int[][] sw, final int[][] btrack, final SWAlignerArguments.OverhangStrategy overhangStrategy) {
        // p holds the position we start backtracking from; we will be assembling a cigar in the backwards order
        int p1 = 0, p2 = 0;

        final int refLength = sw.length-1;
        final int altLength = sw[0].length-1;

        int maxscore = Integer.MIN_VALUE; // sw scores are allowed to be negative
        int segment_length = 0; // length of the segment (continuous matches, insertions or deletions)

        // if we want to consider overhangs as legitimate operators, then just start from the corner of the matrix
        if ( overhangStrategy == SWAlignerArguments.OverhangStrategy.INDEL ) {
            p1 = refLength;
            p2 = altLength;
        } else {
            // look for the largest score on the rightmost column. we use >= combined with the traversal direction
            // to ensure that if two scores are equal, the one closer to diagonal gets picked
            //Note: this is not technically smith-waterman, as by only looking for max values on the right we are
            //excluding high scoring local alignments
            p2=altLength;

            for(int i=1;i<sw.length;i++)  {
               final int curScore = sw[i][altLength];
               if (curScore >= maxscore ) {
                    p1 = i;
                    maxscore = curScore;
               }
            }
            // now look for a larger score on the bottom-most row
            if ( overhangStrategy != SWAlignerArguments.OverhangStrategy.LEADING_INDEL ) {
                final int[] bottomRow=sw[refLength];
                for ( int j = 1 ; j < bottomRow.length; j++) {
                    final int curScore=bottomRow[j];
                    // data_offset is the offset of [n][j]
                    if ( curScore > maxscore ||
                            (curScore == maxscore && Math.abs(refLength - j) < Math.abs(p1 - p2) ) ) {
                        p1 = refLength;
                        p2 = j ;
                        maxscore = curScore;
                        segment_length = altLength - j ; // end of sequence 2 is overhanging; we will just record it as 'M' segment
                    }
                }
            }
        }
        final List<CigarElement> lce = new ArrayList<>(5);
        if ( segment_length > 0 && overhangStrategy == SWAlignerArguments.OverhangStrategy.SOFTCLIP ) {
            lce.add(makeElement(State.CLIP, segment_length));
            segment_length = 0;
        }

        // we will be placing all insertions and deletions into sequence b, so the states are named w/regard
        // to that sequence

        State state = State.MATCH;
        do {
            final int btr = btrack[p1][p2];
            final State new_state;
            int step_length = 1;
            if ( btr > 0 ) {
                new_state = State.DELETION;
                step_length = btr;
            } else if ( btr < 0 ) {
                new_state = State.INSERTION;
                step_length = (-btr);
            } else new_state = State.MATCH; // and step_length =1, already set above

            // move to next best location in the sw matrix:
            switch( new_state ) {
                case MATCH:  p1--; p2--; break; // move back along the diag in the sw matrix
                case INSERTION: p2 -= step_length; break; // move left
                case DELETION:  p1 -= step_length; break; // move up
            }

            // now let's see if the state actually changed:
            if ( new_state == state ) segment_length+=step_length;
            else {
                // state changed, lets emit previous segment, whatever it was (Insertion Deletion, or (Mis)Match).
                lce.add(makeElement(state, segment_length));
                segment_length = step_length;
                state = new_state;
            }
        // next condition is equivalent to  while ( sw[p1][p2] != 0 ) (with modified p1 and/or p2:
        } while ( p1 > 0 && p2 > 0 );

        // post-process the last segment we are still keeping;
        // NOTE: if reads "overhangs" the ref on the left (i.e. if p2>0) we are counting
        // those extra bases sticking out of the ref into the first cigar element if DO_SOFTCLIP is false;
        // otherwise they will be softclipped. For instance,
        // if read length is 5 and alignment starts at offset -2 (i.e. read starts before the ref, and only
        // last 3 bases of the read overlap with/align to the ref), the cigar will be still 5M if
        // DO_SOFTCLIP is false or 2S3M if DO_SOFTCLIP is true.
        // The consumers need to check for the alignment offset and deal with it properly.
        final int alignment_offset;
        if ( overhangStrategy == SWAlignerArguments.OverhangStrategy.SOFTCLIP ) {
            lce.add(makeElement(state, segment_length));
            if ( p2 > 0 ) lce.add(makeElement(State.CLIP, p2));
            alignment_offset = p1;
        } else if ( overhangStrategy == SWAlignerArguments.OverhangStrategy.IGNORE ) {
            lce.add(makeElement(state, segment_length + p2));
            alignment_offset = p1 - p2;
        } else {  // overhangStrategy == OverhangStrategy.INDEL || overhangStrategy == OverhangStrategy.LEADING_INDEL

            // take care of the actual alignment
            lce.add(makeElement(state, segment_length));

            // take care of overhangs at the beginning of the alignment
            if ( p1 > 0 ) {
                lce.add(makeElement(State.DELETION, p1));
            } else if ( p2 > 0 ) {
                lce.add(makeElement(State.INSERTION, p2));
            }

            alignment_offset = 0;
        }

        Collections.reverse(lce);
        return new SWPairwiseAlignmentResult(AlignmentUtils.consolidateCigar(new Cigar(lce)), alignment_offset);
    }

    private static CigarElement makeElement(final State state, final int length) {
        CigarOperator op = null;
        switch (state) {
            case MATCH: op = CigarOperator.M; break;
            case INSERTION: op = CigarOperator.I; break;
            case DELETION: op = CigarOperator.D; break;
            case CLIP: op = CigarOperator.S; break;
        }
        return new CigarElement(length, op);
    }

    @VisibleForTesting
    void printAlignment(final byte[] ref, final byte[] read, final SmithWatermanAlignment alignment) {
        printAlignment(overhangStrategy, ref, read, 100, alignment);
    }

    private static void printAlignment(final SWAlignerArguments.OverhangStrategy overhangStrategy, final byte[] ref, final byte[] read, final int width, final SmithWatermanAlignment alignment) {
        final StringBuilder bread = new StringBuilder();
        final StringBuilder bref = new StringBuilder();
        final StringBuilder match = new StringBuilder();

        int i = 0;
        int j = 0;

        final int offset = alignment.getAlignmentOffset();

        Cigar cigar = alignment.getCigar();

        if ( overhangStrategy != SWAlignerArguments.OverhangStrategy.SOFTCLIP ) {

            // we need to go through all the hassle below only if we do not do softclipping;
            // otherwise offset is never negative
            if ( offset < 0 ) {
                for (  ; j < (-offset) ; j++ ) {
                    bread.append((char)read[j]);
                    bref.append(' ');
                    match.append(' ');
                }
                // at negative offsets, our cigar's first element carries overhanging bases
                // that we have just printed above. Tweak the first element to
                // exclude those bases. Here we create a new list of cigar elements, so the original
                // list/original cigar are unchanged (they are unmodifiable anyway!)

                final List<CigarElement> tweaked = new ArrayList<>();
                tweaked.addAll(cigar.getCigarElements());
                tweaked.set(0,new CigarElement(cigar.getCigarElement(0).getLength()+offset,
                        cigar.getCigarElement(0).getOperator()));
                cigar = new Cigar(tweaked);
            }
        }

        if ( offset > 0 ) { // note: the way this implementation works, cigar will ever start from S *only* if read starts before the ref, i.e. offset = 0
            for (; i < alignment.getAlignmentOffset() ; i++ ) {
                bref.append((char)ref[i]);
                bread.append(' ');
                match.append(' ');
            }
        }
        
        for ( final CigarElement e : cigar.getCigarElements() ) {
            switch (e.getOperator()) {
                case M :
                    for ( int z = 0 ; z < e.getLength() ; z++, i++, j++  ) {
                        bref.append((i<ref.length)?(char)ref[i]:' ');
                        bread.append((j < read.length)?(char)read[j]:' ');
                        match.append( ( i<ref.length && j < read.length ) ? (ref[i] == read[j] ? '.':'*' ) : ' ' );
                    }
                    break;
                case I :
                    for ( int z = 0 ; z < e.getLength(); z++, j++ ) {
                        bref.append('-');
                        bread.append((char)read[j]);
                        match.append('I');
                    }
                    break;
                case S :
                    for ( int z = 0 ; z < e.getLength(); z++, j++ ) {
                        bref.append(' ');
                        bread.append((char)read[j]);
                        match.append('S');
                    }
                    break;
                case D:
                    for ( int z = 0 ; z < e.getLength(); z++ , i++ ) {
                        bref.append((char)ref[i]);
                        bread.append('-');
                        match.append('D');
                    }
                    break;
                default:
                    throw new GATKException("Unexpected Cigar element:" + e.getOperator());
            }
        }
        for ( ; i < ref.length; i++ ) bref.append((char)ref[i]);
        for ( ; j < read.length; j++ ) bread.append((char)read[j]);

        int pos = 0 ;
        final int maxlength = Math.max(match.length(), Math.max(bread.length(), bref.length()));
        while ( pos < maxlength ) {
            printCautiously(match, pos, width);
            printCautiously(bread, pos, width);
            printCautiously(bref, pos, width);
            System.out.println();
            pos += width;
        }
    }

    /** String builder's substring is extremely stupid: instead of trimming and/or returning an empty
     * string when one end/both ends of the interval are out of range, it crashes with an
     * exception. This utility function simply prints the substring if the interval is within the index range
     * or trims accordingly if it is not.
     */
    private static void printCautiously(final StringBuilder s, final int start, final int width) {
        if ( start >= s.length() ) {
            System.out.println();
            return;
        }
        final int end = Math.min(start + width, s.length());
        System.out.println(s.substring(start,end));
    }

}
