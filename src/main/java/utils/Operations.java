package utils;

import genomicelements.GenomicRegion;
import htsjdk.samtools.SAMRecord;

import java.util.Map;

/**
 * Holder class for functions to operate on genomic elements, such as reads or variants
 *
 * @author Nicolas Gaitan
 */
public class Operations {

    // public static Operations getInstance(){
    //     return new Operations();
    //}

    public static int compare(int seqIdx1, int first1, int last1, int seqIdx2, int first2, int last2) {
        boolean ovp = overlap(first1, first2, last1, last2);
        if (seqIdx1 < seqIdx2) {
            return -3;
        }
        if (seqIdx1 > seqIdx2) {
            return 3;
        }
        // For these cases seqIdx1 == seqIdx2
        if (last1 < last2) {
            return ovp ? -1 : -2;
        }
        if (last2 < last1) {
            return ovp ? 1 : 2;
        }
        // For these cases last1 == last2
        if (first1 < first2) {
            return -1;
        }
        if (first2 < first1) {
            return 1;
        }
        return 0;
    }

    public static boolean overlap(GenomicRegion  region1, GenomicRegion region2) {
        int first1 = region1.getStart();
        int first2 = region2.getStart();
        int last1 = region1.getEnd();
        int last2 = region2.getEnd();
        return overlap(first1, first2, last1, last2);
    }

    public static boolean overlap(GenomicRegion region1, int first2, int last2) {
        int first1 = region1.getStart();
        int last1 = region1.getEnd();
        return overlap(first1, first2, last1, last2);
    }

    public static boolean overlap(int first1, int first2, int last1, int last2) {
        return first2 <= last1 && last2 >= first1;
    }

    public static double computeTwoDimEuclideanDistance(int x1, int y1, int x2, int y2) {
        double xDist = Math.pow(x1 - y1, 2);
        double yDist = Math.pow(x2 - y2, 2);
        double squareSum = xDist + yDist;
        return Math.sqrt(squareSum);
    }

    public static double computeThreeDimEuclideanDistance(int x1, int y1, int z1, int x2, int y2, int z2) {
        double squareSum = Math.pow(x1 - x2, 2) + Math.pow(y1 - y2, 2) + Math.pow(z1 - z2, 2);
        return Math.sqrt(squareSum);
    }

    //Subtracts 32 to an ASCII representation of a nucleotide to get the Uppercase ASCII code
    public static byte nucleotideByteToUpperCase(byte nucleotide) {
        return (byte) (nucleotide & 0xDF);
    }

    /**
     * Generates a unique name for a read alignment.
     *
     * @param record The SAMRecord object representing the read alignment.
     * @return A unique string identifier for the alignment.
     */
    public static String generateUniqueAlignmentName(SAMRecord record) {
        // Start with the read name
        StringBuilder uniqueName = new StringBuilder(record.getReadName());

        // Add pair information (is first of pair or second of pair)
        if (record.getReadPairedFlag()) {
            uniqueName.append("_").append(record.getFirstOfPairFlag() ? "1" : "2");
        }

        // Add supplementary alignment information
        if (record.getSupplementaryAlignmentFlag()) {
            uniqueName.append("_SUP");
        }

        // Add the reference name and position on the genome
        uniqueName.append("_")
                .append(record.getReferenceName())
                .append(":")
                .append(record.getAlignmentStart());

        // Optionally, include the strand information
        uniqueName.append(record.getReadNegativeStrandFlag() ? "_" : "_+");

        // Optionally, include the CIGAR string for further specificity
        uniqueName.append("_").append(record.getCigar().toString());

        return uniqueName.toString();
    }
}
