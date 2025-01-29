package utils;

/**
 * Holder class for functions to operate on genomic elements, such as reads or variants
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

    public static boolean overlap(int first1, int first2, int last1, int last2){
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
}
