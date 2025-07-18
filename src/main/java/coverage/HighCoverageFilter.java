package coverage;

import java.util.LinkedList;
import java.util.Queue;
import java.util.Set;
import java.util.HashSet;

/**
 * Class to manage coverage of genomic regions.
 * It allows adding reads and checking if the coverage exceeds a specified maximum.
 * 
 * @author Rodrigo Martin
 */
public class HighCoverageFilter {
    private final Queue<ReadInfo> reads = new LinkedList<>();
    private final Set<String> keptReads;
    private final int maxCoverage;
    private boolean inHighCoverageRegion = false;

    public HighCoverageFilter(int maxCoverage) {
        this.maxCoverage = maxCoverage;
        this.keptReads = new HashSet<>(maxCoverage*2);
    }

    /**
     * Adds a read to the coverage filter.
     * 
     * Returns true if the read should be excluded due to high coverage.
     */
    public boolean addRead(int readStart, int readEnd, String readName, int pairIndex) {
        // Empty the reads queue of reads that are before the current read
        while (!reads.isEmpty() && reads.peek().getEnd() < readStart) {
            ReadInfo removedRead = reads.poll();
            if (inHighCoverageRegion) {
                // If we are in a high coverage region, we need to remove the read from keptReads
                keptReads.remove(removedRead.getUniqueName());
            }
        }
        // Add the current read to the queue
        ReadInfo currentRead = new ReadInfo(readStart, readEnd, readName, pairIndex);
        reads.add(currentRead);
        boolean highCoverageCheck = reads.size() > maxCoverage;
        if (inHighCoverageRegion) {
            // We are still in a high coverage region
            if (highCoverageCheck) {
                if (keptReads.size() < maxCoverage) {
                    // We are still in a high coverage region, but we keep the read
                    keptReads.add(currentRead.getUniqueName());
                    return false;
                } else {
                    return true;
                }
            } else {
                // We are no longer in a high coverage region
                inHighCoverageRegion = false;
                keptReads.clear();
            }
        } else {
            if (highCoverageCheck) {
                // We have entered a high coverage region
                inHighCoverageRegion = true;
                // Add all the current reads to keptReads
                for (ReadInfo read : reads) {
                    keptReads.add(read.getUniqueName());
                }
            }
        }
        return false;
    }

    /**
     * Internal class representing a read with start, end, and name.
     */
    private final class ReadInfo {
        private final int start;
        private final int end;
        private final String name;
        private final int pairIndex;

        public ReadInfo(int start, int end, String name, int pairIndex) {
            this.start = start;
            this.end = end;
            this.name = name;
            this.pairIndex = pairIndex;
        }

        public int getStart() {
            return start;
        }

        public int getEnd() {
            return end;
        }

        public String getUniqueName() {
            return name + "_" + pairIndex;
        }
    }
}
