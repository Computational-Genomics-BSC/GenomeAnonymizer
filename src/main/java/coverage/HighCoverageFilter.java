package coverage;

import java.util.LinkedList;
import java.util.Queue;

/**
 * Class to manage coverage of genomic regions.
 * It allows adding reads and checking if the coverage exceeds a specified maximum.
 * 
 * @author Rodrigo Martin
 */
public class HighCoverageFilter {
    private final Queue<ReadInfo> reads = new LinkedList<>();
    private final int maxCoverage;

    public HighCoverageFilter(int maxCoverage) {
        this.maxCoverage = maxCoverage;
    }

    /**
     * Adds a read to the coverage filter.
     * 
     * Returns true if the read should be excluded due to high coverage.
     */
    public boolean addRead(int readStart, int readEnd) {
        if (reads.isEmpty()) {
            reads.add(new ReadInfo(readStart, readEnd));
            return false;
        }
        // Empty the reads queue of reads that are before the current read
        while (!reads.isEmpty() && reads.peek().getEnd() < readStart) {
            reads.poll();
        }
        // Add the current read to the queue
        ReadInfo currentRead = new ReadInfo(readStart, readEnd);
        reads.add(currentRead);
        return reads.size() > maxCoverage;
    }

    /**
     * Internal class representing a read with start, end, and name.
     */
    private final class ReadInfo {
        private final int start;
        private final int end;

        public ReadInfo(int start, int end) {
            this.start = start;
            this.end = end;
        }

        public int getStart() {
            return start;
        }

        public int getEnd() {
            return end;
        }
    }
}
