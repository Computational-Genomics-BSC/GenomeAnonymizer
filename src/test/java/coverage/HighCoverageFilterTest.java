package coverage;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class HighCoverageFilterTest {

    @Test
    void keepsEveryReadBelowTheMaximumDepth() {
        HighCoverageFilter filter = new HighCoverageFilter(10);
        for (int i = 0; i < 10; i++) {
            assertFalse(filter.addRead(100 + i, 250 + i, "read" + i, 0));
        }
    }

    @Test
    void excludesReadsOnceTheMaximumDepthIsExceeded() {
        int maxDepth = 3;
        HighCoverageFilter filter = new HighCoverageFilter(maxDepth);
        // The read that makes the depth exceed the maximum is still kept, together with the ones before it
        for (int i = 0; i <= maxDepth; i++) {
            assertFalse(filter.addRead(100 + i, 250 + i, "read" + i, 0), "read" + i);
        }
        // Any other read overlapping the pile is excluded
        assertTrue(filter.addRead(110, 260, "extra1", 0));
        assertTrue(filter.addRead(120, 270, "extra2", 1));
    }

    @Test
    void keepsReadsAgainAfterLeavingTheHighDepthRegion() {
        int maxDepth = 3;
        HighCoverageFilter filter = new HighCoverageFilter(maxDepth);
        for (int i = 0; i <= maxDepth; i++) filter.addRead(100 + i, 250 + i, "read" + i, 0);
        assertTrue(filter.addRead(110, 260, "extra", 0));
        assertFalse(filter.addRead(1000, 1150, "far", 0));
        assertFalse(filter.addRead(1001, 1151, "far2", 0));
    }

    @Test
    void bothReadsOfAPairAreCountedSeparately() {
        HighCoverageFilter filter = new HighCoverageFilter(2);
        assertFalse(filter.addRead(100, 250, "pair", 0));
        assertFalse(filter.addRead(100, 250, "pair", 1));
        assertFalse(filter.addRead(101, 251, "other", 0));
        assertTrue(filter.addRead(102, 252, "another", 0));
    }
}
