package genomicelements;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class GenomicRegionBaseImplTest {

    @Test
    void toStringIsUsedAsSuffixOfPartitionFiles() {
        assertEquals("chr1_1-1000", new GenomicRegionBaseImpl("chr1", 0, 1, 1000).toString());
    }

    @Test
    void comparesBySequenceIndexAndThenByCoordinates() {
        GenomicRegion first = new GenomicRegionBaseImpl("chr1", 0, 1, 1000);
        GenomicRegion second = new GenomicRegionBaseImpl("chr1", 0, 1001, 2000);
        GenomicRegion other = new GenomicRegionBaseImpl("chr2", 1, 1, 10);
        assertTrue(first.compareTo(second) < 0);
        assertTrue(second.compareTo(first) > 0);
        assertTrue(second.compareTo(other) < 0);
        assertEquals(0, first.compareTo(new GenomicRegionBaseImpl("chr1", 0, 1, 1000)));
    }

    @Test
    void sequenceIndexDefaultsToZero() {
        GenomicRegion region = new GenomicRegionBaseImpl("chr5", 10, 20);
        assertEquals(0, region.getSequenceIdx());
        region.setSequenceIdx(4);
        assertEquals(4, region.getSequenceIdx());
    }
}
