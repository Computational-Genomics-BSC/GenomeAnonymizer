package utils;

import genomicelements.GenomicRegion;
import genomicelements.GenomicRegionBaseImpl;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class OperationsTest {

    @Test
    void compareOrdersBySequenceIndexFirst() {
        assertEquals(-3, Operations.compare(0, 100, 200, 1, 1, 50));
        assertEquals(3, Operations.compare(2, 1, 50, 1, 100, 200));
    }

    @Test
    void compareDistinguishesOverlapFromSeparation() {
        // Same sequence, first region ends before: -2 without overlap, -1 with overlap
        assertEquals(-2, Operations.compare(0, 1, 10, 0, 20, 30));
        assertEquals(-1, Operations.compare(0, 1, 25, 0, 20, 30));
        assertEquals(2, Operations.compare(0, 20, 30, 0, 1, 10));
        assertEquals(1, Operations.compare(0, 20, 30, 0, 1, 25));
    }

    @Test
    void compareWithSameEndUsesStart() {
        assertEquals(-1, Operations.compare(0, 5, 30, 0, 20, 30));
        assertEquals(1, Operations.compare(0, 20, 30, 0, 5, 30));
        assertEquals(0, Operations.compare(0, 20, 30, 0, 20, 30));
    }

    @Test
    void compareOfRegionsUsesTheirCoordinates() {
        GenomicRegion first = new GenomicRegionBaseImpl("chr1", 0, 1, 10);
        GenomicRegion second = new GenomicRegionBaseImpl("chr1", 0, 20, 30);
        GenomicRegion other = new GenomicRegionBaseImpl("chr2", 1, 1, 10);
        assertEquals(-2, Operations.compare(first, second));
        assertEquals(-3, Operations.compare(second, other));
        assertEquals(-2, first.compareTo(second));
    }

    @Test
    void overlapIncludesSharedEnds() {
        // The four int version takes both starts first and both ends after
        assertTrue(Operations.overlap(1, 10, 10, 20));
        assertFalse(Operations.overlap(1, 11, 10, 20));
        assertTrue(Operations.overlap(5, 1, 6, 20));
    }

    @Test
    void overlapOfRegions() {
        GenomicRegion first = new GenomicRegionBaseImpl("chr1", 1, 10);
        GenomicRegion second = new GenomicRegionBaseImpl("chr1", 10, 20);
        GenomicRegion third = new GenomicRegionBaseImpl("chr1", 11, 20);
        assertTrue(Operations.overlap(first, second));
        assertFalse(Operations.overlap(first, third));
        assertTrue(Operations.overlap(first, 10, 20));
        assertFalse(Operations.overlap(first, 11, 20));
    }

    @Test
    void overlapOfIntervalAndPoint() {
        assertTrue(Operations.overlap(5, 10, 5));
        assertTrue(Operations.overlap(5, 10, 10));
        assertFalse(Operations.overlap(5, 10, 11));
        assertFalse(Operations.overlap(5, 10, 4));
    }

    @Test
    void isContainedRequiresFirstInsideSecond() {
        GenomicRegion inner = new GenomicRegionBaseImpl("chr1", 5, 10);
        GenomicRegion outer = new GenomicRegionBaseImpl("chr1", 1, 10);
        assertTrue(Operations.isContained(inner, outer));
        assertFalse(Operations.isContained(outer, inner));
        assertTrue(Operations.isContained(inner, inner));
    }

    @Test
    void twoDimEuclideanDistanceTakesEachDimensionAsAPairOfArguments() {
        // Arguments are (x1, x2, y1, y2): distance between (0, 0) and (3, 4)
        assertEquals(5.0, Operations.computeTwoDimEuclideanDistance(0, 3, 0, 4), 1e-9);
        assertEquals(0.0, Operations.computeTwoDimEuclideanDistance(7, 7, 2, 2), 1e-9);
        // Thresholds of AnonymizedReadAlignmentProvider are built this way: 5 of position and 5 of length is 7.07
        assertEquals(7.07, Operations.computeTwoDimEuclideanDistance(100, 105, 3, 8), 0.01);
    }

    @Test
    void threeDimEuclideanDistance() {
        assertEquals(3.0, Operations.computeThreeDimEuclideanDistance(0, 0, 0, 1, 2, 2), 1e-9);
    }

    @Test
    void nucleotideByteToUpperCase() {
        assertEquals((byte) 'A', Operations.nucleotideByteToUpperCase((byte) 'a'));
        assertEquals((byte) 'T', Operations.nucleotideByteToUpperCase((byte) 't'));
        assertEquals((byte) 'G', Operations.nucleotideByteToUpperCase((byte) 'G'));
        assertEquals((byte) 'N', Operations.nucleotideByteToUpperCase((byte) 'n'));
    }
}
