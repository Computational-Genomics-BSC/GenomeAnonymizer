package io;

import genomicelements.GenomicRegion;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

public class GenomicRegionBedReaderTest {

    @TempDir
    File tempDir;

    @Test
    void readsRegionsInFileOrderIgnoringExtraColumns() throws IOException {
        File bed = new File(tempDir, "regions.bed");
        try (PrintWriter writer = new PrintWriter(bed)) {
            writer.println("chr1\t100\t200");
            writer.println("chr2\t5\t50\tgeneA\t0\t+");
        }
        List<GenomicRegion> regions = GenomicRegionBedReader.readGenomicRegionBED(bed.getPath());
        assertEquals(2, regions.size());
        assertEquals("chr1", regions.get(0).getSequenceName());
        // Coordinates are taken as written in the file, without conversion
        assertEquals(100, regions.get(0).getStart());
        assertEquals(200, regions.get(0).getEnd());
        assertEquals("chr2", regions.get(1).getSequenceName());
        assertEquals(5, regions.get(1).getStart());
        assertEquals(50, regions.get(1).getEnd());
    }

    @Test
    void failsWithNonNumericCoordinates() throws IOException {
        File bed = new File(tempDir, "broken.bed");
        try (PrintWriter writer = new PrintWriter(bed)) {
            writer.println("chr1\tstart\tend");
        }
        assertThrows(NumberFormatException.class, () -> GenomicRegionBedReader.readGenomicRegionBED(bed.getPath()));
    }
}
