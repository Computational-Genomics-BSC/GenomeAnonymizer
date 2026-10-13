package io;

import genomicelements.GenomicRegionBaseImpl;
import genomicelements.PairedPileup;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import testutils.TestAlignments;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class SamplePairReadAlignmentReaderTest {

    private static final String CONTIG = "chr1";
    private static final int CONTIG_LENGTH = 1000;
    private static final int READ_LENGTH = 50;

    @TempDir
    File tempDir;

    @Test
    void pairsThePileupsOfBothSamplesByPosition() throws IOException {
        String reference = TestAlignments.randomSequence(CONTIG_LENGTH, 3);
        SAMFileHeader header = TestAlignments.header(new String[]{CONTIG}, new int[]{CONTIG_LENGTH});
        File referenceFile = TestAlignments.writeReference(tempDir, "ref.fa", new String[]{CONTIG}, new String[]{reference});
        // Normal covers 101-150 and 301-350, tumor covers 131-180 and 331-380
        List<SAMRecord> normalRecords = new ArrayList<>(List.of(TestAlignments.properPair(header, "n1", CONTIG, reference, 101, 301, READ_LENGTH)));
        List<SAMRecord> tumorRecords = new ArrayList<>(List.of(TestAlignments.properPair(header, "t1", CONTIG, reference, 131, 331, READ_LENGTH)));
        File normal = TestAlignments.writeIndexedBam(tempDir, "normal.bam", header, normalRecords);
        File tumor = TestAlignments.writeIndexedBam(tempDir, "tumor.bam", header, tumorRecords);
        try (SamplePairReadAlignmentReader reader = new SamplePairReadAlignmentReader(normal.getPath(), tumor.getPath(),
                referenceFile.getPath(), reference.getBytes(StandardCharsets.UTF_8), new GenomicRegionBaseImpl(CONTIG, 0, 1, CONTIG_LENGTH))) {
            int previous = 0;
            int onlyNormal = 0;
            int onlyTumor = 0;
            int both = 0;
            for (PairedPileup pileup : reader) {
                int location = pileup.getLocation();
                assertTrue(location > previous, "Pileups out of order at " + location);
                previous = location;
                boolean inNormal = (location >= 101 && location <= 150) || (location >= 301 && location <= 350);
                boolean inTumor = (location >= 131 && location <= 180) || (location >= 331 && location <= 380);
                if (inNormal) assertNotNull(pileup.getNormalPileup(), "Normal pileup at " + location);
                else assertNull(pileup.getNormalPileup(), "Normal pileup at " + location);
                if (inTumor) assertNotNull(pileup.getTumorPileup(), "Tumor pileup at " + location);
                else assertNull(pileup.getTumorPileup(), "Tumor pileup at " + location);
                if (inNormal && inTumor) both++;
                else if (inNormal) onlyNormal++;
                else onlyTumor++;
                assertEquals(inNormal && !inTumor, pileup.getNormalPileup() != null && pileup.isOnlyNormal());
            }
            assertEquals(60, onlyNormal);
            assertEquals(60, onlyTumor);
            assertEquals(40, both);
        }
    }

    @Test
    void iterationIsRestrictedToTheRegion() throws IOException {
        String reference = TestAlignments.randomSequence(CONTIG_LENGTH, 3);
        SAMFileHeader header = TestAlignments.header(new String[]{CONTIG}, new int[]{CONTIG_LENGTH});
        File referenceFile = TestAlignments.writeReference(tempDir, "ref.fa", new String[]{CONTIG}, new String[]{reference});
        List<SAMRecord> records = new ArrayList<>(List.of(TestAlignments.properPair(header, "r1", CONTIG, reference, 101, 301, READ_LENGTH)));
        File normal = TestAlignments.writeIndexedBam(tempDir, "normal.bam", header, records);
        File tumor = TestAlignments.writeIndexedBam(tempDir, "tumor.bam", header, records);
        try (SamplePairReadAlignmentReader reader = new SamplePairReadAlignmentReader(normal.getPath(), tumor.getPath(),
                referenceFile.getPath(), reference.getBytes(StandardCharsets.UTF_8), new GenomicRegionBaseImpl(CONTIG, 0, 250, 320))) {
            int count = 0;
            for (PairedPileup pileup : reader) {
                assertFalse(pileup.getLocation() < 301 || pileup.getLocation() > 320, "Pileup outside the region at " + pileup.getLocation());
                count++;
            }
            assertEquals(20, count);
        }
    }
}
