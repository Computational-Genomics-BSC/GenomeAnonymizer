package genomicelements;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import testutils.TestAlignments;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class LocusPileupIteratorTest {

    private static final String CONTIG = "chr1";
    private static final int CONTIG_LENGTH = 1000;
    private static final int READ_LENGTH = 50;
    // Two pairs: the first reads overlap between 121 and 150, the second reads between 321 and 350
    private static final int FIRST_START = 101;
    private static final int SECOND_START = 121;
    private static final int MATE_OFFSET = 200;
    private static final int SNV_POSITION = 130;

    @TempDir
    File tempDir;

    private String reference;
    private SAMFileHeader header;
    private List<SAMRecord> records;

    @BeforeEach
    void setUp() {
        reference = TestAlignments.randomSequence(CONTIG_LENGTH, 5);
        header = TestAlignments.header(new String[]{CONTIG}, new int[]{CONTIG_LENGTH});
        records = new ArrayList<>();
        SAMRecord[] firstPair = TestAlignments.properPair(header, "pairA", CONTIG, reference, FIRST_START, FIRST_START + MATE_OFFSET, READ_LENGTH);
        SAMRecord[] secondPair = TestAlignments.properPair(header, "pairB", CONTIG, reference, SECOND_START, SECOND_START + MATE_OFFSET, READ_LENGTH);
        // Only the first read of pairB differs from the reference
        TestAlignments.setBaseAtReferencePosition(secondPair[0], SNV_POSITION, TestAlignments.otherBase(reference.charAt(SNV_POSITION - 1)));
        records.add(firstPair[0]);
        records.add(firstPair[1]);
        records.add(secondPair[0]);
        records.add(secondPair[1]);
    }

    private SamReader openReader() {
        File bam = TestAlignments.writeIndexedBam(tempDir, "reads.bam", header, records);
        return SamReaderFactory.makeDefault().open(bam);
    }

    private TreeMap<Integer, LocusPileUp> collect(LocusPileupIterator pileups) {
        TreeMap<Integer, LocusPileUp> answer = new TreeMap<>();
        for (LocusPileUp pileup : pileups) {
            assertFalse(answer.containsKey(pileup.getLocation()), "Position " + pileup.getLocation() + " returned twice");
            answer.put(pileup.getLocation(), pileup);
        }
        return answer;
    }

    @Test
    void returnsOnePileupPerCoveredPositionInOrder() throws IOException {
        try (SamReader reader = openReader()) {
            LocusPileupIterator pileups = new LocusPileupIterator(reader, CONTIG, 1, CONTIG_LENGTH);
            int previous = 0;
            int count = 0;
            for (LocusPileUp pileup : pileups) {
                assertTrue(pileup.getLocation() > previous, "Pileups out of order at " + pileup.getLocation());
                assertEquals(CONTIG, pileup.getSequenceName());
                previous = pileup.getLocation();
                count++;
            }
            // 101-170 and 301-370
            assertEquals(140, count);
        }
    }

    @Test
    void claimAllModeGivesEveryReadAtEachPosition() throws IOException {
        try (SamReader reader = openReader()) {
            TreeMap<Integer, LocusPileUp> pileups = collect(new LocusPileupIterator(reader, CONTIG, 1, CONTIG_LENGTH));
            assertEquals(1, pileups.get(FIRST_START).size());
            assertEquals(2, pileups.get(SNV_POSITION).size());
            assertEquals(1, pileups.get(160).size());
            assertFalse(pileups.containsKey(200));
            List<PileupRead> reads = pileups.get(SNV_POSITION).getAllPileupReads();
            assertEquals(2, reads.size());
            Set<String> names = new HashSet<>();
            for (PileupRead read : reads) {
                names.add(read.getReadName());
                assertEquals(SNV_POSITION, read.getLocation());
                // 0-based position inside the read
                assertEquals(SNV_POSITION - read.getStart(), read.getReadPosition());
            }
            assertEquals(Set.of("pairA", "pairB"), names);
        }
    }

    @Test
    void claimNewAndDifferingModeGivesNewReadsAndMismatchesOnly() throws IOException {
        try (SamReader reader = openReader()) {
            TreeMap<Integer, LocusPileUp> pileups = collect(new LocusPileupIterator(reader, CONTIG, 1, CONTIG_LENGTH,
                    reference.getBytes(StandardCharsets.UTF_8)));
            // Depth is still counted with every read
            assertEquals(2, pileups.get(SNV_POSITION).size());
            for (LocusPileUp pileup : pileups.values()) {
                List<PileupRead> claimed = pileup.claimReadsOnPileup();
                int location = pileup.getLocation();
                if (location == FIRST_START || location == SECOND_START || location == FIRST_START + MATE_OFFSET
                        || location == SECOND_START + MATE_OFFSET) {
                    assertEquals(1, claimed.size(), "Position " + location);
                    assertTrue(claimed.get(0).isNew(), "Position " + location);
                    assertFalse(claimed.get(0).differsFromReferenceAtPileup(), "Position " + location);
                }
                else if (location == SNV_POSITION) {
                    assertEquals(1, claimed.size());
                    PileupRead read = claimed.get(0);
                    assertEquals("pairB", read.getReadName());
                    assertFalse(read.isNew());
                    assertTrue(read.differsFromReferenceAtPileup());
                    assertEquals(reference.charAt(SNV_POSITION - 1), read.getReferenceBase());
                    assertEquals(TestAlignments.otherBase(reference.charAt(SNV_POSITION - 1)), read.getBaseAtPileup());
                }
                else {
                    assertTrue(claimed.isEmpty(), "Position " + location);
                }
            }
        }
    }

    @Test
    void newAndDifferingModeRequiresTheReference() throws IOException {
        try (SamReader reader = openReader()) {
            assertThrows(IllegalStateException.class, () -> new LocusPileupIterator(reader, CONTIG, 1, CONTIG_LENGTH, (byte[]) null));
        }
    }

    @Test
    void duplicatesAreExcludedUnlessRequested() throws IOException {
        records.get(2).setDuplicateReadFlag(true);
        try (SamReader reader = openReader()) {
            TreeMap<Integer, LocusPileUp> pileups = collect(new LocusPileupIterator(reader, CONTIG, 1, CONTIG_LENGTH));
            assertEquals(1, pileups.get(SNV_POSITION).size());
        }
        try (SamReader reader = openReader()) {
            LocusPileupIterator iterator = new LocusPileupIterator(reader, CONTIG, 1, CONTIG_LENGTH);
            iterator.setIncludeDuplicates(true);
            TreeMap<Integer, LocusPileUp> pileups = collect(iterator);
            assertEquals(2, pileups.get(SNV_POSITION).size());
            // Duplicates count for the size of the pileup but not for its proper size
            assertEquals(1, pileups.get(SNV_POSITION).getPropperSize());
        }
    }

    @Test
    void excludedReadNamesLowMappingQualityAndUnmappedMatesAreSkipped() throws IOException {
        records.get(1).setMappingQuality(5);
        records.get(3).setMateUnmappedFlag(true);
        Set<String> readsToExclude = new HashSet<>();
        readsToExclude.add("pairB");
        try (SamReader reader = openReader()) {
            LocusPileupIterator iterator = new LocusPileupIterator(reader, CONTIG, 1, CONTIG_LENGTH);
            iterator.setMinimumMappingQuality(20);
            iterator.setReadsToExclude(readsToExclude);
            TreeMap<Integer, LocusPileUp> pileups = collect(iterator);
            // Only the first read of pairA is left
            assertEquals(READ_LENGTH, pileups.size());
            assertEquals(FIRST_START, pileups.firstKey());
            assertEquals(FIRST_START + READ_LENGTH - 1, pileups.lastKey());
        }
    }

    @Test
    void deletedAndInsertedBasesDoNotMakePileups() throws IOException {
        records.clear();
        // 20 aligned bases, 5 deleted, 4 inserted and 26 aligned
        String bases = reference.substring(100, 120) + "ACGT" + reference.substring(125, 151);
        records.add(TestAlignments.read(header, "indel", CONTIG, 101, "20M5D4I26M", bases, true, false, CONTIG, 301, true, 250));
        try (SamReader reader = openReader()) {
            TreeMap<Integer, LocusPileUp> pileups = collect(new LocusPileupIterator(reader, CONTIG, 1, CONTIG_LENGTH));
            assertEquals(46, pileups.size());
            assertTrue(pileups.containsKey(120));
            assertFalse(pileups.containsKey(121));
            assertFalse(pileups.containsKey(125));
            assertTrue(pileups.containsKey(126));
            // The first base after the indels is the 25th of the read (0-based 24)
            assertEquals(24, pileups.get(126).getAllPileupReads().get(0).getReadPosition());
        }
    }

    @Test
    void readsWithOtherIUPACCodesAreRejected() throws IOException {
        byte[] bases = records.get(0).getReadBases();
        bases[10] = 'R';
        records.get(0).setReadBases(bases);
        try (SamReader reader = openReader()) {
            LocusPileupIterator pileups = new LocusPileupIterator(reader, CONTIG, 1, CONTIG_LENGTH);
            assertThrows(IllegalArgumentException.class, () -> collect(pileups));
        }
    }
}
