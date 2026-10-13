package analysis;

import analysis.ShortReadAnonymizer.ModifiedPairsRequired;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import testutils.TestAlignments;
import utils.GlobalRandom;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class ShortReadAnonymizerTest {

    @BeforeAll
    static void initializeRandom() {
        // The generator can only be seeded once per JVM, and the constructor of the anonymizer needs it
        try {
            GlobalRandom.getInstance();
        }
        catch (IllegalStateException e) {
            GlobalRandom.setSeed(1);
        }
    }

    @TempDir
    File tempDir;

    private ShortReadAnonymizer anonymizer() {
        return new ShortReadAnonymizer("normal.bam", "tumor.bam", "ref.fa", "out");
    }

    @Test
    void randomGeneratorCannotBeSeededTwice() {
        assertThrows(IllegalStateException.class, () -> GlobalRandom.setSeed(2));
    }

    @Test
    void outputNamesIdentifyTheDataset() {
        ShortReadAnonymizer anonymizer = anonymizer();
        assertEquals("prefix.anonymized.N.bam", anonymizer.getBAMOutputName("prefix", ShortReadAnonymizer.NORMAL_DATASET_IDX));
        assertEquals("prefix.anonymized.T.bam", anonymizer.getBAMOutputName("prefix", ShortReadAnonymizer.TUMORAL_DATASET_IDX));
    }

    @Test
    void pairedReadNameMatchesTheOneOfTheRecords() {
        ShortReadAnonymizer anonymizer = anonymizer();
        SAMFileHeader header = TestAlignments.header(new String[]{"chr1"}, new int[]{1000});
        SAMRecord[] pair = TestAlignments.properPair(header, "read1", "chr1", TestAlignments.randomSequence(1000, 1), 101, 301, 50);
        // The mate positions collected per partition are looked up at merge time with SAMRecord.getPairedReadName
        assertEquals(pair[0].getPairedReadName(), anonymizer.getPairedReadName("read1", 0));
        assertEquals(pair[1].getPairedReadName(), anonymizer.getPairedReadName("read1", 1));
    }

    @Test
    void unmodifiedPairIsKept() {
        ShortReadAnonymizer anonymizer = anonymizer();
        SAMFileHeader header = TestAlignments.header(new String[]{"chr1"}, new int[]{1000});
        SAMRecord[] pair = TestAlignments.properPair(header, "read1", "chr1", TestAlignments.randomSequence(1000, 1), 101, 301, 50);
        assertTrue(anonymizer.keepPair(pair[0], new HashMap<>(), new HashMap<>()));
        assertTrue(anonymizer.keepPair(pair[1], new HashMap<>(), new HashMap<>()));
    }

    @Test
    void readWhoseMateWasDiscardedIsNotKept() {
        ShortReadAnonymizer anonymizer = anonymizer();
        SAMFileHeader header = TestAlignments.header(new String[]{"chr1"}, new int[]{1000});
        SAMRecord[] pair = TestAlignments.properPair(header, "read1", "chr1", TestAlignments.randomSequence(1000, 1), 101, 301, 50);
        pair[0].setMateAlignmentStart(-1);
        assertFalse(anonymizer.keepPair(pair[0], new HashMap<>(), new HashMap<>()));
    }

    @Test
    void modifiedPairNeedsTheNewPairAndTheFirstOriginalMate() {
        ModifiedPairsRequired required = new ModifiedPairsRequired();
        assertFalse(required.collectedAllPairs());
        required.addFirstRighmostPair();
        assertFalse(required.collectedAllPairs());
        assertFalse(required.collectedSecondLeftmost());
        required.addSecondLeftmostPair();
        assertTrue(required.collectedSecondLeftmost());
        assertTrue(required.collectedAllPairs());
    }

    @Test
    void modifiedPairWithTwoGeneratedPairsNeedsBothOriginalMates() {
        ModifiedPairsRequired required = new ModifiedPairsRequired();
        required.addFirstRighmostPair();
        required.addSecondLeftmostPair();
        required.setNeedsBothSecondPairs();
        assertFalse(required.collectedAllPairs());
        required.addSecondRighmostPair();
        assertTrue(required.collectedAllPairs());
    }

    @Test
    void mitochondrialReadsTheirPairsAndSupplementariesAreLeftOutOfTheOutput() throws IOException {
        String contig = "chr1";
        String mitochondrialContig = "chrM";
        int contigLength = 2000;
        int readLength = 50;
        String reference = TestAlignments.randomSequence(contigLength, 3);
        String mitochondrialReference = TestAlignments.randomSequence(contigLength, 4);
        String[] contigNames = {contig, mitochondrialContig};
        SAMFileHeader header = TestAlignments.header(contigNames, new int[]{contigLength, contigLength});
        File referenceFile = TestAlignments.writeReference(tempDir, "ref.fa", contigNames, new String[]{reference, mitochondrialReference});
        List<SAMRecord> records = new ArrayList<>();
        String cigar = readLength + "M";
        // The only pair that must reach the output
        SAMRecord[] nuclearPair = TestAlignments.properPair(header, "nuclear", contig, reference, 101, 301, readLength);
        SAMRecord[] mitochondrialPair = TestAlignments.properPair(header, "mito", mitochondrialContig, mitochondrialReference, 101, 301, readLength);
        // A nuclear pair with a supplementary alignment in the mitochondrial chromosome
        SAMRecord[] nuclearPairWithSupplementary = TestAlignments.properPair(header, "nuclearSupp", contig, reference, 601, 801, readLength);
        SAMRecord mitochondrialSupplementary = TestAlignments.read(header, "nuclearSupp", mitochondrialContig, 601, "20M30S",
                mitochondrialReference.substring(600, 620) + reference.substring(620, 650), true, false, contig, 801, true, 0);
        mitochondrialSupplementary.setSupplementaryAlignmentFlag(true);
        // A mitochondrial pair with a supplementary alignment in a nuclear chromosome
        SAMRecord[] mitochondrialPairWithSupplementary = TestAlignments.properPair(header, "mitoSupp", mitochondrialContig, mitochondrialReference,
                1101, 1301, readLength);
        SAMRecord nuclearSupplementary = TestAlignments.read(header, "mitoSupp", contig, 1101, "20M30S",
                reference.substring(1100, 1120) + mitochondrialReference.substring(1120, 1150), true, false, mitochondrialContig, 1301, true, 0);
        nuclearSupplementary.setSupplementaryAlignmentFlag(true);
        for (SAMRecord[] pair : new SAMRecord[][]{nuclearPair, mitochondrialPair, nuclearPairWithSupplementary, mitochondrialPairWithSupplementary}) {
            records.add(pair[0]);
            records.add(pair[1]);
        }
        records.add(mitochondrialSupplementary);
        records.add(nuclearSupplementary);
        // A pair with one read in each chromosome
        records.add(TestAlignments.read(header, "split", contig, 1501, cigar, reference.substring(1500, 1550), true, false,
                mitochondrialContig, 1501, true, 0));
        records.add(TestAlignments.read(header, "split", mitochondrialContig, 1501, cigar, mitochondrialReference.substring(1500, 1550), false, true,
                contig, 1501, false, 0));
        File normalFile = TestAlignments.writeIndexedBam(tempDir, "normal.bam", header, records);
        File tumorFile = TestAlignments.writeIndexedBam(tempDir, "tumor.bam", header, records);
        String outputPrefix = new File(tempDir, "out").getPath();
        ShortReadAnonymizer anonymizer = new ShortReadAnonymizer(normalFile.getPath(), tumorFile.getPath(), referenceFile.getPath(), outputPrefix);
        anonymizer.setPartitions();
        anonymizer.queryReadsToExclude();
        anonymizer.anonymizeReads();
        for (int datasetIdx : new int[]{ShortReadAnonymizer.NORMAL_DATASET_IDX, ShortReadAnonymizer.TUMORAL_DATASET_IDX}) {
            int count = 0;
            try (SamReader reader = SamReaderFactory.makeDefault().open(new File(anonymizer.getBAMOutputName(outputPrefix, datasetIdx)))) {
                for (SAMRecord record : reader) {
                    assertEquals(contig, record.getReferenceName());
                    assertEquals(contig, record.getMateReferenceName());
                    assertFalse(record.isSecondaryOrSupplementary());
                    assertTrue(record.getAlignmentStart() == 101 || record.getAlignmentStart() == 301,
                            "Unexpected read at " + record.getAlignmentStart());
                    count++;
                }
            }
            assertEquals(2, count);
        }
    }
}
