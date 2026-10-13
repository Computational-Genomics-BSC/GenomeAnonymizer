package analysis;

import analysis.ShortReadAnonymizer.ModifiedPairsRequired;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import testutils.TestAlignments;
import utils.GlobalRandom;

import java.util.HashMap;

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
}
