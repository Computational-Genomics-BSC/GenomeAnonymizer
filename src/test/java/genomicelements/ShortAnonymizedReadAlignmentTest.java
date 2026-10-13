package genomicelements;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import testutils.TestAlignments;

import java.nio.charset.StandardCharsets;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class ShortAnonymizedReadAlignmentTest {

    private static final String CONTIG = "chr1";
    private static final String OTHER_CONTIG = "chr2";
    private static final int CONTIG_LENGTH = 1000;
    private static final int READ_LENGTH = 50;
    // 1-based start of the read under test and of its mate, when the read is the leftmost of the pair
    private static final int START = 101;
    private static final int MATE_START = 301;

    private String reference;
    private SAMFileHeader header;

    @BeforeEach
    void setUp() {
        reference = TestAlignments.randomSequence(CONTIG_LENGTH, 11);
        header = TestAlignments.header(new String[]{CONTIG, OTHER_CONTIG}, new int[]{CONTIG_LENGTH, CONTIG_LENGTH});
    }

    /** Reference bases of the 1-based inclusive interval */
    private String ref(int first, int last) {
        return reference.substring(first - 1, last);
    }

    private SAMRecord leftmostRead(int start, String cigar, String bases) {
        return TestAlignments.read(header, "read1", CONTIG, start, cigar, bases, true, false, CONTIG, MATE_START, true, 250);
    }

    private SAMRecord rightmostRead(int start, String cigar, String bases) {
        return TestAlignments.read(header, "read1", CONTIG, start, cigar, bases, false, true, CONTIG, 1, false, -250);
    }

    private ShortAnonymizedReadAlignment anonymizedRead(SAMRecord record) {
        ShortAnonymizedReadAlignment answer = new ShortAnonymizedReadAlignment(record, false);
        answer.setReferenceContigSequence(reference.getBytes(StandardCharsets.UTF_8));
        return answer;
    }

    private Signal indelSignal(ShortAnonymizedReadAlignment read, int location, int cigarIdx, int length, Signal.IndelSignalType type) {
        Signal signal = new Signal(CONTIG, location, read.getReadAlignmentId(), cigarIdx, length, Signal.Source.INDEL);
        signal.setIndelSignalType(type);
        return signal;
    }

    @Test
    void readWithoutSignalsIsAlreadyAnonymizedAndUnchanged() {
        String bases = ref(START, START + READ_LENGTH - 1);
        ShortAnonymizedReadAlignment read = anonymizedRead(leftmostRead(START, "50M", bases));
        assertTrue(read.isAnonymized());
        SAMRecord record = read.getAnonymizedSamRecord();
        assertEquals(bases, record.getReadString());
        assertEquals("50M", record.getCigarString());
        assertEquals(START, record.getAlignmentStart());
    }

    @Test
    void recordIsNotAvailableUntilPendingSignalsAreApplied() {
        ShortAnonymizedReadAlignment read = anonymizedRead(leftmostRead(START, "50M", ref(START, START + READ_LENGTH - 1)));
        Signal snv = new Signal(CONTIG, START + 10, read.getReadAlignmentId(), 10, 1, Signal.Source.SNV);
        snv.setRefAllele(new byte[]{(byte) reference.charAt(START + 9)});
        read.addSignalToAnonymize(snv);
        assertFalse(read.isAnonymized());
        assertThrows(IllegalStateException.class, read::getAnonymizedSamRecord);
        assertThrows(IllegalStateException.class, read::getPairUpdatedPos);
    }

    @Test
    void anonymizingWithoutReferenceFails() {
        ShortAnonymizedReadAlignment read = new ShortAnonymizedReadAlignment(leftmostRead(START, "50M", ref(START, START + READ_LENGTH - 1)), false);
        assertThrows(IllegalStateException.class, read::anonymizeRead);
    }

    @Test
    void germlineSNVIsReplacedByTheReferenceBase() {
        int snvReadIdx = 10;
        int snvPosition = START + snvReadIdx;
        char refBase = reference.charAt(snvPosition - 1);
        char altBase = TestAlignments.otherBase(refBase);
        SAMRecord record = leftmostRead(START, "50M", ref(START, START + READ_LENGTH - 1));
        TestAlignments.setBaseAtReferencePosition(record, snvPosition, altBase);
        // A second mismatch that is not marked as germline has to survive
        int keptReadIdx = 30;
        char keptAltBase = TestAlignments.otherBase(reference.charAt(START + keptReadIdx - 1));
        TestAlignments.setBaseAtReferencePosition(record, START + keptReadIdx, keptAltBase);
        ShortAnonymizedReadAlignment read = anonymizedRead(record);
        Signal snv = new Signal(CONTIG, snvPosition, read.getReadAlignmentId(), snvReadIdx, 1, Signal.Source.SNV);
        snv.setRefAllele(new byte[]{(byte) refBase});
        snv.setAltAllele(new byte[]{(byte) altBase});
        read.addSignalToAnonymize(snv);
        read.anonymizeRead();
        assertTrue(read.isAnonymized());
        SAMRecord anonymized = read.getAnonymizedSamRecord();
        assertEquals(refBase, anonymized.getReadString().charAt(snvReadIdx));
        assertEquals(keptAltBase, anonymized.getReadString().charAt(keptReadIdx));
        assertEquals("50M", anonymized.getCigarString());
        assertEquals(START, anonymized.getAlignmentStart());
        assertEquals(1, anonymized.getIntegerAttribute("NM"));
        assertEquals(TestAlignments.DEFAULT_BASE_QUALITY, anonymized.getBaseQualities()[snvReadIdx]);
        assertFalse(read.updateInfoForMate());
    }

    @Test
    void germlineDeletionIsFilledWithReferenceAndTheReadKeepsItsLength() {
        // 20 bases, 3 deleted from the reference, 30 bases
        String bases = ref(START, START + 19) + ref(START + 23, START + 52);
        ShortAnonymizedReadAlignment read = anonymizedRead(leftmostRead(START, "20M3D30M", bases));
        read.addSignalToAnonymize(indelSignal(read, START + 19, 1, 3, Signal.IndelSignalType.DELETION));
        read.anonymizeRead();
        SAMRecord anonymized = read.getAnonymizedSamRecord();
        assertEquals("50M", anonymized.getCigarString());
        assertEquals(START, anonymized.getAlignmentStart());
        assertEquals(ref(START, START + READ_LENGTH - 1), anonymized.getReadString());
        assertEquals(READ_LENGTH, anonymized.getBaseQualities().length);
        assertEquals(0, anonymized.getIntegerAttribute("NM"));
        assertFalse(read.updateInfoForMate());
    }

    @Test
    void deletionFilledOverAmbiguousReferenceBaseGetsN() {
        StringBuilder ambiguous = new StringBuilder(reference);
        ambiguous.setCharAt(START + 20, 'R');
        String bases = ref(START, START + 19) + ref(START + 23, START + 52);
        ShortAnonymizedReadAlignment read = new ShortAnonymizedReadAlignment(leftmostRead(START, "20M3D30M", bases), false);
        read.setReferenceContigSequence(ambiguous.toString().getBytes(StandardCharsets.UTF_8));
        read.addSignalToAnonymize(indelSignal(read, START + 19, 1, 3, Signal.IndelSignalType.DELETION));
        read.anonymizeRead();
        assertEquals('N', read.getAnonymizedSamRecord().getReadString().charAt(21));
    }

    @Test
    void germlineInsertionIsRemovedAndTheLeftmostReadIsExtendedToTheLeft() {
        String bases = ref(START, START + 19) + "ACGT" + ref(START + 20, START + 49);
        ShortAnonymizedReadAlignment read = anonymizedRead(leftmostRead(START, "20M4I30M", bases));
        read.addSignalToAnonymize(indelSignal(read, START + 19, 1, 4, Signal.IndelSignalType.INSERTION));
        read.anonymizeRead();
        SAMRecord anonymized = read.getAnonymizedSamRecord();
        assertEquals("54M", anonymized.getCigarString());
        assertEquals(START - 4, anonymized.getAlignmentStart());
        assertEquals(ref(START - 4, START + 49), anonymized.getReadString());
        assertEquals(0, anonymized.getIntegerAttribute("NM"));
        // The mate has to be told that this read moved
        assertTrue(read.updateInfoForMate());
        assertEquals(START - 4, read.getPairUpdatedPos());
    }

    @Test
    void germlineInsertionIsRemovedAndTheRightmostReadIsExtendedToTheRight() {
        String bases = ref(START, START + 19) + "ACGT" + ref(START + 20, START + 49);
        ShortAnonymizedReadAlignment read = anonymizedRead(rightmostRead(START, "20M4I30M", bases));
        read.addSignalToAnonymize(indelSignal(read, START + 19, 1, 4, Signal.IndelSignalType.INSERTION));
        read.anonymizeRead();
        SAMRecord anonymized = read.getAnonymizedSamRecord();
        assertEquals("54M", anonymized.getCigarString());
        assertEquals(START, anonymized.getAlignmentStart());
        assertEquals(ref(START, START + 53), anonymized.getReadString());
        assertFalse(read.updateInfoForMate());
    }

    @Test
    void germlineSoftClipAtTheStartOfTheLeftmostReadIsFilledWithReference() {
        int clip = 6;
        String bases = "GGGGGG" + ref(START + clip, START + 49);
        ShortAnonymizedReadAlignment read = anonymizedRead(leftmostRead(START + clip, "6S44M", bases));
        read.addSignalToAnonymize(new Signal(CONTIG, START + clip, read.getReadAlignmentId(), 0, clip, Signal.Source.SOFT_CLIP));
        read.anonymizeRead();
        SAMRecord anonymized = read.getAnonymizedSamRecord();
        assertEquals("50M", anonymized.getCigarString());
        assertEquals(START, anonymized.getAlignmentStart());
        assertEquals(ref(START, START + 49), anonymized.getReadString());
        assertTrue(read.updateInfoForMate());
    }

    @Test
    void germlineSoftClipAtTheEndOfTheRightmostReadIsFilledWithReference() {
        int clip = 6;
        String bases = ref(START, START + 43) + "GGGGGG";
        ShortAnonymizedReadAlignment read = anonymizedRead(rightmostRead(START, "44M6S", bases));
        read.addSignalToAnonymize(new Signal(CONTIG, START + 43, read.getReadAlignmentId(), 1, clip, Signal.Source.SOFT_CLIP));
        read.anonymizeRead();
        SAMRecord anonymized = read.getAnonymizedSamRecord();
        assertEquals("50M", anonymized.getCigarString());
        assertEquals(START, anonymized.getAlignmentStart());
        assertEquals(ref(START, START + 49), anonymized.getReadString());
        assertFalse(read.updateInfoForMate());
    }

    @Test
    void alternativeAlignmentTagsAreRemovedFromAnonymizedReads() {
        SAMRecord record = leftmostRead(START, "50M", ref(START, START + READ_LENGTH - 1));
        TestAlignments.setBaseAtReferencePosition(record, START + 10, TestAlignments.otherBase(reference.charAt(START + 9)));
        record.setAttribute("XA", "chr2,+500,50M,1;");
        record.setAttribute("XS", 40);
        ShortAnonymizedReadAlignment read = anonymizedRead(record);
        Signal snv = new Signal(CONTIG, START + 10, read.getReadAlignmentId(), 10, 1, Signal.Source.SNV);
        snv.setRefAllele(new byte[]{(byte) reference.charAt(START + 9)});
        read.addSignalToAnonymize(snv);
        read.anonymizeRead();
        assertNull(read.getAnonymizedSamRecord().getAttribute("XA"));
        assertNull(read.getAnonymizedSamRecord().getAttribute("XS"));
    }

    @Test
    void strandOrientationSignalTurnsThePairIntoForwardReverse() {
        // Both reads of the pair in the forward strand
        SAMRecord record = TestAlignments.read(header, "read1", CONTIG, START, "50M", ref(START, START + READ_LENGTH - 1),
                true, false, CONTIG, MATE_START, false, 250);
        ShortAnonymizedReadAlignment read = anonymizedRead(record);
        read.addSignalToAnonymize(new Signal(CONTIG, START, read.getReadAlignmentId(), 0, 250, Signal.Source.STRAND_ORIENTATION));
        read.anonymizeRead();
        assertTrue(read.fixOrientation());
        assertFalse(read.hasDestructiveSignal());
        SAMRecord anonymized = read.getAnonymizedSamRecord();
        assertFalse(anonymized.getReadNegativeStrandFlag());
        assertTrue(anonymized.getMateNegativeStrandFlag());
        assertEquals(250, anonymized.getInferredInsertSize());
    }

    @Test
    void correctOrientationOfTheRightmostRead() {
        SAMRecord record = TestAlignments.read(header, "read1", CONTIG, MATE_START, "50M", ref(MATE_START, MATE_START + READ_LENGTH - 1),
                false, false, CONTIG, START, false, 250);
        ShortAnonymizedReadAlignment.correctOrientation(record);
        assertTrue(record.getReadNegativeStrandFlag());
        assertFalse(record.getMateNegativeStrandFlag());
        assertEquals(-250, record.getInferredInsertSize());
    }

    @Test
    void correctOrientationOfReadsWithTheSameStartDependsOnWhichIsFoundFirst() {
        SAMRecord first = TestAlignments.read(header, "read1", CONTIG, START, "50M", ref(START, START + READ_LENGTH - 1),
                true, true, CONTIG, START, true, 50);
        SAMRecord second = TestAlignments.read(header, "read1", CONTIG, START, "50M", ref(START, START + READ_LENGTH - 1),
                false, true, CONTIG, START, true, 50);
        ShortAnonymizedReadAlignment.correctOrientation(first, true);
        ShortAnonymizedReadAlignment.correctOrientation(second);
        assertFalse(first.getReadNegativeStrandFlag());
        assertTrue(first.getMateNegativeStrandFlag());
        assertTrue(second.getReadNegativeStrandFlag());
        assertFalse(second.getMateNegativeStrandFlag());
    }

    @Test
    void chromosomeChangeSignalReplacesTheMateWithASyntheticOneFromTheReference() {
        int insertSizeMedian = 300;
        SAMRecord record = TestAlignments.read(header, "read1", CONTIG, START, "50M", ref(START, START + READ_LENGTH - 1),
                true, false, OTHER_CONTIG, 200, true, 0);
        record.setProperPairFlag(false);
        ShortAnonymizedReadAlignment read = anonymizedRead(record);
        read.addSignalToAnonymize(new Signal(CONTIG, START, read.getReadAlignmentId(), 0, 0, Signal.Source.CHROM_CHANGE));
        read.anonymizeRead();
        assertTrue(read.hasDestructiveSignal());
        SAMRecord newMate = read.getNewPair(insertSizeMedian);
        SAMRecord anonymized = read.getAnonymizedSamRecord();
        // The new mate is placed downstream so the fragment has the median insert size
        int expectedMateStart = START + insertSizeMedian - READ_LENGTH;
        assertEquals("read1", newMate.getReadName());
        assertEquals(CONTIG, newMate.getContig());
        assertEquals(expectedMateStart, newMate.getAlignmentStart());
        assertEquals("50M", newMate.getCigarString());
        assertEquals(ref(expectedMateStart, expectedMateStart + READ_LENGTH - 1), newMate.getReadString());
        assertEquals(TestAlignments.DEFAULT_BASE_QUALITY, newMate.getBaseQualities()[0]);
        assertTrue(newMate.getSecondOfPairFlag());
        assertTrue(newMate.getReadNegativeStrandFlag());
        assertTrue(newMate.getProperPairFlag());
        assertEquals(START, newMate.getMateAlignmentStart());
        assertEquals(-insertSizeMedian, newMate.getInferredInsertSize());
        assertEquals(0, newMate.getIntegerAttribute("NM"));
        // The original read now points to the new mate
        assertEquals(CONTIG, anonymized.getMateReferenceName());
        assertEquals(expectedMateStart, anonymized.getMateAlignmentStart());
        assertTrue(anonymized.getProperPairFlag());
        assertFalse(anonymized.getReadNegativeStrandFlag());
        assertTrue(anonymized.getMateNegativeStrandFlag());
        assertEquals(insertSizeMedian, anonymized.getInferredInsertSize());
        // Both carry the temporary tag used to choose which copy of the pair is kept when merging
        assertEquals(ShortAnonymizedReadAlignment.PAIR_1_IDX, newMate.getIntegerAttribute(ShortAnonymizedReadAlignment.ORIGIN_PAIR_TAG));
        assertEquals(ShortAnonymizedReadAlignment.PAIR_1_IDX, anonymized.getIntegerAttribute(ShortAnonymizedReadAlignment.ORIGIN_PAIR_TAG));
    }

    @Test
    void supplementaryAlignmentsAreMarkedForEliminationWithDestructiveSignals() {
        SAMRecord record = TestAlignments.read(header, "read1", CONTIG, START, "50M", ref(START, START + READ_LENGTH - 1),
                true, false, OTHER_CONTIG, 200, true, 0);
        record.setAttribute("SA", "chr2,700,+,20S30M,60,0;");
        ShortAnonymizedReadAlignment read = anonymizedRead(record);
        read.addSignalToAnonymize(new Signal(CONTIG, START, read.getReadAlignmentId(), 0, 0, Signal.Source.CHROM_CHANGE));
        read.anonymizeRead();
        assertTrue(read.eliminateSupplementaries());
        assertNull(read.getAnonymizedSamRecord().getAttribute("SA"));
    }

    @Test
    void readAlignmentIdDistinguishesTheReadsOfAPair() {
        SAMRecord[] pair = TestAlignments.properPair(header, "read1", CONTIG, reference, START, MATE_START, READ_LENGTH);
        String firstId = ShortAnonymizedReadAlignment.generateReadAlnId(pair[0]);
        String secondId = ShortAnonymizedReadAlignment.generateReadAlnId(pair[1]);
        assertNotEquals(firstId, secondId);
        assertTrue(firstId.startsWith("read1;0;" + START));
        assertTrue(secondId.startsWith("read1;1;" + MATE_START));
        assertEquals(ShortAnonymizedReadAlignment.PAIR_1_IDX, new ShortAnonymizedReadAlignment(pair[0]).getPairIdx());
        assertEquals(ShortAnonymizedReadAlignment.PAIR_2_IDX, new ShortAnonymizedReadAlignment(pair[1]).getPairIdx());
    }
}
