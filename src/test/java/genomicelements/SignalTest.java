package genomicelements;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class SignalTest {

    private Signal signal(Signal.Source source) {
        return new Signal("chr1", 100, "read;0;90", 3, 5, source);
    }

    @Test
    void intraAlignmentSignals() {
        for (Signal.Source source : new Signal.Source[]{Signal.Source.SNV, Signal.Source.INDEL, Signal.Source.SOFT_CLIP}) {
            Signal signal = signal(source);
            assertTrue(signal.isIntraAlignmentSignal(), source.getName());
            assertFalse(signal.isInterAlignmentSignal(), source.getName());
            assertFalse(signal.isDestructiveSignal(), source.getName());
        }
        assertTrue(signal(Signal.Source.SNV).isSNV());
        assertTrue(signal(Signal.Source.INDEL).isIndel());
        assertTrue(signal(Signal.Source.SOFT_CLIP).isSoftClip());
    }

    @Test
    void interAlignmentSignals() {
        for (Signal.Source source : new Signal.Source[]{Signal.Source.INSERT_SIZE, Signal.Source.STRAND_ORIENTATION, Signal.Source.CHROM_CHANGE}) {
            Signal signal = signal(source);
            assertTrue(signal.isInterAlignmentSignal(), source.getName());
            assertFalse(signal.isIntraAlignmentSignal(), source.getName());
        }
    }

    @Test
    void onlyInsertSizeAndChromChangeAreDestructive() {
        assertTrue(signal(Signal.Source.INSERT_SIZE).isDestructiveSignal());
        assertTrue(signal(Signal.Source.CHROM_CHANGE).isDestructiveSignal());
        assertFalse(signal(Signal.Source.STRAND_ORIENTATION).isDestructiveSignal());
    }

    @Test
    void signalIsAGenomicRegionOfOneBase() {
        Signal signal = signal(Signal.Source.SNV);
        assertEquals("chr1", signal.getSequenceName());
        assertEquals(100, signal.getLocation());
        assertEquals(100, signal.getStart());
        assertEquals(100, signal.getEnd());
    }

    @Test
    void keyOfIntraAlignmentSignalIsSharedByReadsWithTheSameEvent() {
        Signal first = new Signal("chr1", 100, "readA;0;90", 1, 5, Signal.Source.INDEL);
        first.setIndelSignalType(Signal.IndelSignalType.DELETION);
        Signal second = new Signal("chr1", 100, "readB;1;60", 2, 5, Signal.Source.INDEL);
        second.setIndelSignalType(Signal.IndelSignalType.DELETION);
        assertEquals("chr1:100:5:INDEL:DELETION", first.getSignalKey());
        assertEquals(first.getSignalKey(), second.getSignalKey());
        second.setIndelSignalType(Signal.IndelSignalType.INSERTION);
        assertNotEquals(first.getSignalKey(), second.getSignalKey());
    }

    @Test
    void keyOfInterAlignmentSignalIsTheReadAlignment() {
        Signal signal = new Signal("chr1", 100, "readA;0;90", 0, 900, Signal.Source.INSERT_SIZE);
        assertEquals("readA;0;90", signal.getSignalKey());
    }

    @Test
    void defaultState() {
        Signal signal = signal(Signal.Source.INDEL);
        assertFalse(signal.isGermline());
        assertFalse(signal.isHandled());
        assertFalse(signal.isClassifiedByDistance());
        assertFalse(signal.comesFromDuplicateRead());
        assertEquals(Signal.IndelSignalType.NOT_INDEL, signal.getIndelSignalType());
        signal.setIsFromNormalDataset(true);
        assertTrue(signal.isFromNormalDataset());
        assertFalse(signal.isFromTumoralDataset());
    }

    @Test
    void equalsComparesLocationReadAndSource() {
        Signal first = signal(Signal.Source.INDEL);
        Signal same = signal(Signal.Source.INDEL);
        Signal otherSource = signal(Signal.Source.SOFT_CLIP);
        Signal otherRead = new Signal("chr1", 100, "other;0;90", 3, 5, Signal.Source.INDEL);
        assertEquals(first, same);
        assertNotEquals(first, otherSource);
        assertNotEquals(first, otherRead);
    }

    @Test
    void cloneCanBeAssignedToAnotherRead() throws CloneNotSupportedException {
        Signal original = signal(Signal.Source.INDEL);
        original.setIndelSignalType(Signal.IndelSignalType.INSERTION);
        Signal clone = (Signal) original.clone();
        assertNotSame(original, clone);
        clone.setReadAlnName("other;1;200");
        clone.setInReadPosition(7);
        assertEquals("read;0;90", original.getReadAlnName());
        assertEquals(3, original.getInReadPosition());
        assertEquals(original.getSignalKey(), clone.getSignalKey());
        assertEquals(Signal.IndelSignalType.INSERTION, clone.getIndelSignalType());
    }
}
