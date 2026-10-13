package genomicelements;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class MultiSourceCIGARSignalTest {

    private Signal deletion() {
        Signal signal = new Signal("chr1", 100, "readA;0;90", 1, 5, Signal.Source.INDEL);
        signal.setIndelSignalType(Signal.IndelSignalType.DELETION);
        signal.setIsFromNormalDataset(true);
        return signal;
    }

    @Test
    void startsWithTheReadOfTheWrappedSignal() {
        MultiSourceCIGARSignal multiSource = new MultiSourceCIGARSignal(deletion());
        assertEquals(1, multiSource.getSources().size());
        assertEquals(1, multiSource.getSources().get("readA;0;90"));
        assertEquals(Signal.IndelSignalType.DELETION, multiSource.getIndelSignalType());
        assertTrue(multiSource.isIndel());
        assertTrue(multiSource.isIntraAlignmentSignal());
        assertFalse(multiSource.isInterAlignmentSignal());
        assertFalse(multiSource.isDestructiveSignal());
        assertTrue(multiSource.isFromNormalDataset());
    }

    @Test
    void addSourceKeepsTheCigarIndexOfEachRead() {
        MultiSourceCIGARSignal multiSource = new MultiSourceCIGARSignal(deletion());
        multiSource.addSource("readB;1;60", 3);
        assertEquals(2, multiSource.getSources().size());
        assertEquals(3, multiSource.getSources().get("readB;1;60"));
    }

    @Test
    void readSpecificSignalIsACopyForThatRead() {
        Signal wrapped = deletion();
        MultiSourceCIGARSignal multiSource = new MultiSourceCIGARSignal(wrapped);
        multiSource.addSource("readB;1;60", 3);
        Signal readSignal = multiSource.getReadAlnSignal("readB;1;60", 3);
        assertEquals("readB;1;60", readSignal.getReadAlnName());
        assertEquals(3, readSignal.getInReadPosition());
        assertEquals(wrapped.getLocation(), readSignal.getLocation());
        assertEquals(wrapped.getLength(), readSignal.getLength());
        assertEquals(Signal.IndelSignalType.DELETION, readSignal.getIndelSignalType());
        // The wrapped signal is not modified
        assertEquals("readA;0;90", wrapped.getReadAlnName());
        assertEquals(1, wrapped.getInReadPosition());
    }

    @Test
    void classificationStateIsSharedWithTheWrappedSignal() {
        Signal wrapped = deletion();
        MultiSourceCIGARSignal multiSource = new MultiSourceCIGARSignal(wrapped);
        multiSource.setIsGermline(true);
        multiSource.setHandled(true);
        multiSource.setClassifiedByDistance(true);
        assertTrue(wrapped.isGermline());
        assertTrue(wrapped.isHandled());
        assertTrue(wrapped.isClassifiedByDistance());
        assertEquals(wrapped.getSignalKey(), multiSource.getSignalKey());
    }
}
