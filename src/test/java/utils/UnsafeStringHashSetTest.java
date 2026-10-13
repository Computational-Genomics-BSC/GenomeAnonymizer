package utils;

import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class UnsafeStringHashSetTest {

    @Test
    void addAndContains() {
        UnsafeStringHashSet set = new UnsafeStringHashSet();
        assertTrue(set.isEmpty());
        assertTrue(set.add("A00741:213:HSIMNDSX2:1:1101:1000:2000"));
        assertFalse(set.add("A00741:213:HSIMNDSX2:1:1101:1000:2000"));
        assertTrue(set.contains("A00741:213:HSIMNDSX2:1:1101:1000:2000"));
        assertFalse(set.contains("A00741:213:HSIMNDSX2:1:1101:1000:2001"));
        assertFalse(set.contains(42));
        assertEquals(1, set.size());
    }

    @Test
    void growsBeyondInitialCapacityWithoutLosingElements() {
        UnsafeStringHashSet set = new UnsafeStringHashSet(4);
        int n = 10_000;
        for (int i = 0; i < n; i++) {
            assertTrue(set.add("read_" + i));
        }
        assertEquals(n, set.size());
        for (int i = 0; i < n; i++) {
            assertTrue(set.contains("read_" + i), "Missing read_" + i);
        }
        assertFalse(set.contains("read_" + n));
    }

    @Test
    void removeKeepsOtherElementsReachable() {
        UnsafeStringHashSet set = new UnsafeStringHashSet(8);
        for (int i = 0; i < 100; i++) set.add("read_" + i);
        for (int i = 0; i < 100; i += 2) assertTrue(set.remove("read_" + i));
        assertFalse(set.remove("read_0"));
        assertEquals(50, set.size());
        for (int i = 0; i < 100; i++) {
            assertEquals(i % 2 == 1, set.contains("read_" + i), "read_" + i);
        }
    }

    @Test
    void addAllFromAnotherSetAndFromCollection() {
        UnsafeStringHashSet first = new UnsafeStringHashSet();
        first.add("a");
        first.add("b");
        UnsafeStringHashSet second = new UnsafeStringHashSet();
        second.add("b");
        second.add("c");
        assertTrue(first.addAll(second));
        assertFalse(first.addAll(second));
        assertEquals(3, first.size());
        assertTrue(first.containsAll(Arrays.asList("a", "b", "c")));
        assertTrue(first.addAll(Arrays.asList("c", "d")));
        assertEquals(4, first.size());
    }

    @Test
    void clearEmptiesTheSet() {
        UnsafeStringHashSet set = new UnsafeStringHashSet();
        set.add("a");
        set.clear();
        assertTrue(set.isEmpty());
        assertFalse(set.contains("a"));
    }

    @Test
    void unsupportedOperationsAndInvalidArguments() {
        UnsafeStringHashSet set = new UnsafeStringHashSet();
        assertThrows(UnsupportedOperationException.class, set::iterator);
        assertThrows(UnsupportedOperationException.class, set::toArray);
        assertThrows(NullPointerException.class, () -> set.add(null));
        assertThrows(IllegalArgumentException.class, () -> new UnsafeStringHashSet(0));
    }
}
