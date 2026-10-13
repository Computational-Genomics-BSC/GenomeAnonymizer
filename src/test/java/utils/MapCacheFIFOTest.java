package utils;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class MapCacheFIFOTest {

    @Test
    void removesTheOldestEntryWhenFull() {
        MapCacheFIFO<Integer, String> cache = new MapCacheFIFO<>(3);
        for (int i = 1; i <= 4; i++) cache.put(i, "v" + i);
        assertEquals(3, cache.size());
        assertFalse(cache.containsKey(1));
        assertTrue(cache.containsKey(2));
        assertTrue(cache.containsKey(4));
    }

    @Test
    void readingAnEntryDoesNotProtectItFromRemoval() {
        MapCacheFIFO<Integer, String> cache = new MapCacheFIFO<>(2);
        cache.put(1, "v1");
        cache.put(2, "v2");
        cache.get(1);
        cache.put(3, "v3");
        assertFalse(cache.containsKey(1));
        assertTrue(cache.containsKey(2));
    }

    @Test
    void putEntryRemembersTheLastAddition() {
        MapCacheFIFO<Integer, String> cache = new MapCacheFIFO<>();
        cache.putEntry(7, "seven");
        cache.putEntry(8, "eight");
        assertEquals(8, cache.getLastAddedKey());
        assertEquals("eight", cache.getLastAddedValue());
        assertEquals("seven", cache.get(7));
    }
}
