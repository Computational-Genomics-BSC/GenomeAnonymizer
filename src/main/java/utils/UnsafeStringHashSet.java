package utils;

import java.util.Set;
import java.util.Collection;
import java.util.Arrays;
import java.nio.charset.StandardCharsets;

/**
 * A simple hash set implementation for strings using longs for hashing.
 * It does not store the actual strings, only their hashes.
 * 
 * This implementation is designed for memory efficiency and speed,
 * but it does not handle hash collisions (i.e., two different strings producing the same hash).
 * 
 * This implementation is not thread-safe.
 * 
 * @author Rodrigo Martin
 */
public class UnsafeStringHashSet implements Set<String> {
    private static final int DEFAULT_CAPACITY = 16;
    private static final double LOAD_FACTOR = 0.75;
    
    private long[] table;
    private int size;
    private int threshold;
    private final int initialCapacity;
    
    public UnsafeStringHashSet() {
        this(DEFAULT_CAPACITY);
    }

    public UnsafeStringHashSet(int initialCapacity) {
        if (initialCapacity <= 0) {
            throw new IllegalArgumentException("Initial capacity must be greater than 0");
        }
        this.table = new long[initialCapacity];
        this.initialCapacity = initialCapacity;
        this.size = 0;
        this.threshold = (int) (initialCapacity * LOAD_FACTOR);
    }

    public boolean add(String s) {
        if (s == null) {
            throw new NullPointerException("Cannot add null to UnsafeStringHashSet");
        }
        long hash = computeHash(s);
        return addRaw(hash);
    }

    public boolean addRaw(long hash) {
        int index = indexFor(hash, table.length);

        while (table[index] != 0) {
            if (table[index] == hash) {
                return false; // Already present
            }
            index = (index + 1) % table.length; // Linear probing
        }
        table[index] = hash;
        size++;
        if (size >= threshold) {
            resize();
        }
        return true;
    }

    public boolean addAll(UnsafeStringHashSet other) {
        boolean modified = false;
        for (long hash : other.getRawHashes()) {
            if (hash != 0 && addRaw(hash)) {
                modified = true;
            }
        }
        return modified;
    }

    public boolean addAll(Collection<? extends String> c) {
        boolean modified = false;
        for (String s : c) {
            if (add(s)) {
                modified = true;
            }
        }
        return modified;
    }

    public void clear() {
        table = new long[initialCapacity];
        size = 0;
        threshold = (int) (initialCapacity * LOAD_FACTOR);
    }

    public boolean contains(Object o) {
        if (!(o instanceof String)) {
            return false;
        }
        String s = (String) o;
        long hash = computeHash(s);
        int index = indexFor(hash, table.length);
        
        while (table[index] != 0) {
            if (table[index] == hash) {
                return true; // Found
            }
            index = (index + 1) % table.length; // Linear probing
        }
        return false; // Not found
    }

    public boolean containsAll(Collection<?> c) {
        for (Object o : c) {
            if (!contains(o)) {
                return false;
            }
        }
        return true;
    }

    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Set)) return false;
        Set<?> other = (Set<?>) o;
        if (size() != other.size()) return false;
        for (Object item : other) {
            if (!contains(item)) {
                return false; // If any item is not contained, sets are not equal
            }
        }
        return true;
    }

    public int hashCode() {
        int hash = 0;
        for (long h : table) {
            if (h != 0) {
                hash += Long.hashCode(h);
            }
        }
        return hash;
    }

    public boolean isEmpty() {
        return size == 0;
    }

    public java.util.Iterator<String> iterator() {
        throw new UnsupportedOperationException("iterator not supported");
    }

    public boolean remove(Object o) {
        if (!(o instanceof String)) {
            return false;
        }
        String s = (String) o;
        long hash = computeHash(s);
        int index = indexFor(hash, table.length);
        
        while (table[index] != 0) {
            if (table[index] == hash) {
                table[index] = 0; // Remove
                size--;
                return true;
            }
            index = (index + 1) % table.length; // Linear probing
        }
        return false; // Not found
    }

    public boolean removeAll(Collection<?> c) {
        boolean modified = false;
        for (Object o : c) {
            if (remove(o)) {
                modified = true;
            }
        }
        return modified;
    }
    
    public boolean retainAll(Collection<?> c) {
        throw new UnsupportedOperationException("retainAll not supported");
    }

    public int size() {
        return size;
    }

    public Object[] toArray() {
        throw new UnsupportedOperationException("toArray() not supported");
    }

    public <T> T[] toArray(T[] a) {
        throw new UnsupportedOperationException("toArray(T[] a) not supported");
    }

    public long[] getRawHashes() {
        return Arrays.stream(table).filter(h -> h != 0).toArray();
    }

    private long computeHash(String s) {
        // MurmurHash3-like hash function for strings
        byte[] data = s.getBytes(StandardCharsets.UTF_8);
        final long seed = 0x7f3a21eaL; // Can be any seed value
        final long m = 0xc6a4a7935bd1e995L;
        final int r = 47;
        long h = seed ^ (data.length * m);

        for (int i = 0; i < data.length; i++) {
            long k = data[i];
            k *= m;
            k ^= k >>> r;
            k *= m;
            h ^= k;
            h *= m;
        }

        h ^= h >>> r;
        h *= m;
        h ^= h >>> r;
        // Ensure the hash is non-zero
        if (h == 0) {
            h = 1; // Avoid zero hash
        }
        return h;
    }

    private int indexFor(long hash, int length) {
        return (int) Math.floorMod(hash, length);
    }

    private void resize() {
        int newCapacity = table.length * 2;
        long[] newTable = new long[newCapacity];
        threshold = (int) (newCapacity * LOAD_FACTOR);
        
        for (long hash : table) {
            if (hash != 0) {
                int index = indexFor(hash, newCapacity);
                while (newTable[index] != 0) {
                    index = (index + 1) % newCapacity; // Linear probing
                }
                newTable[index] = hash;
            }
        }
        table = newTable;
    }
}
