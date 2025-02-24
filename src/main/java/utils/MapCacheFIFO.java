package utils;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * @author Nicolas Gaitan
 * Simple Cache extending LinkedHashMap to hold a certain amount of objects and provide First in First Out(FIFO) insertion-removal operations
 * when the default amount of entris is reached
 */
public class MapCacheFIFO<K,V> extends LinkedHashMap<K,V> {
    private static final int DEFAULT_MAX_ENTRIES = 1000;
    private static final int DEFAULT_INITIAL_CAPACITY = 2000;
    private static final float DEFAULT_LOAD_FACTOR = 0.75f;
    private static final boolean DEFAULT_ACCESS_ORDER = false;

    private K lastAddedKey;
    private V lastAddedValue;
    private int maxEntries = DEFAULT_MAX_ENTRIES;

    public MapCacheFIFO(){
        super(DEFAULT_INITIAL_CAPACITY, DEFAULT_LOAD_FACTOR, DEFAULT_ACCESS_ORDER);
    }

    public MapCacheFIFO(int maxEntries){
        super(maxEntries+1, DEFAULT_LOAD_FACTOR, DEFAULT_ACCESS_ORDER);
        this.maxEntries = maxEntries;
    }

    public K getLastAddedKey(){
        return this.lastAddedKey;
    }

    public V getLastAddedValue(){
        return this.lastAddedValue;
    }

    public void putEntry(K key, V value){
        this.put(key, value);
        lastAddedKey = key;
        lastAddedValue = value;
    }

    @Override
    protected boolean removeEldestEntry(Map.Entry eldest){
        return size() > maxEntries;
    }
}
