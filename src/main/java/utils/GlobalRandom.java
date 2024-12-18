package utils;

import java.util.Random;

/**
 * Singleton class to provide a global random number generator.
 */
public class GlobalRandom {
    private static Random RANDOM;

    public static Random getInstance() {
        if (RANDOM == null) {
            throw new IllegalStateException("Random not initialized");
        }
        return RANDOM;
    }

    public static void setSeed(int seed) {
        if (RANDOM != null) {
            throw new IllegalStateException("Random already initialized");
        }
        if (seed == -1) {
            seed = (int) System.currentTimeMillis();
        }
        RANDOM = new Random(seed);
    }
}
