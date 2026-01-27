package utils;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

public class DNAalphabet {
    public static final Set<Character> ALPHABET = new HashSet<>(
            Arrays.asList(
                    'A', 'T', 'C', 'G', 'N'
            )
    );

    public static boolean isValidBase(char base) {
        return ALPHABET.contains(Character.toUpperCase(base));
    }

    public static boolean isValidByteBase(byte base) {
        return ALPHABET.contains((char) Character.toUpperCase(base));
    }
}
