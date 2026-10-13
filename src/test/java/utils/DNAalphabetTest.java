package utils;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class DNAalphabetTest {

    @Test
    void acceptsTheFourBasesAndNInAnyCase() {
        for (char base : "ACGTNacgtn".toCharArray()) {
            assertTrue(DNAalphabet.isValidBase(base), "Base " + base);
            assertTrue(DNAalphabet.isValidByteBase((byte) base), "Byte base " + base);
        }
    }

    @Test
    void rejectsOtherIUPACCodesAndSymbols() {
        for (char base : "RYKMSWBDHVU*-.".toCharArray()) {
            assertFalse(DNAalphabet.isValidBase(base), "Base " + base);
            assertFalse(DNAalphabet.isValidByteBase((byte) base), "Byte base " + base);
        }
    }
}
