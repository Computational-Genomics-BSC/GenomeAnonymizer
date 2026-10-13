package simulated;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Runs the GenomeAnonymizer JAR over the malformed simulated datasets. Each case has exactly one defect, and
 * src/test/resources/simulated/malformed_cases.tsv declares for each one whether the run must fail and the text
 * that the log must contain. The file records the behaviour of the build that was accepted last: when the
 * handling of a defect is changed on purpose (a clearer message, a new input check), the file has to be edited.
 * @author Nicolas Gaitan
 */
@Tag(SimulatedTestConfig.SIMULATED_TAG)
public class MalformedDatasetTest {

    public static final String CASES_FILE = "malformed_cases.tsv";
    public static final String EXPECTED_FAIL = "FAIL";
    public static final String EXPECTED_PASS = "PASS";
    public static final String FATAL_ERROR_MESSAGE = "Fatal error happened";

    static List<String> caseNames() throws IOException {
        return new ArrayList<>(readCases().getRows().keySet());
    }

    private static Table readCases() throws IOException {
        return Table.read(new File(SimulatedTestConfig.getDataDir(), CASES_FILE));
    }

    @Test
    void everyMalformedDatasetIsDeclared() throws IOException {
        Set<String> declared = new TreeSet<>(readCases().getRows().keySet());
        Set<String> present = new TreeSet<>();
        File[] directories = SimulatedTestConfig.getMalformedDir().listFiles(File::isDirectory);
        assertTrue(directories != null && directories.length > 0, "No malformed datasets in " + SimulatedTestConfig.getMalformedDir());
        for (File directory : directories) present.add(directory.getName());
        assertEquals(present, declared, "The cases of " + CASES_FILE + " and the directories of malformed datasets differ");
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("caseNames")
    void malformedDatasetBehavesAsDeclared(String caseName) throws Exception {
        Map<String, String> expectations = readCases().getRow(caseName);
        File caseDir = new File(SimulatedTestConfig.getMalformedDir(), caseName);
        String outputPrefix = new File(new File(SimulatedTestConfig.getOutputDir(), "malformed"), caseName).getPath();
        AnonymizerRunner.RunResult run = AnonymizerRunner.run(new File(caseDir, SimulatedDataset.NORMAL_FILE),
                new File(caseDir, SimulatedDataset.TUMOR_FILE), new File(caseDir, SimulatedDataset.REFERENCE_FILE), outputPrefix);
        String context = "Case " + caseName + " (" + expectations.get("defect") + "). Log: " + run.logFile();
        if (EXPECTED_PASS.equals(expectations.get("expected"))) {
            assertEquals(0, run.exitCode(), "The run should finish correctly. " + context);
            assertTrue(run.normalOutput().exists() && run.tumorOutput().exists(), "The outputs should be written. " + context);
            return;
        }
        assertEquals(1, run.exitCode(), "The run should fail with exit code 1. " + context);
        assertTrue(run.log().contains(FATAL_ERROR_MESSAGE), "The log should report a fatal error. " + context);
        assertTrue(run.log().contains(expectations.get("expected_message")),
                "The log should contain '" + expectations.get("expected_message") + "'. " + context);
        // A failed run must not leave anything that can be mistaken for a result
        assertFalse(run.normalOutput().exists(), "A failed run should not write the normal output. " + context);
        assertFalse(run.tumorOutput().exists(), "A failed run should not write the tumor output. " + context);
    }
}
