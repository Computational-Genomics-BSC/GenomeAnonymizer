package simulated;

import analysis.GenomeAnonymizer;

import java.io.File;

/**
 * Configuration of the simulated validation tests, read from system properties set by the simulatedTest gradle task.
 * Every property has a default relative to the project directory, so the tests can also be launched from an IDE
 * after running ./gradlew jar
 * @author Nicolas Gaitan
 */
public class SimulatedTestConfig {

    public static final String SIMULATED_TAG = "simulated";

    public static final String DEFAULT_THREADS = "8";
    public static final String DEFAULT_MAX_HEAP = "6g";
    public static final String DEFAULT_SEED = "42";
    // Upper limit for a single run of GenomeAnonymizer on a simulated dataset
    public static final int RUN_TIMEOUT_MINUTES = 15;

    public static File getJar() {
        return new File(System.getProperty("simulated.jar", "build/libs/GenomeAnonymizer-" + GenomeAnonymizer.VERSION + ".jar"));
    }

    public static File getDataDir() {
        return new File(System.getProperty("simulated.dataDir", "src/test/resources/simulated"));
    }

    public static File getWellformedDir() {
        return new File(getDataDir(), "wellformed");
    }

    public static File getMalformedDir() {
        return new File(getDataDir(), "malformed");
    }

    public static File getBaselineDir() {
        return new File(getDataDir(), "baseline");
    }

    public static File getOutputDir() {
        return new File(System.getProperty("simulated.outputDir", "build/simulated-test"));
    }

    public static File getReportDir() {
        return new File(System.getProperty("simulated.reportDir", "build/reports/simulated"));
    }

    public static int getThreads() {
        return Integer.parseInt(System.getProperty("simulated.threads", DEFAULT_THREADS));
    }

    public static String getMaxHeap() {
        return System.getProperty("simulated.maxHeap", DEFAULT_MAX_HEAP);
    }

    public static boolean updateBaseline() {
        return Boolean.parseBoolean(System.getProperty("simulated.updateBaseline", "false"));
    }

    public static boolean strictBaseline() {
        return Boolean.parseBoolean(System.getProperty("simulated.strictBaseline", "false"));
    }
}
