package simulated;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;

/**
 * Runs the GenomeAnonymizer JAR as an external process, exactly as a user would, and collects its exit code and log
 * @author Nicolas Gaitan
 */
public class AnonymizerRunner {

    public static final String NORMAL_OUTPUT_SUFFIX = ".anonymized.N.bam";
    public static final String TUMOR_OUTPUT_SUFFIX = ".anonymized.T.bam";
    public static final String NORMAL_FIXVAF_OUTPUT_SUFFIX = ".anonymized.N.fixVAF.bam";
    public static final String TUMOR_FIXVAF_OUTPUT_SUFFIX = ".anonymized.T.fixVAF.bam";

    public record RunResult(int exitCode, String log, double seconds, String outputPrefix, File logFile) {
        public File normalOutput() {
            return new File(outputPrefix + NORMAL_OUTPUT_SUFFIX);
        }

        public File tumorOutput() {
            return new File(outputPrefix + TUMOR_OUTPUT_SUFFIX);
        }

        public File normalFixVAFOutput() {
            return new File(outputPrefix + NORMAL_FIXVAF_OUTPUT_SUFFIX);
        }

        public File tumorFixVAFOutput() {
            return new File(outputPrefix + TUMOR_FIXVAF_OUTPUT_SUFFIX);
        }
    }

    public static RunResult run(File normal, File tumor, File reference, String outputPrefix, String... extraArgs)
            throws IOException, InterruptedException {
        File jar = SimulatedTestConfig.getJar();
        if (!jar.exists()) {
            throw new IllegalStateException("The GenomeAnonymizer JAR was not found at " + jar.getAbsolutePath()
                    + ". Run ./gradlew jar first, or launch the tests with ./gradlew simulatedTest");
        }
        File outputDir = new File(outputPrefix).getAbsoluteFile().getParentFile();
        outputDir.mkdirs();
        // Remove the outputs of previous runs with the same prefix, so stale files are never validated
        String prefixName = new File(outputPrefix).getName();
        File[] previous = outputDir.listFiles((dir, name) -> name.startsWith(prefixName + ".") || name.startsWith(prefixName + "_"));
        if (previous != null) {
            for (File file : previous) file.delete();
        }
        List<String> command = new ArrayList<>();
        command.add(new File(System.getProperty("java.home"), "bin/java").getAbsolutePath());
        command.add("-Xmx" + SimulatedTestConfig.getMaxHeap());
        command.add("-jar");
        command.add(jar.getAbsolutePath());
        command.add("-in");
        command.add(normal.getAbsolutePath());
        command.add("-it");
        command.add(tumor.getAbsolutePath());
        command.add("-r");
        command.add(reference.getAbsolutePath());
        command.add("-o");
        command.add(new File(outputPrefix).getAbsolutePath());
        command.add("-t");
        command.add(String.valueOf(SimulatedTestConfig.getThreads()));
        command.add("-s");
        command.add(SimulatedTestConfig.DEFAULT_SEED);
        command.addAll(Arrays.asList(extraArgs));
        File logFile = new File(outputPrefix + ".log");
        ProcessBuilder builder = new ProcessBuilder(command);
        builder.redirectErrorStream(true);
        builder.redirectOutput(logFile);
        long start = System.currentTimeMillis();
        Process process = builder.start();
        boolean finished = process.waitFor(SimulatedTestConfig.RUN_TIMEOUT_MINUTES, TimeUnit.MINUTES);
        if (!finished) {
            process.destroyForcibly();
            throw new IllegalStateException("GenomeAnonymizer did not finish in " + SimulatedTestConfig.RUN_TIMEOUT_MINUTES
                    + " minutes. Log: " + logFile.getAbsolutePath());
        }
        long end = System.currentTimeMillis();
        String log = new String(Files.readAllBytes(logFile.toPath()), StandardCharsets.UTF_8);
        return new RunResult(process.exitValue(), log, (double) (end - start) / 1000, new File(outputPrefix).getAbsolutePath(), logFile);
    }
}
