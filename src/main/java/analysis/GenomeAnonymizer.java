package analysis;

import utils.GlobalRandom;
import org.apache.commons.cli.*;
import java.io.File;
import java.io.IOException;
import java.util.logging.ConsoleHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

/**
 * Main class that executes the anonymization method on sequencing data
 * @author Nicolas Gaitan
 */
public class GenomeAnonymizer {

    public static final String VERSION = "1.0.0";
    private static final Logger LOGGER = logConfigure();

    public static final String DEFAULT_RUN_MODE_FUNCTIONALITY = "default";

    public final static int DEFAULT_MIN_MAPPING_QUALITY = 0;
    public static final int DEFAULT_MAX_READS_IN_RAM = 500_000;
    public static final int DEFAULT_MIN_MAX_DEPTH = 100;

    public final static String BAM_FILE = ".bam";
    public final static String SAM_FILE = ".sam";
    public final static String CRAM_FILE = ".cram";

    //Optional arguments attributes
    private int minMappingQuality = DEFAULT_MIN_MAPPING_QUALITY;
    private File tmpDir;
    private int randomSeed;
    private int maxReadsInMemory;
    private int maxDepth = -1;

    /**
     * Run anonymizer with the benchmark of somatic variants functionality. Any somatic variant will be sparred from anonymization
     * @param normalPath
     * @param tumorPath
     * @param refGenome
     * @param outputPrefix
     * @param algorithm
     * @param mode
     * @param nThreads
     * @throws IOException
     */
    public void run(String normalPath, String tumorPath, String refGenome, String outputPrefix,
                    String algorithm, String mode, int nThreads) throws Exception {
        LOGGER.info("Beginning anonymization in " + mode + " mode");
        GlobalRandom.setSeed(randomSeed);
        AnonymizerAlgorithm anonymizer = getAnonymizer(algorithm, normalPath, tumorPath, refGenome, outputPrefix);
        if (tmpDir != null) anonymizer.setTmpDir(tmpDir);
        if(maxDepth > 100) anonymizer.setMaxDepth(maxDepth);
        anonymizer.setThreadNumber(nThreads);
        anonymizer.setMaxReadsInRam(maxReadsInMemory);
        anonymizer.setMinimumMappingQuality(minMappingQuality);
        anonymizer.setPartitions();
        long start1 = System.currentTimeMillis();
        anonymizer.queryReadsToExclude();
        long end1 = System.currentTimeMillis();
        LOGGER.info("Initial read query phase finished in: "+ (double) (end1-start1)/1000 + " seconds");
        long start2 = System.currentTimeMillis();
        anonymizer.anonymizeReads();
        long end2 = System.currentTimeMillis();
        LOGGER.info("Read Anonymization phase finished in: "+ (double) (end2-start2)/1000 + " seconds");
    }

    public void setTmpDir(File tmpDir) {
        this.tmpDir = tmpDir;
    }

    private void setMinimumMappingQuality(int minMQ) {
        this.minMappingQuality = minMQ;
    }

    public void setRandomSeed(int randomSeed) {
        this.randomSeed = randomSeed;
    }

    public void setMaxReadsInMemory(int maxReadsInMemory) {
        this.maxReadsInMemory = maxReadsInMemory;
    }

    private void setMaxDepth(int maxDepth) {this.maxDepth = maxDepth;}

    private static AnonymizerAlgorithm getAnonymizer(String algorithm, String normalPath, String tumorPath,
                                                     String refGenome, String outputPrefix) {
        AnonymizerAlgorithm anonymizer = null;
        if (AnonymizerAlgorithm.SHORT_READ_ALGORITHM.equals(algorithm)) anonymizer =
                new ShortReadAnonymizer(normalPath, tumorPath, refGenome, outputPrefix);
        return anonymizer;
    }

    public static void main(String[] args) {
        GenomeAnonymizer appInstance = new GenomeAnonymizer();
        long start2 = System.currentTimeMillis();
        Options options = buildCommandLineArguments();
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine commandLine = parser.parse(options, args);
            String appName = "GenomeAnonymizer" + " v" + VERSION;
            LOGGER.info("Running " + appName);
            if(args.length==0 || commandLine.hasOption("h")){
                HelpFormatter formatter = new HelpFormatter();
                String header = "OPTIONS";
                String footer = "";
                String cmdLineSyntax = "java -jar build/libs/GenomeAnonymizer-" + VERSION + ".jar";
                if(args.length==0) footer = "No arguments provided. Displaying default help message.";
                formatter.setOptionComparator(null);
                formatter.printHelp(cmdLineSyntax, header, options, footer, true);
                return;
            }
            String normalPath = commandLine.getOptionValue("in");
            String tumorPath = commandLine.getOptionValue("it");
            String refGenome = commandLine.getOptionValue("r");
            String outputPrefix = commandLine.getOptionValue("o", removeSuffixIfExists(normalPath, BAM_FILE));
            int nThreads = Integer.parseInt(commandLine.getOptionValue("t", "12"));
            int maxReadsInMemory = Integer.parseInt(commandLine.getOptionValue("maxReadsInMemory", String.valueOf(DEFAULT_MAX_READS_IN_RAM)));
            appInstance.setMaxReadsInMemory(maxReadsInMemory);
            String mode = commandLine.getOptionValue("m", DEFAULT_RUN_MODE_FUNCTIONALITY);
            int minMQ = Integer.parseInt(commandLine.getOptionValue("minMQ", String.valueOf(DEFAULT_MIN_MAPPING_QUALITY)));
            int randomSeed = Integer.parseInt(commandLine.getOptionValue("s", "-1"));
            appInstance.setMinimumMappingQuality(minMQ);
            appInstance.setRandomSeed(randomSeed);
            if (commandLine.hasOption("tmpDir")) {
                String tmpDirPath = commandLine.getOptionValue("tmpDir");
                appInstance.setTmpDir(new File(tmpDirPath));
            }
            if (commandLine.hasOption("maxDepth")) {
                int maxDepth = Integer.parseInt(commandLine.getOptionValue("maxDepth"));
                if (maxDepth <= DEFAULT_MIN_MAX_DEPTH) {
                    throw new IllegalArgumentException("Maximum depth must be higher than " + DEFAULT_MIN_MAX_DEPTH +
                            ". If you want to disable this feature, use -1 as value.");
                }
                appInstance.setMaxDepth(maxDepth);
            }
            LOGGER.info("Running with parameters - \n normalPath: " + normalPath + "\n tumorPath: " + tumorPath +
                        "\n refGenome: " + refGenome + "\n outputPrefix: " + outputPrefix +
                        "\n mode: " + mode + "\n minMQ: " + minMQ + "\n nThreads: " + nThreads);
            if(DEFAULT_RUN_MODE_FUNCTIONALITY.equals(mode)){
                if(commandLine.hasOption("v")) LOGGER.warning("Mode is not set to benchmark, but vcf file was provided," +
                        " default mode will be run normally," +
                        " but variants recorded in the vcf will not be kept");
                appInstance.run(normalPath, tumorPath, refGenome, outputPrefix, AnonymizerAlgorithm.SHORT_READ_ALGORITHM,
                        DEFAULT_RUN_MODE_FUNCTIONALITY, nThreads);
            }
        }
        catch (Exception e){
            LOGGER.log(Level.SEVERE,
                    "Fatal error happened: " + e.getMessage() + "\n halting execution and" +
                            " printing stacktrace",
                    e);
            System.exit(1);
        }
        long end2 = System.currentTimeMillis();
        LOGGER.info("Completed execution in: "+ (double) (end2-start2)/1000 + " seconds");
    }

    private static Options buildCommandLineArguments(){
        Options options = new Options();
        options.addOption("h", "Prints help message with information about the arguments");
        options.addOption(Option.builder("in")
                .desc("Input mappings coming from the normal sample (SAM/BAM/CRAM)")
                .argName("FILE")
                .hasArg(true)
                .build());
        options.addOption(Option.builder("it")
                .desc("Input mappings coming from the tumoral sample (SAM/BAM/CRAM)")
                .argName("FILE")
                .hasArg(true)
                .build());
        options.addOption(Option.builder("r")
                .desc("Reference genome against which the samples are aligned (.fasta, .fa)")
                .argName("FILE")
                .hasArg(true)
                .build());
        options.addOption(Option.builder("o")
                .desc("Prefix to name the output files")
                .argName("STRING")
                .hasArg(true)
                .build());
        options.addOption(Option.builder("t")
                .desc("Number of threads to run the anonymizer (default=12)")
                .hasArg(true)
                .argName("INTEGER")
                .type(Integer.class)
                .build());
        options.addOption(Option.builder("maxReadsInMemory")
                .desc("Limit of reads to be kept in memory. Trades memory consumption for IO effort, set lower for systems with low memory capacity "
                        + "(default=" + DEFAULT_MAX_READS_IN_RAM + ")")
                .hasArg(true)
                .argName("INTEGER")
                .type(Integer.class)
                .build());
        options.addOption(Option.builder("m")
                .desc("""
                        Mode that defines the functionality of the anonymizer, between: \
                        'default': Anonymizes all possible germline variants\
                        """)
                .argName("STRING")
                .hasArg(true)
                .build());
        options.addOption(Option.builder("s")
                .desc("Seed for random number generation. Use -1 for random seed. Default is -1")
                .argName("INTEGER")
                .hasArg(true)
                .build());
        options.addOption(Option.builder("minMQ")
                .desc("Minimum mapping quality for reads to be considered in the anonymization process (default=0; scale=0-60 PHRED)")
                .argName("INTEGER")
                .hasArg(true)
                .build());
        options.addOption(Option.builder("tmpDir")
                .desc("Temporary directory to store intermediate files")
                .argName("DIRECTORY")
                .hasArg(true)
                .build());
        options.addOption( Option.builder("maxDepth")
                .desc("Maximum depth of reads to be considered in the anonymization process. If set, reads in regions with depth " +
                        "higher than this value will be excluded (default=10000)" +
                        ". Use -1 to disable this feature" +
                        ". If set to a value lower than 100, it will be ignored.")
                .argName("INTEGER")
                .hasArg(true)
                .build());
        return options;
    }

    private static Logger logConfigure(){
        Logger answer = Logger.getLogger(GenomeAnonymizer.class.getName());
        answer.setLevel(Level.ALL);
        ConsoleHandler consoleHandler = new ConsoleHandler();
        consoleHandler.setFormatter(new SimpleFormatter());
        answer.addHandler(consoleHandler);
        answer.setUseParentHandlers(false);
        return answer;
    }

    private static String removeSuffixIfExists(String key, String suffix) {
        return key.endsWith(suffix)
                ? key.substring(0, key.length() - suffix.length())
                : key;
    }
}
