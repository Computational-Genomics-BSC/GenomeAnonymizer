package analysis;

import genomicelements.PairCalledVariation;
import genomicelements.GenomicRegion;
import genomicelements.GenomicRegionBaseImpl;
import genomicelements.Signal;
import utils.GlobalRandom;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.FastaSequenceIndexEntry;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.apache.commons.cli.*;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.logging.ConsoleHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;
import static io.VCFReader.readVCF;

/**
 * Main class that executes the anonymization method on sequencing data
 * @author Nicolas Gaitan
 */
public class GenomeAnonymizer {

    public static final String VERSION = "0.0.5";
    private static final Logger LOGGER = logConfigure();

    public static final String DEFAULT_RUN_MODE_FUNCTIONALITY = "default";
    public static final String SOMATIC_BENCHMARK_RUN_MODE_FUNCTIONALITY = "benchmark";

    public final static int DEFAULT_MIN_MAPPING_QUALITY = 1;
    public final static int DEFAULT_QUERY_REGION_LENGTH = 2000;
    
    public final static String BAM_FILE = ".bam";
    public final static String SAM_FILE = ".sam";
    public final static String CRAM_FILE = ".cram";
    private static final int MAX_QUERY_REGION_SIZE = 100000;

    //Optional arguments as attributes
    //String bedFile;
    private int minMappingQuality = DEFAULT_MIN_MAPPING_QUALITY;
    private int queryRegionLength = DEFAULT_QUERY_REGION_LENGTH;
    private String canvasN;
    private String canvasT;
    private Map<String, Map<Integer, PairCalledVariation>> somaticVariantsToKeep = new HashMap<>();
    private List<GenomicRegion> queryRegions = new ArrayList<>();
    private int randomSeed;

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
                    String algorithm, String mode, boolean merge, int nThreads) throws Exception {
        LOGGER.info("Beginning anonymization in " + mode + " mode");
        GlobalRandom.setSeed(randomSeed);
        AnonymizerAlgorithm anonymizer = getAnonymizer(algorithm, normalPath, tumorPath, refGenome, outputPrefix);
        anonymizer.setThreadNumber(nThreads);
        anonymizer.setMinimumMappingQuality(minMappingQuality);
        anonymizer.setQueryRegions(queryRegions);
        anonymizer.setVCFVariantsToKeep(somaticVariantsToKeep);
        anonymizer.setPartitions();
        long start1 = System.currentTimeMillis();
        anonymizer.queryReadsToExclude();
        long end1 = System.currentTimeMillis();
        LOGGER.info("Initial read query phase finished in: "+ (double) (end1-start1)/1000 + " seconds");
        long start2 = System.currentTimeMillis();
        if(merge){
            anonymizer.setCanvasFiles(this.canvasN, this.canvasT);
        }
        anonymizer.anonymizeReads();
        long end2 = System.currentTimeMillis();
        LOGGER.info("Read Anonymization phase finished in: "+ (double) (end2-start2)/1000 + " seconds");
//        anonymizer.mergeAnonymizedReads();
    }

    public void setCanvasN(String canvasN) {
        this.canvasN = canvasN;
    }

    public void setCanvasT(String canvasT) {
        this.canvasT = canvasT;
    }

    public void setSomaticVariantsToKeep(Map<String, Map<Integer, PairCalledVariation>> somaticVariantsToKeep) {
        this.somaticVariantsToKeep = somaticVariantsToKeep;
    }

    private void setMinimumMappingQuality(int minMQ) {
        this.minMappingQuality = minMQ;
    }

    public void setRandomSeed(int randomSeed) {
        this.randomSeed = randomSeed;
    }

    private static AnonymizerAlgorithm getAnonymizer(String algorithm, String normalPath, String tumorPath,
                                                     String refGenome, String outputPrefix) {
        AnonymizerAlgorithm anonymizer = null;
        if (AnonymizerAlgorithm.SHORT_READ_ALGORITHM.equals(algorithm)) anonymizer =
                new ShortReadAnonymizer(normalPath, tumorPath, refGenome, outputPrefix);
        return anonymizer;
    }

    /**
     * Computes query regions such that all intersecting reads will be included in the anonymized output.
     * @param mode Run mode to determine the source of the query regions
     * @pos queryRegions may be unsorted, they must be sorted later by sequence and coordinate
     */
    private void computeQueryRegions(String mode) {
        queryRegions = new ArrayList<>();
        int padding = queryRegionLength/2;
        if(SOMATIC_BENCHMARK_RUN_MODE_FUNCTIONALITY.equals(mode)){
            for (String sequenceName : somaticVariantsToKeep.keySet()){
                Map<Integer, PairCalledVariation> varsPerPos = somaticVariantsToKeep.get(sequenceName);
                for (PairCalledVariation var : varsPerPos.values()){
                    PairCalledVariation.VariantType variantType = var.getVariantType();
                    int start = var.getPos();
                    int end = var.getEnd();
                    int length = var.getLength();
                    String altEndsequenceName = sequenceName;
                    if (var.getBreakendRecord()!=null){
                        altEndsequenceName = var.getBreakendRecord().contig();
                        if(!sequenceName.equals(altEndsequenceName)){
                            end = var.getBreakendRecord().pos();
                        }
                    }
                    GenomicRegion queryRegion;
                    if (variantType == PairCalledVariation.VariantType.SNV){
                        queryRegion = new GenomicRegionBaseImpl(sequenceName, start - padding,
                                start + padding);
                    } else if (variantType == PairCalledVariation.VariantType.INS) {
                        queryRegion = new GenomicRegionBaseImpl(sequenceName, start  - padding,
                                start + length + padding);
                    } else if (variantType == PairCalledVariation.VariantType.TRA) {
                        queryRegion = new GenomicRegionBaseImpl(sequenceName, start - padding,
                                start + padding);
                        GenomicRegion secondQueryRegion = new GenomicRegionBaseImpl(altEndsequenceName, end - padding,
                                end + padding);
                        queryRegions.add(secondQueryRegion);
                    }
                    else if (variantType == PairCalledVariation.VariantType.INV){
                        if(start + padding > end - padding){
                            queryRegion = new GenomicRegionBaseImpl(sequenceName, start - padding,
                                    end + padding);
                        }
                        else{
                            queryRegion = new GenomicRegionBaseImpl(sequenceName, start - padding,
                                    start + padding);
                            GenomicRegion secondQueryRegion = new GenomicRegionBaseImpl(altEndsequenceName, end - padding,
                                    end + padding);
                            queryRegions.add(secondQueryRegion);
                        }
                    }
                    else{
                        if (length < MAX_QUERY_REGION_SIZE){
                            queryRegion = new GenomicRegionBaseImpl(sequenceName, start - padding,
                                    end + padding);
                        }
                        else{
                            queryRegion = new GenomicRegionBaseImpl(sequenceName, start - padding,
                                    start + padding);
                            GenomicRegion secondQueryRegion = new GenomicRegionBaseImpl(sequenceName, end - padding,
                                    end + padding);
                            queryRegions.add(secondQueryRegion);
                        }
                    }
                    queryRegions.add(queryRegion);
                }
            }
        }
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
            int nThreads = Integer.parseInt(commandLine.getOptionValue("t", "4"));
            String mode = commandLine.getOptionValue("m", DEFAULT_RUN_MODE_FUNCTIONALITY);
            int minMQ = Integer.parseInt(commandLine.getOptionValue("minMQ", "1"));
            int randomSeed = Integer.parseInt(commandLine.getOptionValue("s", "-1"));
            appInstance.setMinimumMappingQuality(minMQ);
            appInstance.setRandomSeed(randomSeed);
            boolean merge = false;
            if (SOMATIC_BENCHMARK_RUN_MODE_FUNCTIONALITY.equals(mode)){
                if(!commandLine.hasOption("v")) throw new ParseException("benchmark mode requires VCF file, " +
                        "but none was provided");
                String vcfFilePath = commandLine.getOptionValue("v");
                Map<String, Map<Integer, PairCalledVariation>> somaticVariantsToKeep = readVCF(vcfFilePath);
                appInstance.setSomaticVariantsToKeep(somaticVariantsToKeep);
                appInstance.computeQueryRegions(mode);
                if(commandLine.hasOption("merge")){
//                    if(!commandLine.hasOption("bed")) throw new ParseException("benchmark mode with merge requires BED file, " +
//                            "but none was provided");
                    if(!commandLine.hasOption("canvasN")) throw new ParseException("benchmark mode with merge requires" +
                            " canvas normal file, but none was provided");
                    if(!commandLine.hasOption("canvasT")) throw new ParseException("benchmark mode with merge requires" +
                            " canvas tumoral file, but none was provided");
                    merge = true;
//                    String bedFilePath = commandLine.getOptionValue("bed");
//                    appInstance.setBedFile(bedFilePath);
                    String canvasNFilePath = commandLine.getOptionValue("canvasN");
                    appInstance.setCanvasN(canvasNFilePath);
                    String canvasTFilePath = commandLine.getOptionValue("canvasT");
                    appInstance.setCanvasT(canvasTFilePath);
                }
                appInstance.run(normalPath, tumorPath, refGenome, outputPrefix, AnonymizerAlgorithm.SHORT_READ_ALGORITHM,
                        SOMATIC_BENCHMARK_RUN_MODE_FUNCTIONALITY, merge, nThreads);
            }
            else{
                if(commandLine.hasOption("v")) LOGGER.warning("Mode is not set to benchmark, but vcf file was provided," +
                        " default mode will be run normally," +
                        " but variants recorded in the vcf will not be kept");
                appInstance.run(normalPath, tumorPath, refGenome, outputPrefix, AnonymizerAlgorithm.SHORT_READ_ALGORITHM,
                        DEFAULT_RUN_MODE_FUNCTIONALITY, false, nThreads);
            }
        }
        catch (Exception e){
            //System.err.println(e.getMessage());
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
                //.required(true)
                .build());
        options.addOption(Option.builder("it")
                .desc("Input mappings coming from the tumoral sample (SAM/BAM/CRAM)")
                .argName("FILE")
                .hasArg(true)
                //.required(true)
                .build());
        options.addOption(Option.builder("r")
                .desc("Reference genome against which the samples are aligned (.fasta, .fa)")
                .argName("FILE")
                .hasArg(true)
                //.required(true)
                .build());
        options.addOption(Option.builder("o")
                .desc("Prefix to name the output files")
                .argName("STRING")
                .hasArg(true)
                //.required(true)
                .build());
        options.addOption(Option.builder("t")
                .desc("Number of threads to run the anonymizer (default=4)")
                .hasArg(true)
                .argName("INTEGER")
                //.required(false)
                .type(Integer.class)
                .build());
        options.addOption(Option.builder("m")
                .desc("""
                        Mode that defines the functionality of the anonymizer, between: \
                        
                        'default': Anonymizes all possible germline variants\
                        
                        'benchmark': Keeps marked variants for somatic calling benchmarking extracting sample reads
                         into public mappings (requires a VCF file, and a canvas BAM pair matching the coverage of the samples)""")
                        // + "'functional:' ")
                .argName("STRING")
                .hasArg(true)
                //.required(false)
                .build());
        options.addOption(Option.builder("v")
                .desc("VCF file containing the variants to be kept  (.VCF), required for running the Anonymizer in benchmark mode")
                .argName("FILE")
                .hasArg(true)
                //.required(false)
                .build());
        options.addOption(Option.builder("s")
                .desc("Seed for random number generation. Use -1 for random seed. Default is -1")
                .argName("INTEGER")
                .hasArg(true)
                //.required(false)
                .build());
        options.addOption(Option.builder("merge")
                .desc("Option to cover the empty mapped regions of the anonymized output with reads from a canvas dataset" +
                        "Caution: Only use this option with the benchmark mode" +
                        "Requires: -bed: BED file with the sample regions" +
                        " -canvasN && -canvasT: Read mapping datasets to cover the normal and tumoral outputs")
                .argName("OPTION")
                .hasArg(false)
                //.required(false)
                .build());
        options.addOption(Option.builder("minMQ")
                .desc("Minimum mapping quality for reads to be considered in the anonymization process (default=1; scale=0-60 PHRED)")
                .argName("INTEGER")
                .hasArg(true)
                //.required(false)
                .build());
//        options.addOption(Option.builder("bed")
//                .desc("BED file where the regions from the original sample are defined, " +
//                        "and therefore will not be covered with the canvas")
//                .argName("FILE")
//                .hasArg(true)
//                //.required(false)
//                .build());
        options.addOption(Option.builder("canvasN")
                .desc("Dataset with read mappings to cover the normal dataset")
                .argName("FILE")
                .hasArg(true)
                //.required(false)
                .build());
        options.addOption(Option.builder("canvasT")
                .desc("Dataset with read mappings to cover the tumoral dataset")
                .argName("FILE")
                .hasArg(true)
                //.required(false)
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
