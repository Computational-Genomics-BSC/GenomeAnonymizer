package analysis;

import genomicelements.CalledVariation;
import genomicelements.GenomicRegion;
import genomicelements.GenomicRegionBaseImpl;
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

    public static final String VERSION = "0.0.3";
    private static final Logger LOGGER = logConfigure();


    public static final String DEFAULT_RUN_MODE_FUNCTIONALITY = "default";
    public static final String SOMATIC_BENCHMARK_RUN_MODE_FUNCTIONALITY = "benchmark";

    public final static String BAM_FILE = ".bam";
    public final static String SAM_FILE = ".sam";
    public final static String CRAM_FILE = ".cram";

    //Optional arguments as attributes
    String vcfFile;
    String bedFile;
    String canvasN;
    String canvasT;


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
        AnonymizerAlgorithm anonymizer = getAnonymizer(algorithm);
        List<GenomicRegion> partitions = getPartitions(refGenome, nThreads);
        anonymizer.setPartitions(partitions);
        anonymizer.queryReadsToExclude(normalPath, tumorPath, nThreads);
        Set<String> readsToExclude = anonymizer.getReadsToExclude();
        long start1 = System.currentTimeMillis();
        Map<String, List<CalledVariation>> readGermlinesToAnonymize =
                callVariationInParallel(normalPath, tumorPath, refGenome, mode, vcfFile, partitions, readsToExclude, nThreads);
        long end1 = System.currentTimeMillis();
        LOGGER.info("Variation calling phase finished in: "+ (double) (end1-start1)/1000 + " seconds");
        long start2 = System.currentTimeMillis();
        anonymizer.setReadGermlinesToAnonymize(readGermlinesToAnonymize);
        if(merge){
            anonymizer.setRegions(this.bedFile);
            anonymizer.setCanvasFiles(this.canvasN, this.canvasT);
        }
        anonymizer.anonymizeReads(normalPath, tumorPath, refGenome, outputPrefix, false);
        long end2 = System.currentTimeMillis();
        LOGGER.info("Anonymization phase finished in: "+ (double) (end2-start2)/1000 + " seconds");
    }

    private Map<String, List<CalledVariation>> callVariationInParallel(String normalPath, String tumorPath, String refGenome,
                                                                                     String mode, String vcfFile, List<GenomicRegion> partitions,
                                                                                     Set<String> readsToExclude, int nThreads) throws IOException {
        Map<String, List<CalledVariation>> readGermlinesToAnonymize = new HashMap<>();
        MultithreadClassifier[] mClassifiers = new MultithreadClassifier[partitions.size()];
        Map<String, Map<Integer,CalledVariation>> somaticVariantsToKeep = new HashMap<>();
        if(SOMATIC_BENCHMARK_RUN_MODE_FUNCTIONALITY.equals(mode)) somaticVariantsToKeep = readVCF(vcfFile);
        for(int i = 0; i < partitions.size(); i++){
            GenomicRegion region = partitions.get(i);
            mClassifiers[i] = new MultithreadClassifier(normalPath, tumorPath, refGenome, mode, vcfFile, region, readsToExclude);
            mClassifiers[i].setVCFVariantsToKeep(somaticVariantsToKeep);
        }
        ExecutorService executorService = Executors.newFixedThreadPool(nThreads);
        List<CompletableFuture<Void>> futures = new ArrayList<>();
        for(int i = 0; i < mClassifiers.length; i++){
            MultithreadClassifier runnableClassifier = mClassifiers[i];
            CompletableFuture<Void> future = CompletableFuture.runAsync(runnableClassifier, executorService);
            futures.add(future);
        }
        CompletableFuture.allOf(futures.toArray(new CompletableFuture[0])).join();
        executorService.shutdown();
        for(MultithreadClassifier classifier : mClassifiers){
            Map<String, List<CalledVariation>> germlinesInPartitionReads = classifier.getAnswer();
            for(var entryByReadAlnId : germlinesInPartitionReads.entrySet()){
                String readAlnId = entryByReadAlnId.getKey();
                List<CalledVariation> newVariationsInReadAln = entryByReadAlnId.getValue();
                if(readGermlinesToAnonymize.containsKey(readAlnId)){
                    List<CalledVariation> currentVariationsInReadAln = readGermlinesToAnonymize.get(readAlnId);
                    currentVariationsInReadAln.addAll(newVariationsInReadAln);
                }
                else{
                    readGermlinesToAnonymize.put(readAlnId, newVariationsInReadAln);
                }
            }
        }
        return readGermlinesToAnonymize;
    }

    private List<GenomicRegion> getPartitions(String refGenome, int nThreads) throws IOException {
        IndexedFastaSequenceFile reference = new IndexedFastaSequenceFile(new File(refGenome));
        FastaSequenceIndex refIndexes = reference.getIndex();
        long genomeSize = 0;
        List<FastaSequenceIndexEntry> sequences = new ArrayList<>();
        Map<String, Integer> refSequenceOrder = new HashMap<>();
        int idx = 0;
        for (FastaSequenceIndexEntry refEntry : refIndexes){
            sequences.add(refEntry);
            genomeSize += refEntry.getSize();
            refSequenceOrder.put(refEntry.getContig(), idx);
            idx++;
        }
        reference.close();
        // long genomeSize = seqDict.getReferenceLength();
        List<GenomicRegion> regions = new ArrayList<>();
        Collections.sort(sequences, Comparator.comparing(FastaSequenceIndexEntry::getSize));
        //long basesPerThread = sequences.size() % 2 == 0 ? (sequences.get(sequences.size()/2).getSize() + sequences.get((sequences.size()/2)-1).getSize())/2 :
        //        sequences.get(sequences.size()/2).getSize();
        Collections.reverse(sequences);
        //int nNewPartitions = nThreads-currentPartitions;
        //int availableThreads = nThreads-sequences.size();
        //System.out.println("init avail=" + availableThreads);
        int[] partitionsPerSequence = new int[sequences.size()];
        Arrays.fill(partitionsPerSequence, 1);
        long basesPerThread = genomeSize/(nThreads* 10L);
        //long basesPerThread = availableThreads < sequences.size() ? sequences.get(availableThreads).getSize() :
        //        sequences.get(sequences.size() - 1).getSize();
        // Estimate threads to be assigned to each contig
        for(int i = 0; i < sequences.size(); i++){
            //if(availableThreads <= 0) break;
            partitionsPerSequence[i] += (int) (sequences.get(i).getSize() / basesPerThread);
            //availableThreads -= partitionsPerSequence[i];
        }
        // Make propper contig partitions into the regions, based on n assigned threads
        for(int i = 0; i < sequences.size(); i++){
            FastaSequenceIndexEntry currentContig = sequences.get(i);
            String contig = currentContig.getContig();
            int contigIdx = refSequenceOrder.get(contig);
            int contigLength = (int) currentContig.getSize();
            int partitionSize = contigLength / partitionsPerSequence[i];
            int currentFirst = 1;
            for(int j = 0; j < partitionsPerSequence[i]; j++){
                //Be careful with very large chromosomes, with humans there should not be a problem
                int currentLast = j == partitionsPerSequence[i]-1 ? (int) currentContig.getSize() : currentFirst + partitionSize;
                GenomicRegion region = new GenomicRegionBaseImpl(contig, contigIdx, currentFirst, currentLast);
                regions.add(region);
                currentFirst += partitionSize + 1;
            }
        }
        //DEBUG
        //regions.forEach(r-> System.out.println("# " + r.getContig() + " " + r.getStart() + " " + r.getEnd()));
        //System.exit(0);
        //DEBUG
        return regions;
    }

    public String getVcfFile() {
        return vcfFile;
    }

    public void setVcfFile(String vcfFile) {
        this.vcfFile = vcfFile;
    }

    public String getBedFile() {
        return bedFile;
    }

    public void setBedFile(String bedFile) {
        this.bedFile = bedFile;
    }

    public String getCanvasN() {
        return canvasN;
    }

    public void setCanvasN(String canvasN) {
        this.canvasN = canvasN;
    }

    public String getCanvasT() {
        return canvasT;
    }

    public void setCanvasT(String canvasT) {
        this.canvasT = canvasT;
    }

//    /**
//     * Run anonymization in default mode
//     * @param normalPath
//     * @param tumorPath
//     * @param refGenome
//     * @param outputPrefix
//     * @param compressed
//     * @param algorithm
//     * @param nThreads
//     * @throws IOException
//     */
//    public void run(String normalPath, String tumorPath, String refGenome, String outputPrefix, boolean compressed,
//                    String algorithm, int nThreads) throws Exception {
//        run(normalPath, tumorPath, refGenome, outputPrefix, algorithm, DEFAULT_RUN_MODE_FUNCTIONALITY, nThreads);
//    }

    private static AnonymizerAlgorithm getAnonymizer(String algorithm) {
        AnonymizerAlgorithm anonymizer = null;
        if (AnonymizerAlgorithm.SHORT_READ_ALGORITHM.equals(algorithm)) anonymizer = new ShortReadAnonymizer();
        assert anonymizer != null: "No Anonymizer class impl was instantiated";
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
                //System.exit(0);
            }
            String normalPath = commandLine.getOptionValue("in");
            String tumorPath = commandLine.getOptionValue("it");
            String refGenome = commandLine.getOptionValue("r");
            String outputPrefix = commandLine.getOptionValue("o", removeSuffixIfExists(normalPath, BAM_FILE));
            int nThreads = Integer.parseInt(commandLine.getOptionValue("t", "4"));
            String mode = commandLine.getOptionValue("m", DEFAULT_RUN_MODE_FUNCTIONALITY);
            boolean merge = false;
            if (SOMATIC_BENCHMARK_RUN_MODE_FUNCTIONALITY.equals(mode)){
                if(!commandLine.hasOption("v")) throw new ParseException("benchmark mode requires VCF file, " +
                        "but none was provided");
                String vcfFilePath = commandLine.getOptionValue("v");
                appInstance.setVcfFile(vcfFilePath);
                if(commandLine.hasOption("merge")){
                    if(!commandLine.hasOption("bed")) throw new ParseException("benchmark mode with merge requires BED file, " +
                            "but none was provided");
                    if(!commandLine.hasOption("canvasN")) throw new ParseException("benchmark mode with merge requires" +
                            " canvas normal file, but none was provided");
                    if(!commandLine.hasOption("canvasT")) throw new ParseException("benchmark mode with merge requires" +
                            " canvas tumoral file, but none was provided");
                    merge = true;
                    String bedFilePath = commandLine.getOptionValue("bed");
                    appInstance.setBedFile(bedFilePath);
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
            LOGGER.severe("Exception happened during execution: " + e.getMessage() + "\n halting execution and" +
                    " printing stacktrace");
            e.printStackTrace();
            System.exit(1);
        }
        long end2 = System.currentTimeMillis();
        LOGGER.info("Completed execution in: "+ (double) (end2-start2)/1000 + " seconds");
//        catch (ParseException e){
//            System.err.println("Invalid arguments: " + e.getMessage());
//        }
//        catch (IOException e){
//            System.err.println("Invalid arguments: " + e.getMessage());
//        }
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
                .desc("Number of threads to run the anonymizer (default=4)")// + "'functional:' ")
                .hasArg(true)
                .argName("INT")
                //.required(false)
                .type(Integer.class)
                .build());
        options.addOption(Option.builder("m")
                .desc("""
                        Mode that defines the functionality of the anonymizer, between: \
                        
                        'default': Anonymizes all possible germline variants\
                        
                        'benchmark': Keeps marked variants for somatic calling benchmarking (requires a VCF file)""")// + "'functional:' ")
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
        options.addOption(Option.builder("merge")
                .desc("Option to cover the empty mapped regions of the anonymized output with reads from a canvas dataset" +
                        "Caution: Only use this option with the benchmark mode" +
                        "Requires: -bed: BED file with the sample regions" +
                        " -canvasN && -canvasT: Read mapping datasets to cover the normal and tumoral outputs")
                .argName("OPTION")
                .hasArg(false)
                //.required(false)
                .build());
        options.addOption(Option.builder("bed")
                .desc("BED file where the regions from the original sample are defined, " +
                        "and therefore will not be covered with the canvas")
                .argName("FILE")
                .hasArg(true)
                //.required(false)
                .build());
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
        //consoleHandler.setOutputStream(System.out);
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
