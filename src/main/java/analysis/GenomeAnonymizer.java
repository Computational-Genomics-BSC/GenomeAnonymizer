package analysis;

import genomicelements.CalledVariation;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.FastaSequenceIndexEntry;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.SimpleFeature;
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

/**
 * Main class that executes the anonymization method on sequencing data
 * @author Nicolas Gaitan
 */
public class GenomeAnonymizer {

    public static final String VERSION = "0.0.2";
    private static final Logger LOGGER = logConfigure();


    public static final String DEFAULT_RUN_MODE_FUNCTIONALITY = "default";
    public static final String SOMATIC_BENCHMARK_RUN_MODE_FUNCTIONALITY = "benchmark";

    public final static String BAM_FILE = ".bam";
    public final static String SAM_FILE = ".sam";
    public final static String CRAM_FILE = ".cram";

    /**
     * Run anonymizer with the benchmark of somatic variants functionality. Any somatic variant will be sparred from anonymization
     * @param normalPath
     * @param tumorPath
     * @param refGenome
     * @param outputPrefix
     * @param compressed
     * @param algorithm
     * @param mode
     * @param vcfFile
     * @param nThreads
     * @throws IOException
     */
    public void run(String normalPath, String tumorPath, String refGenome, String outputPrefix, boolean compressed,
                    String algorithm, String mode, String vcfFile, int nThreads) throws Exception {
        LOGGER.info("Beginning anonymization in " + mode + " mode");
        AnonymizerAlgorithm anonymizer = getAnonymizer(algorithm);
        List<SimpleFeature> partitions = getPartitions(refGenome, nThreads);
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
        anonymizer.anonymizeReads(normalPath, tumorPath, refGenome, outputPrefix, false);
        long end2 = System.currentTimeMillis();
        LOGGER.info("Anonymization phase finished in: "+ (double) (end2-start2)/1000 + " seconds");
    }

    private Map<String, List<CalledVariation>> callVariationInParallel(String normalPath, String tumorPath, String refGenome,
                                                                                     String mode, String vcfFile, List<SimpleFeature> partitions,
                                                                                     Set<String> readsToExclude, int nThreads) {
        Map<String, List<CalledVariation>> readGermlinesToAnonymize = new HashMap<>();
        MultithreadClassifier[] mClassifiers = new MultithreadClassifier[partitions.size()];
        for(int i = 0; i < partitions.size(); i++){
            SimpleFeature region = partitions.get(i);
            mClassifiers[i] = new MultithreadClassifier(normalPath, tumorPath, refGenome, mode, vcfFile, region, readsToExclude);
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

    private List<SimpleFeature> getPartitions(String refGenome, int nThreads) throws IOException {
        IndexedFastaSequenceFile reference = new IndexedFastaSequenceFile(new File(refGenome));
        FastaSequenceIndex refIndexes = reference.getIndex();
        long genomeSize = 0;
        List<FastaSequenceIndexEntry> sequences = new ArrayList<>();
        for (FastaSequenceIndexEntry refEntry : refIndexes){
            sequences.add(refEntry);
            genomeSize += refEntry.getSize();
        }
        reference.close();
        // long genomeSize = seqDict.getReferenceLength();
        List<SimpleFeature> regions = new ArrayList<>();
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
            int contigLength = (int) currentContig.getSize();
            int partitionSize = contigLength / partitionsPerSequence[i];
            int currentFirst = 1;
            for(int j = 0; j < partitionsPerSequence[i]; j++){
                //Be careful with very large chromosomes, with humans there should not be a problem
                int currentLast = j == partitionsPerSequence[i]-1 ? (int) currentContig.getSize() : currentFirst + partitionSize;
                SimpleFeature region = new SimpleFeature(contig, currentFirst, currentLast);
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

    /**
     * Run anonymization in default mode
     * @param normalPath
     * @param tumorPath
     * @param refGenome
     * @param outputPrefix
     * @param compressed
     * @param algorithm
     * @param nThreads
     * @throws IOException
     */
    public void run(String normalPath, String tumorPath, String refGenome, String outputPrefix, boolean compressed,
                    String algorithm, int nThreads) throws Exception {
        run(normalPath, tumorPath, refGenome, outputPrefix, compressed, algorithm, DEFAULT_RUN_MODE_FUNCTIONALITY, "", nThreads);
    }

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
                // TODO: Change when it is set to be run as a jar, or container
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
            String vcfFilePath;
            if (SOMATIC_BENCHMARK_RUN_MODE_FUNCTIONALITY.equals(mode)){
                if(!commandLine.hasOption("v")) throw new ParseException("benchmark mode requires VCF file, but none was provided");
                vcfFilePath = commandLine.getOptionValue("v");
                appInstance.run(normalPath, tumorPath, refGenome, outputPrefix, false,
                        AnonymizerAlgorithm.SHORT_READ_ALGORITHM, SOMATIC_BENCHMARK_RUN_MODE_FUNCTIONALITY, vcfFilePath, nThreads);
            }
            else{
                if(commandLine.hasOption("v")) LOGGER.warning("Mode is not set to benchmark, but vcf file was provided," +
                        " default mode will be run normally," +
                        " but variants recorded in the vcf will not be kept");
                appInstance.run(normalPath, tumorPath, refGenome, outputPrefix, false,
                        AnonymizerAlgorithm.SHORT_READ_ALGORITHM, nThreads);
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
