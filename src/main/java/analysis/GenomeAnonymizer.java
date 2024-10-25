package analysis;

import genomicelements.CalledVariation;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.FastaSequenceIndexEntry;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.SimpleFeature;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class GenomeAnonymizer {

    public static final String DEFAULT_RUN_MODE_FUNCTIONALITY = "DEFAULT";
    public static final String SOMATIC_BENCHMARK_RUN_MODE_FUNCTIONALITY = "SOMATIC_BENCHMARK";

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
                    String algorithm, String mode, String vcfFile, int nThreads) throws IOException {
        //anonymizer.setRemoveUnmapped(false);
        long start1 = System.currentTimeMillis();
        // TODO: Parallelize per chromosome, and then, per reads to anonymize
        Map<String, Map<Integer, List<CalledVariation>>> readGermlinesToAnonymize = new HashMap<>();
        if (nThreads==1){
            VariationClassifier classifier = new VariationClassifier();
            classifier.callVariation(normalPath, tumorPath, refGenome, mode, vcfFile);
            readGermlinesToAnonymize = classifier.getPotentialGermlinesPerRead();
        }
        else{
            List<SimpleFeature> partitions = getPartitions(refGenome, nThreads);
            readGermlinesToAnonymize = callVariationInParallel(normalPath, tumorPath, refGenome, mode, vcfFile, partitions, nThreads);
        }
        long end1 = System.currentTimeMillis();
        System.out.println("Elapsed Time in seconds for variation calling: "+ (double) (end1-start1)/1000);
        long start2 = System.currentTimeMillis();
        AnonymizerAlgorithm anonymizer = getAnonymizer(algorithm, readGermlinesToAnonymize);
        anonymizer.anonymizeReads(normalPath, tumorPath, refGenome, outputPrefix, false);
        long end2 = System.currentTimeMillis();
        System.out.println("Elapsed Time in seconds for anonymization: "+ (double) (end2-start2)/1000);
    }

    private Map<String, Map<Integer, List<CalledVariation>>> callVariationInParallel(String normalPath, String tumorPath, String refGenome,
                                                                                     String mode, String vcfFile, List<SimpleFeature> partitions, int nThreads) {
        assert (nThreads==partitions.size()): "Partition number is not equal to the number of threads";
        Map<String, Map<Integer, List<CalledVariation>>> readGermlinesToAnonymize = new HashMap<>();
        MultithreadClassifier[] mClassifiers = new MultithreadClassifier[partitions.size()];
        for(int i = 0; i < partitions.size(); i++){
            SimpleFeature region = partitions.get(i);
            mClassifiers[i] = new MultithreadClassifier(normalPath, tumorPath, refGenome, mode, vcfFile, region);
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
            String contig = classifier.getContig();
            Map<Integer, List<CalledVariation>> germlinesInPartition = classifier.getAnswer();
            Map<Integer, List<CalledVariation>> germlinesInContig = readGermlinesToAnonymize.computeIfAbsent(contig, v -> new HashMap<>());
            germlinesInContig.putAll(germlinesInPartition);
        }
        return readGermlinesToAnonymize;
    }

    private List<SimpleFeature> getPartitions(String refGenome, int nThreads) throws IOException {
        IndexedFastaSequenceFile reference = new IndexedFastaSequenceFile(new File(refGenome));
        //SAMSequenceDictionary seqDict = reference.getSequenceDictionary();
        //System.out.println(seqDict.isEmpty());
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
        if(nThreads>sequences.size()) {
            Collections.sort(sequences, Comparator.comparing(FastaSequenceIndexEntry::getSize));
            Collections.reverse(sequences);
            //int nNewPartitions = nThreads-currentPartitions;
            int availableThreads = nThreads-sequences.size();
            int[] partitionsPerSequence = new int[sequences.size()];
            Arrays.fill(partitionsPerSequence, 1);
            long basesPerThread = genomeSize/nThreads;
            // Estimate threads to be assigned to each contig
            for(int i = 0; i < sequences.size(); i++){
                partitionsPerSequence[i] += (int) (sequences.get(i).getSize() / basesPerThread);
                availableThreads -= partitionsPerSequence[i];
                assert (availableThreads>=0): "Available threads are lower than 0, this should not happen";
                if(availableThreads==0) break;
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
                    int currentLast = j == partitionsPerSequence[i]-1 ? (int) currentContig.getSize() : currentFirst + 1 + partitionSize;
                    SimpleFeature region = new SimpleFeature(contig, currentFirst, currentLast);
                    regions.add(region);
                    currentFirst += partitionSize;
                }
            }
        }
        else{
            regions = sequences.stream()
                    .map(sequenceIndexEntry -> new SimpleFeature(sequenceIndexEntry.getContig(), 1, (int) sequenceIndexEntry.getSize()))
                    .toList();
        }
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
                    String algorithm, int nThreads) throws IOException {
        run(normalPath, tumorPath, refGenome, outputPrefix, compressed, algorithm, DEFAULT_RUN_MODE_FUNCTIONALITY, "", nThreads);
    }

    private static AnonymizerAlgorithm getAnonymizer(String algorithm, Map<String, Map<Integer, List<CalledVariation>>> readGermlinesToAnonymize) {
        AnonymizerAlgorithm anonymizer = null;
        if (AnonymizerAlgorithm.SHORT_READ_ALGORITHM.equals(algorithm)) anonymizer = new ShortReadAnonymizer(readGermlinesToAnonymize);
        assert anonymizer != null: "No Anonymizer class impl was instantiated";
        return anonymizer;
    }

    public static void main(String[] args) throws IOException {
        // Test functionalities temporarily
        String normalPath = args[0];
        System.out.println(normalPath);
        String tumorPath = args[1];
        System.out.println(tumorPath);
        String refGenome = args[2];
        System.out.println(refGenome);
        String vcfFilePath = args[3];
        System.out.println(vcfFilePath);
        int nThreads = Integer.parseInt(args[4]);
        long start2 = System.currentTimeMillis();
        GenomeAnonymizer appInstance = new GenomeAnonymizer();
        appInstance.run(normalPath, tumorPath, refGenome, removeSuffixIfExists(normalPath, BAM_FILE), false,
                AnonymizerAlgorithm.SHORT_READ_ALGORITHM, SOMATIC_BENCHMARK_RUN_MODE_FUNCTIONALITY, vcfFilePath, nThreads);
        long end2 = System.currentTimeMillis();
        System.out.println("Total execution Time in seconds: "+ (double) (end2-start2)/1000);
    }

    private static String removeSuffixIfExists(String key, String suffix) {
        return key.endsWith(suffix)
                ? key.substring(0, key.length() - suffix.length())
                : key;
    }
}
