package analysis;

import genomicelements.CalledVariation;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.SimpleFeature;

import java.io.File;
import java.io.IOException;
import java.util.*;

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
            classifier.callVariation(normalPath, tumorPath, refGenome, SOMATIC_BENCHMARK_RUN_MODE_FUNCTIONALITY, vcfFile);
            readGermlinesToAnonymize = classifier.getPotentialGermlinesPerRead();
        }
        else{
            IndexedFastaSequenceFile reference = new IndexedFastaSequenceFile(new File(refGenome));
            SAMSequenceDictionary seqDict = reference.getSequenceDictionary();
            reference.close();
            long genomeSize = seqDict.getReferenceLength();
            List<SAMSequenceRecord> sequences = seqDict.getSequences();
            List<SimpleFeature> regions = new ArrayList<>();
            if(nThreads>sequences.size()) {
                Collections.sort(sequences, Comparator.comparing(SAMSequenceRecord::getSequenceLength));
                Collections.reverse(sequences);
                //int nNewPartitions = nThreads-currentPartitions;
                int availableThreads = nThreads-sequences.size();
                int[] partitionsPerSequence = new int[sequences.size()];
                Arrays.fill(partitionsPerSequence, 1);
                long basesPerThread = genomeSize/nThreads;
                for(int i = 0; i < sequences.size(); i++){
                    partitionsPerSequence[i] += (int) (sequences.get(i).getSequenceLength() / basesPerThread);
                    availableThreads -= partitionsPerSequence[i];
                    assert (availableThreads>=0): "Available threads are lower than 0, this shouldnt happen";
                    if(availableThreads==0) break;
                }
            }
            else{
                regions = sequences.stream()
                        .map(samSeq -> new SimpleFeature(samSeq.getContig(), samSeq.getStart(), samSeq.getEnd()))
                        .toList();
            }
        }
        long end1 = System.currentTimeMillis();
        System.out.println("Elapsed Time in seconds for variation calling: "+ (double) (end1-start1)/1000);
        long start2 = System.currentTimeMillis();
        AnonymizerAlgorithm anonymizer = getAnonymizer(algorithm, readGermlinesToAnonymize);
        anonymizer.anonymizeReads(normalPath, tumorPath, refGenome, outputPrefix, false);
        long end2 = System.currentTimeMillis();
        System.out.println("Elapsed Time in seconds for anonymization: "+ (double) (end2-start2)/1000);
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
        long start2 = System.currentTimeMillis();
        GenomeAnonymizer appInstance = new GenomeAnonymizer();
        appInstance.run(normalPath, tumorPath, refGenome, removeSuffixIfExists(normalPath, BAM_FILE), false,
                AnonymizerAlgorithm.SHORT_READ_ALGORITHM, SOMATIC_BENCHMARK_RUN_MODE_FUNCTIONALITY, vcfFilePath, 1);
        long end2 = System.currentTimeMillis();
        System.out.println("Total execution Time in seconds: "+ (double) (end2-start2)/1000);
    }

    private static String removeSuffixIfExists(String key, String suffix) {
        return key.endsWith(suffix)
                ? key.substring(0, key.length() - suffix.length())
                : key;
    }
}
