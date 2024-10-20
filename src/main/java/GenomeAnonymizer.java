import analysis.AnonymizerAlgorithm;
import analysis.ShortReadAnonymizer;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import io.SamplePairReadAlignmentReader;

import java.io.File;
import java.io.IOException;

public class GenomeAnonymizer {

    public final static String BAM_FILE = ".bam";
    public final static String SAM_FILE = ".sam";
    public final static String CRAM_FILE = ".cram";

    public void run(String normalPath, String tumorPath, String refGenome, String outputPrefix, boolean compressed, String algorithm, int nThreads) throws IOException {
        AnonymizerAlgorithm anonymizer = getAnonymizer(algorithm);
        //anonymizer.setRemoveUnmapped(false);
        long start1 = System.currentTimeMillis();
        // TODO: Parallelize per chromosome, and then, per reads to anonymize
        anonymizer.callVariation(normalPath, tumorPath, refGenome);
        long end1 = System.currentTimeMillis();
        System.out.println("Elapsed Time in seconds for variation calling: "+ (double) (end1-start1)/1000);
        long start2 = System.currentTimeMillis();
        anonymizer.anonymizeReads(normalPath, tumorPath, refGenome, outputPrefix, false);
        long end2 = System.currentTimeMillis();
        System.out.println("Elapsed Time in seconds for anonymization: "+ (double) (end2-start2)/1000);
        //long start3 = System.currentTimeMillis();
        //anonymizer.writeUnmodifiedReads(normalPath, tumorPath, outputPrefix, false);
        //long end3 = System.currentTimeMillis();
        //System.out.println("Elapsed Time in seconds for retrieving missing reads: "+ (double) (end2-start2)/1000);
    }

    private static AnonymizerAlgorithm getAnonymizer(String algorithm) {
        AnonymizerAlgorithm anonymizer = null;
        if (AnonymizerAlgorithm.SHORT_READ_ALGORITHM.equals(algorithm)) anonymizer = new ShortReadAnonymizer();
        //Try with intervals as chr for parallelization, then maybe chr chunks too
        // try(SamplePairReadAlignmentReader pairPileupReader = new SamplePairReadAlignmentReader(normalPath, tumorPath, refGenome, intervals)
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
//        try(SamplePairReadAlignmentReader pairPileupReader = new SamplePairReadAlignmentReader(normalPath, tumorPath, refGenome);) {
//            pairPileupReader.forEach(x -> System.out.println("Pileup seq=" + x.getRefenceSequenceName() + " pos=" + x.getReferencePos()
//                                        + " n_reads="+x.totalSize()));
//        }
        long start2 = System.currentTimeMillis();
        GenomeAnonymizer appInstance = new GenomeAnonymizer();
        appInstance.run(normalPath, tumorPath, refGenome, removeSuffixIfExists(normalPath, BAM_FILE), false, AnonymizerAlgorithm.SHORT_READ_ALGORITHM, 1);
        long end2 = System.currentTimeMillis();
        System.out.println("Total execution Time in seconds: "+ (double) (end2-start2)/1000);
    }

    private static String removeSuffixIfExists(String key, String suffix) {
        return key.endsWith(suffix)
                ? key.substring(0, key.length() - suffix.length())
                : key;
    }
}
