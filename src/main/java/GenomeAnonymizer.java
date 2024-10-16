import analysis.AnonymizerAlgorithm;
import analysis.ShortReadAnonymizer;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import io.SamplePairReadAlignmentReader;

import java.io.File;
import java.io.IOException;

public class GenomeAnonymizer {


    public void run(String normalPath, String tumorPath, String refGenome, String algorithm, int nThreads) throws IOException {
        AnonymizerAlgorithm anonymizer = null;
        if (AnonymizerAlgorithm.SHORT_READ_ALGORITHM.equals(algorithm)) anonymizer = new ShortReadAnonymizer();
        //Try with intervals as chr for parallelization, then maybe chr chunks too
        // try(SamplePairReadAlignmentReader pairPileupReader = new SamplePairReadAlignmentReader(normalPath, tumorPath, refGenome, intervals)
        try(SamplePairReadAlignmentReader pairPileupReader = new SamplePairReadAlignmentReader(normalPath, tumorPath, refGenome);
            IndexedFastaSequenceFile referenceWalker = new IndexedFastaSequenceFile(new File(refGenome));) {
            assert anonymizer != null: "No Anonymizer class impl was instantiated";
            anonymizer.anonymizeReads(pairPileupReader, referenceWalker);
        }
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
        appInstance.run(normalPath, tumorPath, refGenome, AnonymizerAlgorithm.SHORT_READ_ALGORITHM, 1);
        long end2 = System.currentTimeMillis();
        System.out.println("Elapsed Time in seconds: "+ (double) (end2-start2)/1000);
    }
}
