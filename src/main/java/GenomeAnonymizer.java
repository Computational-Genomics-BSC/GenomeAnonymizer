import io.SamplePairReadAlignmentReader;
import java.io.IOException;

public class GenomeAnonymizer {

    public static void main(String[] args) throws IOException {
        // Test functionalities temporarily
        String normalPath = args[0];
        System.out.println(normalPath);
        String tumorPath = args[1];
        System.out.println(tumorPath);
        String refGenome = args[2];
        System.out.println(refGenome);
        try(SamplePairReadAlignmentReader pairPileupReader = new SamplePairReadAlignmentReader(normalPath, tumorPath, refGenome);) {
            pairPileupReader.forEach(x -> System.out.println("Pileup seq=" + x.getRefenceSequenceName() + " pos=" + x.getReferencePos()
                                        + " n_reads="+x.totalSize()));
        }
    }
}
