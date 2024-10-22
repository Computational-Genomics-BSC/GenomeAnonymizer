import analysis.AnonymizerAlgorithm;
import analysis.ShortReadAnonymizer;
import analysis.VariationClassifier;
import genomicelements.CalledVariation;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import io.SamplePairReadAlignmentReader;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static analysis.AnonymizerAlgorithm.DEFAULT_RUN_MODE_FUNCTIONALITY;
import static analysis.AnonymizerAlgorithm.SOMATIC_BENCHMARK_RUN_MODE_FUNCTIONALITY;

public class GenomeAnonymizer {

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
        VariationClassifier classifier = new VariationClassifier();
        Map<String, Map<Integer, List<CalledVariation>>> readGermlinesToAnonymize = classifier
                .callVariation(normalPath, tumorPath, refGenome, SOMATIC_BENCHMARK_RUN_MODE_FUNCTIONALITY, vcfFile);
        long end1 = System.currentTimeMillis();
        System.out.println("Elapsed Time in seconds for variation calling: "+ (double) (end1-start1)/1000);
        long start2 = System.currentTimeMillis();
        AnonymizerAlgorithm anonymizer = getAnonymizer(algorithm, readGermlinesToAnonymize);
        anonymizer.anonymizeReads(normalPath, tumorPath, refGenome, outputPrefix, false);
        long end2 = System.currentTimeMillis();
        System.out.println("Elapsed Time in seconds for anonymization: "+ (double) (end2-start2)/1000);
        //long start3 = System.currentTimeMillis();
        //anonymizer.writeUnmodifiedReads(normalPath, tumorPath, outputPrefix, false);
        //long end3 = System.currentTimeMillis();
        //System.out.println("Elapsed Time in seconds for retrieving missing reads: "+ (double) (end2-start2)/1000);
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
//        try(SamplePairReadAlignmentReader pairPileupReader = new SamplePairReadAlignmentReader(normalPath, tumorPath, refGenome);) {
//            pairPileupReader.forEach(x -> System.out.println("Pileup seq=" + x.getRefenceSequenceName() + " pos=" + x.getReferencePos()
//                                        + " n_reads="+x.totalSize()));
//        }
        long start2 = System.currentTimeMillis();
        GenomeAnonymizer appInstance = new GenomeAnonymizer();
        //appInstance.run(normalPath, tumorPath, refGenome, removeSuffixIfExists(normalPath, BAM_FILE), false, AnonymizerAlgorithm.SHORT_READ_ALGORITHM, 1);
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
