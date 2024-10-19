package analysis;

import genomicelements.*;
import genomicelements.CalledVariation.SomaticVariationType;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import io.SamplePairReadAlignmentReader;

import java.io.File;
import java.io.IOException;
import java.util.*;

import static genomicelements.ShortAnonymizedReadPair.*;

public class ShortReadAnonymizer implements AnonymizerAlgorithm{

    public static final int NORMAL_DATASET_IDX = 0;
    public static final int TUMORAL_DATASET_IDX = 1;
    public static final int OUTPUT_FILE_NUMBER = 4;

    public static final int SLIDING_WINDOW_LIMIT = 200;

    Iterable<PairedPileup> reader;
    // Map containing all potential germlines (value: List), per pair (nested key, 0 or 1), per read (key)
    Map<String, Map<Integer, List<CalledVariation>>> potentialGermlinesPerAnonRead;
    // Map collecting to-be-anonymized reads while they can be masked and written
    Map<String, ShortAnonymizedReadPair> anonReadContainer;
    Set<String> anonReadNames;

    public ShortReadAnonymizer() {
        anonReadContainer = new HashMap<>();
        potentialGermlinesPerAnonRead = new HashMap<>();
        anonReadNames = new HashSet<>();
    }

    @Override
    public void callVariation(String normalPath, String tumorPath, String refGenome) throws IOException {
        try(SamplePairReadAlignmentReader pairPileupReader = new SamplePairReadAlignmentReader(normalPath, tumorPath, refGenome);
            IndexedFastaSequenceFile referenceWalker = new IndexedFastaSequenceFile(new File(refGenome));){
            Map<Integer, List<CalledVariation>> variationPerPos = new HashMap<>();
            Set<String> seenReads = new HashSet<>();
            VariationClassifier classifier = new VariationClassifier();
            int p = 1;
            for (PairedPileup pileup : pairPileupReader){
                int pos = pileup.getReferencePos();
                classifier.classifyVariationInPairedPileup(variationPerPos, pileup, seenReads, referenceWalker);
                processPotentialGermlines(variationPerPos.get(pos));
                // DEBUG
                // testingClassifier(variationPerPos, pos);
                // DEBUG
                variationPerPos.remove(pos-SLIDING_WINDOW_LIMIT);
                if (p==SLIDING_WINDOW_LIMIT) p = 0;
                p++;
            }
        }
    }

    private void processPotentialGermlines(List<CalledVariation> variationInPos) {
        for (CalledVariation var : variationInPos){
            if (!SomaticVariationType.TUMORAL_NORMAL_VARIANT.equals(var.getSomaticVariationType())) continue;
            // DEBUG
            System.out.print(var.toString());
            System.out.print("\tpresent in reads:\t");
            // DEBUG
            Map<String, Integer> supportingReads = var.getSupportingReads();
            for (Map.Entry<String, Integer> entry : supportingReads.entrySet()){
                String[] keyElems = entry.getKey().split(READ_PAIR_NAME_SEPARATOR);
                String readName = keyElems[0];
                int pairIdx = Integer.parseInt(keyElems[0]);
                Map<Integer, List<CalledVariation>> potentialGermlinesInReadPair = potentialGermlinesPerAnonRead.computeIfAbsent(readName, v -> new HashMap<>());
                List<CalledVariation> potentialGermlinesInRead = potentialGermlinesInReadPair.computeIfAbsent(pairIdx, v -> new ArrayList<>());
                potentialGermlinesInRead.add(var);
                // DEBUG
                System.out.println(entry.getKey() + " in_read_pos=" + entry.getValue());
                // DEBUG
            }
        }
    }

    @Override
    public void anonymizeReads(String normalPath, String tumorPath, String refGenome, String outputPrefix, boolean compressed) throws IOException {
        FastqWriter[][] writers = null;
        IndexedFastaSequenceFile referenceWalker = null;
        SamReader normalSamReader = null;
        SamReader tumoralSamReader = null;
        try {
            //Open streams and anonymize
            referenceWalker = new IndexedFastaSequenceFile(new File(refGenome));
            SamReaderFactory factory = SamReaderFactory.makeDefault();
            normalSamReader = factory.open(new File(normalPath));
            tumoralSamReader = factory.open(new File(tumorPath));
            FastqWriterFactory[][] factories = new FastqWriterFactory[OUTPUT_FILE_NUMBER / 2][OUTPUT_FILE_NUMBER / 2];
            writers = new FastqWriter[OUTPUT_FILE_NUMBER / 2][OUTPUT_FILE_NUMBER / 2];
            openOutputStreams(outputPrefix, compressed, factories, writers);
            anonymizeReadsInFile(normalSamReader, writers,true);
            anonymizeReadsInFile(tumoralSamReader, writers,false);
        } finally {
            // Close every stream
            assert referenceWalker != null;
            assert writers != null;
            referenceWalker.close();
            normalSamReader.close();
            tumoralSamReader.close();
            closeOutputStreams(writers);
        }
    }

    /**
     *
     * @param samReader
     * @param writers
     * @param isNormalDataset
     * @param removeUnmapped
     */
    private void anonymizeReadsInFile(SamReader samReader, FastqWriter[][] writers, boolean isNormalDataset, boolean removeUnmapped) {
        for (SAMRecord samRecord : samReader){
            String readName = samRecord.getReadName();
            if ((removeUnmapped && samRecord.getReadUnmappedFlag()) || !potentialGermlinesPerAnonRead.containsKey(readName)) continue;
            ShortAnonymizedRead anonRead = new ShortAnonymizedRead(samRecord);
            ShortAnonymizedReadPair pair;
            if (!anonReadContainer.containsKey(readName)){
                pair = new ShortAnonymizedReadPair(anonRead);
                anonReadContainer.put(readName, pair);
            }
            else{
                pair = anonReadContainer.get(readName);
                anonRead = pair.addOrUpdatePair(anonRead);
            }
            int pairIdx = anonRead.getPairIdx();
            List<CalledVariation> variantsToAnonymizeInRead = potentialGermlinesPerAnonRead.get(readName).getOrDefault(pairIdx, new ArrayList<>());
            if (anonRead.variantsToAnonymizeIsEmpty() && !variantsToAnonymizeInRead.isEmpty()) anonRead.addAllVariantsToAnonymize(variantsToAnonymizeInRead);
            anonRead.updateIfPossible(samRecord);
            if(!anonRead.isSupplementaryOrSecondary()) {
                anonRead.anonymizeVariantsInRead();
            }
            if (pair.isWriteable()){
                writeFastqRecord(pair, writers, isNormalDataset);
                potentialGermlinesPerAnonRead.remove(readName);
                anonReadContainer.remove(readName);
                anonReadNames.add(readName);
            }
        }
    }

    private void anonymizeReadsInFile(SamReader samReader, FastqWriter[][] writers, boolean isNormalDataset) {
        anonymizeReadsInFile(samReader, writers, isNormalDataset, true);
    }

    @Override
    public Map<String, ShortAnonymizedReadPair> getAnonymizedReadContainer() {
        return anonReadContainer;
    }

//    public void addAnonReadContainerIfAbsent(AnonymizedReadContainer newAnonRead){
//        String readName = newAnonRead.getReadName();
//        anonReadContainer.putIfAbsent(readName, newAnonRead);
//    }

    private void writeFastqRecord(ShortAnonymizedReadPair pairToWrite, FastqWriter[][] writers, boolean isNormalDataset){
        int datasetIdx = isNormalDataset ? NORMAL_DATASET_IDX : TUMORAL_DATASET_IDX;
        FastqRecord[] fqRecords = pairToWrite.getFastqRecords();
        writers[datasetIdx][PAIR_1_IDX].write(fqRecords[PAIR_1_IDX]);
        writers[datasetIdx][PAIR_2_IDX].write(fqRecords[PAIR_2_IDX]);
    }

    private void openOutputStreams(String prefix, boolean compressed, FastqWriterFactory[][] factories, FastqWriter[][] writers) throws IOException{
        factories[NORMAL_DATASET_IDX][PAIR_1_IDX] = new FastqWriterFactory();
        factories[NORMAL_DATASET_IDX][PAIR_2_IDX] = new FastqWriterFactory();
        factories[TUMORAL_DATASET_IDX][PAIR_1_IDX] = new FastqWriterFactory();
        factories[TUMORAL_DATASET_IDX][PAIR_2_IDX] = new FastqWriterFactory();
        File normalPair1File = new File(getFastqOutputName(prefix, NORMAL_DATASET_IDX, PAIR_1_IDX, compressed));
        File normalPair2File = new File(getFastqOutputName(prefix, NORMAL_DATASET_IDX, PAIR_2_IDX, compressed));
        File tumoralPair1File = new File(getFastqOutputName(prefix, TUMORAL_DATASET_IDX, PAIR_1_IDX, compressed));
        File tumoralPair2File = new File(getFastqOutputName(prefix, TUMORAL_DATASET_IDX, PAIR_2_IDX, compressed));
        writers[NORMAL_DATASET_IDX][PAIR_1_IDX] = factories[NORMAL_DATASET_IDX][PAIR_1_IDX].newWriter(normalPair1File);
        writers[NORMAL_DATASET_IDX][PAIR_2_IDX] = factories[NORMAL_DATASET_IDX][PAIR_2_IDX].newWriter(normalPair2File);
        writers[TUMORAL_DATASET_IDX][PAIR_1_IDX] = factories[TUMORAL_DATASET_IDX][PAIR_1_IDX].newWriter(tumoralPair1File);
        writers[TUMORAL_DATASET_IDX][PAIR_2_IDX] = factories[TUMORAL_DATASET_IDX][PAIR_2_IDX].newWriter(tumoralPair2File);
    }

    private void closeOutputStreams(FastqWriter[][] writers){
        writers[NORMAL_DATASET_IDX][PAIR_1_IDX].close();
        writers[NORMAL_DATASET_IDX][PAIR_2_IDX].close();
        writers[TUMORAL_DATASET_IDX][PAIR_1_IDX].close();
        writers[TUMORAL_DATASET_IDX][PAIR_2_IDX].close();
    }

    public String getFastqOutputName(String outputPrefix, int datasetIdx, int pairIdx, boolean compressed){
        String datasetIdStr = datasetIdx == NORMAL_DATASET_IDX ? ".N" : ".T";
        String pairIdxStr = pairIdx == PAIR_1_IDX ? ".1" : ".2";
        String extension = compressed ? ".fastq.gz" : ".fastq";
        return outputPrefix + datasetIdStr + pairIdxStr + extension;
    }

    // DEBUG
    public void testingClassifier(Map<Integer, List<CalledVariation>> variationPerPos, int pos){
        List<CalledVariation> variationInPos = variationPerPos.get(pos);
        // System.out.println("pos " + pos);
        // if (variationInPos == null) {
        //     System.out.println("#POS_BUG: " + pos);
        // }
        for (CalledVariation var : variationInPos){
            System.out.print(var.toString());
            System.out.print("\tpresent in reads:\t");
            for (Map.Entry<String, Integer> entry : var.getSupportingReads().entrySet()){
                System.out.println(entry.getKey() + " in_read_pos=" + entry.getValue());
            }
        }
    }
    // DEBUG
}
