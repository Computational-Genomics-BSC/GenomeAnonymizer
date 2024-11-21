package analysis;

import genomicelements.*;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.SimpleFeature;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.logging.Logger;
import static genomicelements.ShortAnonymizedReadPair.*;

/**
 * AnonymizerAlgorithm implementation for short read data
 * @author Nicolas Gaitan
 */
public class ShortReadAnonymizer implements AnonymizerAlgorithm{

    private static final Logger LOGGER = Logger.getLogger(ShortReadAnonymizer.class.getName());

    public static final int NORMAL_DATASET_IDX = 0;
    public static final int TUMORAL_DATASET_IDX = 1;
    public static final int OUTPUT_FILE_NUMBER = 4;

    //Set that contaains all the reads that will be excluded from the result (e.g. Unmapped and MAPQ=0)
    Set<String> readsToExclude;
    // Map containing all potential germlines (value: List), per pair (nested key, 0 or 1), per read (key)
    Map<String, List<CalledVariation>> readGermlinesToAnonymize;
    // Map collecting to-be-anonymized reads while they can be masked and written
    Map<String, ShortAnonymizedReadPair> anonReadContainer;
    Set<String> anonReadNames;
    File[][] outputFiles;
    FastqWriter[][] writers;
    boolean removeUnmapped;
    boolean writersAreOpen;

    public ShortReadAnonymizer() {
        readsToExclude = new HashSet<>();
        anonReadContainer = new HashMap<>();
        readGermlinesToAnonymize = new HashMap<>();
        anonReadNames = new HashSet<>();
        outputFiles = new File[OUTPUT_FILE_NUMBER/2][OUTPUT_FILE_NUMBER/2];
        writers = new FastqWriter[OUTPUT_FILE_NUMBER/2][OUTPUT_FILE_NUMBER/2];
        removeUnmapped = true;
        writersAreOpen = false;
    }

    public void setReadGermlinesToAnonymize(Map<String, List<CalledVariation>> readGermlinesToAnonymize) {
        this.readGermlinesToAnonymize = readGermlinesToAnonymize;
    }

    public void queryReadsToExclude(String normalPath, String tumorPath, List<SimpleFeature> partitions, int threads) throws Exception{
        String[] paths = new String[2];
        paths[0] = normalPath;
        paths[1] = tumorPath;
        SamReaderFactory factory = SamReaderFactory.makeDefault();
        ExecutorService exec = Executors.newFixedThreadPool(threads);
        List<CompletableFuture<Set<String>>> futures = new ArrayList<>();
        for(SimpleFeature partition : partitions){
            for(String path : paths){
                CompletableFuture<Set<String>> future = CompletableFuture.supplyAsync (() -> {
                    Set<String> answer;
                    try {
                        answer = queryReadsToExcludeInPartition(path, factory, partition);
                    }
                    catch (IOException e) {
                        LOGGER.severe("Exception in thread querying reads to exclude in region: "
                                + partition.getContig() + " " + partition.getStart() + " " + partition.getEnd());
                        throw new RuntimeException(e);
                    }
                    return answer;
                }, exec);
                futures.add(future);
            }
        }
        CompletableFuture<Void> allFutures = CompletableFuture.allOf(futures.toArray(new CompletableFuture[0]));
        try{
            for (CompletableFuture<Set<String>> future : futures) {
                Set<String> answer = future.get();
                readsToExclude.addAll(answer);
            }
        }
        catch (Exception e){
            LOGGER.severe("Exception when retrieving excluded reads from thread halting execution prematurely");
            LOGGER.severe(e.getMessage());
            throw new RuntimeException(e);
        }
        allFutures.join();
        exec.shutdown();
    }

    private Set<String> queryReadsToExcludeInPartition(String filePath, SamReaderFactory factory,  SimpleFeature partition) throws IOException{
        Set<String> readsToExcludeInPartition = new HashSet<>();
        try(SamReader reader = factory.open(new File(filePath))){
            SAMRecordIterator it = reader.query(partition.getContig(), partition.getStart(), partition.getEnd(), false);
            while (it.hasNext()) {
                SAMRecord samRecord = it.next();
                if(samRecord.getReadUnmappedFlag() ||
                        (samRecord.getMappingQuality()==0 && !samRecord.isSecondaryOrSupplementary())){
                    readsToExcludeInPartition.add(samRecord.getReadName());
                }
            }
        }
        return readsToExcludeInPartition;
    }

    @Override
    public Set<String> getReadsToExclude() {
        return readsToExclude;
    }

    @Override
    public void anonymizeReads(String normalPath, String tumorPath, String refGenome, String outputPrefix, boolean compressed) throws IOException {
        //IndexedFastaSequenceFile referenceWalker = null;
        SamReader normalSamReader = null;
        SamReader tumoralSamReader = null;
        try {
            //Open streams and anonymize
            //referenceWalker = new IndexedFastaSequenceFile(new File(refGenome));
            SamReaderFactory factory = SamReaderFactory.makeDefault();
            normalSamReader = factory.open(new File(normalPath));
            tumoralSamReader = factory.open(new File(tumorPath));
            writers = new FastqWriter[OUTPUT_FILE_NUMBER / 2][OUTPUT_FILE_NUMBER / 2];
            createOutputStreams(outputPrefix, compressed);
            anonymizeReadsInFile(normalSamReader, true);
            anonymizeReadsInFile(tumoralSamReader, false);
            writeUnmodifiedReads(normalPath, tumorPath, outputPrefix, compressed);
        } finally {
            // Close every stream
            //assert referenceWalker != null;
            assert writers != null;
            //referenceWalker.close();
            normalSamReader.close();
            tumoralSamReader.close();
            closeOutputStreams();
        }
    }

    /**
     *
     * @param samReader
     * @param isNormalDataset
     */
    private void anonymizeReadsInFile(SamReader samReader, boolean isNormalDataset) {
        for (SAMRecord samRecord : samReader){
            String readName = samRecord.getReadName();
            if ((removeUnmapped && samRecord.getReadUnmappedFlag()) || !readGermlinesToAnonymize.containsKey(readName)) continue;
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
            List<CalledVariation> variantsToAnonymizeInRead = readGermlinesToAnonymize.get(readName).getOrDefault(pairIdx, new ArrayList<>());
            if (anonRead.variantsToAnonymizeIsEmpty() && !variantsToAnonymizeInRead.isEmpty()) anonRead.addAllVariantsToAnonymize(variantsToAnonymizeInRead);
            anonRead.updateIfPossible(samRecord);
            if(!anonRead.isSupplementaryOrSecondary() && !anonRead.isAnonymized()) {
                anonRead.anonymizeVariantsInRead();
            }
            if (pair.isWriteable()){
                writeFastqRecord(pair, isNormalDataset);
                readGermlinesToAnonymize.remove(readName);
                anonReadContainer.remove(readName);
                anonReadNames.add(readName);
            }
        }
    }


//    @Override
//    public void anonymizeReads(String normalPath, String tumorPath, String refGenome, String outputPrefix, boolean compressed) throws IOException {
//        //IndexedFastaSequenceFile referenceWalker = null;
//        SamReader normalSamReader = null;
//        SamReader tumoralSamReader = null;
//        try {
//            //Open streams and anonymize
//            //referenceWalker = new IndexedFastaSequenceFile(new File(refGenome));
//            SamReaderFactory factory = SamReaderFactory.makeDefault();
//            normalSamReader = factory.open(new File(normalPath));
//            tumoralSamReader = factory.open(new File(tumorPath));
//            writers = new FastqWriter[OUTPUT_FILE_NUMBER / 2][OUTPUT_FILE_NUMBER / 2];
//            createOutputStreams(outputPrefix, compressed);
//            anonymizeReadsInFile(normalSamReader, true);
//            anonymizeReadsInFile(tumoralSamReader, false);
//            writeUnmodifiedReads(normalPath, tumorPath, outputPrefix, compressed);
//        } finally {
//            // Close every stream
//            //assert referenceWalker != null;
//            assert writers != null;
//            //referenceWalker.close();
//            normalSamReader.close();
//            tumoralSamReader.close();
//            closeOutputStreams();
//        }
//    }


//
//    /**
//     *
//     * @param samReader
//     * @param isNormalDataset
//     */
//    private void anonymizeReadsInFile(SamReader samReader, boolean isNormalDataset) {
//        for (SAMRecord samRecord : samReader){
//            String readName = samRecord.getReadName();
//            if ((removeUnmapped && samRecord.getReadUnmappedFlag()) || !readGermlinesToAnonymize.containsKey(readName)) continue;
//            ShortAnonymizedRead anonRead = new ShortAnonymizedRead(samRecord);
//            ShortAnonymizedReadPair pair;
//            if (!anonReadContainer.containsKey(readName)){
//                pair = new ShortAnonymizedReadPair(anonRead);
//                anonReadContainer.put(readName, pair);
//            }
//            else{
//                pair = anonReadContainer.get(readName);
//                anonRead = pair.addOrUpdatePair(anonRead);
//            }
//            int pairIdx = anonRead.getPairIdx();
//            List<CalledVariation> variantsToAnonymizeInRead = readGermlinesToAnonymize.get(readName).getOrDefault(pairIdx, new ArrayList<>());
//            if (anonRead.variantsToAnonymizeIsEmpty() && !variantsToAnonymizeInRead.isEmpty()) anonRead.addAllVariantsToAnonymize(variantsToAnonymizeInRead);
//            anonRead.updateIfPossible(samRecord);
//            if(!anonRead.isSupplementaryOrSecondary() && !anonRead.isAnonymized()) {
//                anonRead.anonymizeVariantsInRead();
//            }
//            if (pair.isWriteable()){
//                writeFastqRecord(pair, isNormalDataset);
//                readGermlinesToAnonymize.remove(readName);
//                anonReadContainer.remove(readName);
//                anonReadNames.add(readName);
//            }
//        }
//    }

    public void writeUnmodifiedReads(String normalPath, String tumorPath, String outputPrefix, boolean compressed) throws IOException{
        if(!writersAreOpen) openOutputStreams();
        SamReaderFactory factory = SamReaderFactory.makeDefault();
        try(SamReader normalSamReader = factory.open(new File(normalPath));
            SamReader tumoralSamReader = factory.open(new File(tumorPath));){
            handleUnmodifiedReadsInFile(normalSamReader, true);
            handleUnmodifiedReadsInFile(tumoralSamReader, false);
        }
        finally {
            closeOutputStreams();
        }
    }

    private void handleUnmodifiedReadsInFile(SamReader samReader, boolean isNormalDataset){
        Map<String, SAMRecord[]> unmodifiedReadsContainer = new HashMap<>();
        int datasetIdx = isNormalDataset ? NORMAL_DATASET_IDX : TUMORAL_DATASET_IDX;
        for (SAMRecord samRecord : samReader){
            String readName = samRecord.getReadName();
            if ((removeUnmapped && samRecord.getReadUnmappedFlag()) || anonReadNames.contains(readName) || samRecord.isSecondaryOrSupplementary()) continue;
            SAMRecord[] arrayPair = unmodifiedReadsContainer.computeIfAbsent(readName, v -> new SAMRecord[2]);
            int pairIdx = samRecord.getFirstOfPairFlag() ? PAIR_1_IDX : PAIR_2_IDX;
            arrayPair[pairIdx] = samRecord;
            if(arrayPair[PAIR_1_IDX] != null && arrayPair[PAIR_2_IDX] != null){
                writeFastqRecord(arrayPair, isNormalDataset);
                unmodifiedReadsContainer.remove(readName);
            }
        }
    }

    public Map<String, ShortAnonymizedReadPair> getAnonymizedReadContainer() {
        return anonReadContainer;
    }

    private void writeFastqRecord(ShortAnonymizedReadPair pairToWrite, boolean isNormalDataset){
        int datasetIdx = isNormalDataset ? NORMAL_DATASET_IDX : TUMORAL_DATASET_IDX;
        FastqRecord[] fqRecords = pairToWrite.getFastqRecords();
        writers[datasetIdx][PAIR_1_IDX].write(fqRecords[PAIR_1_IDX]);
        writers[datasetIdx][PAIR_2_IDX].write(fqRecords[PAIR_2_IDX]);
    }

    private void writeFastqRecord(SAMRecord[] pairToWrite, boolean isNormalDataset){
        int datasetIdx = isNormalDataset ? NORMAL_DATASET_IDX : TUMORAL_DATASET_IDX;
        String comment = "";
        FastqRecord[] fqRecords = new FastqRecord[2];
        fqRecords[PAIR_1_IDX] = getFastqRecordFromSamRecord(pairToWrite[PAIR_1_IDX]);
        fqRecords[PAIR_2_IDX] = getFastqRecordFromSamRecord(pairToWrite[PAIR_2_IDX]);
        writers[datasetIdx][PAIR_1_IDX].write(fqRecords[PAIR_1_IDX]);
        writers[datasetIdx][PAIR_2_IDX].write(fqRecords[PAIR_2_IDX]);
    }

    private FastqRecord getFastqRecordFromSamRecord(SAMRecord samRec){
        byte[] sequenceArray = samRec.getReadBases().clone();
        byte[] qualitiesArray = samRec.getBaseQualities().clone();
        if (samRec.getReadNegativeStrandFlag()){
            SequenceUtil.reverseComplement(sequenceArray);
            SequenceUtil.reverseQualities(qualitiesArray);
        }
        String name = samRec.getReadName();
        //String readSequence = new String(sequenceArray, StandardCharsets.UTF_8);
        //String qualitySequence = new String(qualitiesArray, StandardCharsets.UTF_8);
        String comment = "";
        return new FastqRecord(name, sequenceArray, comment, qualitiesArray);
    }

    private void createOutputStreams(String prefix, boolean compressed) throws IOException{
        File normalPair1File = new File(getFastqOutputName(prefix, NORMAL_DATASET_IDX, PAIR_1_IDX, compressed));
        File normalPair2File = new File(getFastqOutputName(prefix, NORMAL_DATASET_IDX, PAIR_2_IDX, compressed));
        File tumoralPair1File = new File(getFastqOutputName(prefix, TUMORAL_DATASET_IDX, PAIR_1_IDX, compressed));
        File tumoralPair2File = new File(getFastqOutputName(prefix, TUMORAL_DATASET_IDX, PAIR_2_IDX, compressed));
        outputFiles[NORMAL_DATASET_IDX][PAIR_1_IDX] = normalPair1File;
        outputFiles[NORMAL_DATASET_IDX][PAIR_2_IDX] = normalPair2File;
        outputFiles[TUMORAL_DATASET_IDX][PAIR_1_IDX] = tumoralPair1File;
        outputFiles[TUMORAL_DATASET_IDX][PAIR_2_IDX] = tumoralPair2File;
        openOutputStreams();
    }

    private void openOutputStreams() throws IOException {
        if(writersAreOpen) return;
        FastqWriterFactory factory = new FastqWriterFactory();
        writers[NORMAL_DATASET_IDX][PAIR_1_IDX] = factory.newWriter(outputFiles[NORMAL_DATASET_IDX][PAIR_1_IDX]);
        writers[NORMAL_DATASET_IDX][PAIR_2_IDX] = factory.newWriter(outputFiles[NORMAL_DATASET_IDX][PAIR_2_IDX]);
        writers[TUMORAL_DATASET_IDX][PAIR_1_IDX] = factory.newWriter(outputFiles[TUMORAL_DATASET_IDX][PAIR_1_IDX]);
        writers[TUMORAL_DATASET_IDX][PAIR_2_IDX] = factory.newWriter(outputFiles[TUMORAL_DATASET_IDX][PAIR_2_IDX]);
        writersAreOpen = true;
    }

//    private void createOutputStreams() throws IOException{
//        writers[NORMAL_DATASET_IDX][PAIR_1_IDX].;
//        writers[NORMAL_DATASET_IDX][PAIR_2_IDX] = factory.newWriter(normalPair2File);
//        writers[TUMORAL_DATASET_IDX][PAIR_1_IDX] = factory.newWriter(tumoralPair1File);
//        writers[TUMORAL_DATASET_IDX][PAIR_2_IDX] = factory.newWriter(tumoralPair2File);
//    }

    private void closeOutputStreams(){
        writers[NORMAL_DATASET_IDX][PAIR_1_IDX].close();
        writers[NORMAL_DATASET_IDX][PAIR_2_IDX].close();
        writers[TUMORAL_DATASET_IDX][PAIR_1_IDX].close();
        writers[TUMORAL_DATASET_IDX][PAIR_2_IDX].close();
        writersAreOpen = false;
    }

    public String getFastqOutputName(String outputPrefix, int datasetIdx, int pairIdx, boolean compressed){
        String anonTag = ".anonymized";
        String datasetIdStr = datasetIdx == NORMAL_DATASET_IDX ? ".N" : ".T";
        String pairIdxStr = pairIdx == PAIR_1_IDX ? ".1" : ".2";
        String extension = compressed ? ".fastq.gz" : ".fastq";
        return outputPrefix + anonTag + datasetIdStr + pairIdxStr + extension;
    }

    /**
     * Must be set before calling any public method
     * @param removeUnmapped
     */
    public void setRemoveUnmapped(boolean removeUnmapped){
        this.removeUnmapped = removeUnmapped;
    }

    // DEBUG
//    public void testingClassifier(Map<Integer, List<CalledVariation>> variationPerPos, int pos){
//        List<CalledVariation> variationInPos = variationPerPos.get(pos);
//        // System.out.println("pos " + pos);
//        // if (variationInPos == null) {
//        //     System.out.println("#POS_BUG: " + pos);
//        // }
//        for (CalledVariation var : variationInPos){
//            System.out.print(var.toString());
//            System.out.print("\tpresent in reads:\t");
//            for (Map.Entry<String, Integer> entry : var.getSupportingReads().entrySet()){
//                System.out.println(entry.getKey() + " in_read_pos=" + entry.getValue());
//            }
//        }
//    }
    // DEBUG
}
