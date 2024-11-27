package analysis;

import genomicelements.*;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.SimpleFeature;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.logging.Logger;

import static analysis.GenomeAnonymizer.BAM_FILE;
import static genomicelements.ShortReadAlignment.generateReadId;

/**
 * AnonymizerAlgorithm implementation for short read data
 * @author Nicolas Gaitan
 */
public class ShortReadAnonymizer implements AnonymizerAlgorithm{

    private static final Logger LOGGER = Logger.getLogger(ShortReadAnonymizer.class.getName());

    public static final int NORMAL_DATASET_IDX = 0;
    public static final int TUMORAL_DATASET_IDX = 1;
    public static final int OUTPUT_FILE_NUMBER = 4;

    List<SimpleFeature> partitions;
    //Set that contaains all the reads that will be excluded from the result (e.g. Unmapped and MAPQ=0)
    Set<String> readsToExclude;
    // Map containing all potential germlines (value: List), per pair (nested key, 0 or 1), per read (key)
    Map<String, List<CalledVariation>> readGermlinesToAnonymize;
    // Map collecting to-be-anonymized reads while they can be masked and written
    Map<String, ShortAnonymizedReadPair> anonReadContainer;
    Set<String> anonReadNames;
    SAMFileHeader normalFileHeader;
    SAMFileHeader tumoralFileHeader;
    File normalOutputFile;
    File tumoralOutputFile;
    SAMFileWriter normalWriter;
    SAMFileWriter tumoralWriter;
    boolean removeUnmapped;
    boolean writersAreOpen;

    public ShortReadAnonymizer() {
        readsToExclude = new HashSet<>();
        anonReadContainer = new HashMap<>();
        readGermlinesToAnonymize = new HashMap<>();
        anonReadNames = new HashSet<>();
        //outputFiles = new File[OUTPUT_FILE_NUMBER/2][OUTPUT_FILE_NUMBER/2];
        //writers = new FastqWriter[OUTPUT_FILE_NUMBER/2][OUTPUT_FILE_NUMBER/2];
        removeUnmapped = true;
        writersAreOpen = false;
    }

    public void setPartitions(List<SimpleFeature> partitions) {
        this.partitions = partitions;
    }

    public void setReadGermlinesToAnonymize(Map<String, List<CalledVariation>> readGermlinesToAnonymize) {
        this.readGermlinesToAnonymize = readGermlinesToAnonymize;
    }

    public void queryReadsToExclude(String normalPath, String tumorPath, int threads) throws Exception{
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

    private Set<String> queryReadsToExcludeInPartition(String filePath, SamReaderFactory factory,
                                                       SimpleFeature partition) throws IOException{
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
    public void anonymizeReads(String normalPath, String tumorPath, String refGenome, String outputPrefix,
                               boolean compressed) throws IOException {
        IndexedFastaSequenceFile referenceWalker = null;
        SamReader normalSamReader = null;
        SamReader tumoralSamReader = null;
        try {
            //Open streams and anonymize
            referenceWalker = new IndexedFastaSequenceFile(new File(refGenome));
            SamReaderFactory factory = SamReaderFactory.makeDefault();
            normalSamReader = factory.open(new File(normalPath));
            tumoralSamReader = factory.open(new File(tumorPath));
            normalFileHeader = normalSamReader.getFileHeader();
            tumoralFileHeader = tumoralSamReader.getFileHeader();
            openOutputStreams(outputPrefix);
            anonymizeReadsInFile(normalSamReader, true, referenceWalker);
            anonymizeReadsInFile(tumoralSamReader, false, referenceWalker);
        }
        finally {
            // Close every stream
            referenceWalker.close();
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
    private void anonymizeReadsInFile(SamReader samReader, boolean isNormalDataset,
                                      IndexedFastaSequenceFile referenceGenome) {
        String currentContig = "";
        byte[] referenceContigSequence = new byte[0];
        for (SimpleFeature partition : partitions){
            if(!currentContig.equals(partition.getContig())){
                currentContig = partition.getContig();
                referenceContigSequence = referenceGenome.getSubsequenceAt(partition.getContig(),
                        partition.getStart(), partition.getEnd()).getBases();
            }
            SAMRecordIterator it = samReader.query(partition.getContig(), partition.getStart(), partition.getEnd(), false);
            while (it.hasNext()) {
                SAMRecord samRecord = it.next();
                String readName = samRecord.getReadName();
                if(readsToExclude.contains(readName)) continue;
                String readId = generateReadId(samRecord);
                SAMRecord newSamRecord = samRecord;
                if(readGermlinesToAnonymize.containsKey(readId)){
                    ShortReadAlignment readAlignment = new ShortReadAlignment(samRecord);
                    ShortAnonymizedReadAlignment anonymizedReadAlignment = new ShortAnonymizedReadAlignment(readAlignment);
                    anonymizedReadAlignment.setReferenceContigSequence(referenceContigSequence);
                    anonymizedReadAlignment.setVariantsToAnonymize(readGermlinesToAnonymize.get(readId));
                    anonymizedReadAlignment.anonymizeVariants();
                    newSamRecord = anonymizedReadAlignment.getAnonymizedSamRecord();
                }
                if(isNormalDataset) normalWriter.addAlignment(newSamRecord);
                else tumoralWriter.addAlignment(newSamRecord);
            }
        }
    }

    public Map<String, ShortAnonymizedReadPair> getAnonymizedReadContainer() {
        return anonReadContainer;
    }

    private void openOutputStreams(String prefix){
        normalOutputFile = new File(getBAMOutputName(prefix, NORMAL_DATASET_IDX));
        tumoralOutputFile = new File(getBAMOutputName(prefix, TUMORAL_DATASET_IDX));
        SAMFileWriterFactory factory = new SAMFileWriterFactory();
        normalWriter = factory.makeBAMWriter(normalFileHeader, true, normalOutputFile);
        tumoralWriter = factory.makeBAMWriter(tumoralFileHeader, true, tumoralOutputFile);
        writersAreOpen = true;
    }

    private void closeOutputStreams(){
        normalWriter.close();
        tumoralWriter.close();
        writersAreOpen = false;
    }

    public String getBAMOutputName(String outputPrefix, int datasetIdx){
        String anonTag = ".anonymized";
        String datasetIdStr = datasetIdx == NORMAL_DATASET_IDX ? ".N" : ".T";
        //String extension = compressed ? ".BAM" : ".SAM";
        String extension = BAM_FILE;
        return outputPrefix + anonTag + datasetIdStr + extension;
    }
}
