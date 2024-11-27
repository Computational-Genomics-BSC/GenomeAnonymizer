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

    List<SimpleFeature> partitions;
    //Set that contaains all the reads that will be excluded from the result (e.g. Unmapped and MAPQ=0)
    Set<String> readsToExclude;
    // Map containing all potential germlines (value: List), per pair (nested key, 0 or 1), per read (key)
    Map<String, List<CalledVariation>> readGermlinesToAnonymize;
    SAMFileHeader normalFileHeader;
    SAMFileHeader tumoralFileHeader;
    File normalOutputFile;
    File tumoralOutputFile;
    SamReaderFactory factory;
    SAMFileWriter normalWriter;
    SAMFileWriter tumoralWriter;
    boolean removeUnmapped;
    boolean writersAreOpen;

    public ShortReadAnonymizer() {
        readsToExclude = new HashSet<>();
        readGermlinesToAnonymize = new HashMap<>();
        factory = SamReaderFactory.makeDefault();
        removeUnmapped = true;
        writersAreOpen = false;
    }

    public void setPartitions(List<SimpleFeature> partitions) {
        this.partitions = partitions;
    }

    public void setReadGermlinesToAnonymize(Map<String, List<CalledVariation>> readGermlinesToAnonymize) {
        this.readGermlinesToAnonymize = readGermlinesToAnonymize;
    }

    private void retrieveFileHeaders(String bamFile, boolean isNormalDataset) throws IOException {
        try (SamReader samReader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .open(new File(bamFile))) {
            // Retrieve the SAMFileHeader
            SAMFileHeader header = samReader.getFileHeader();
            if (isNormalDataset) normalFileHeader = header;
            else tumoralFileHeader = header;
        }
    }

    public void queryReadsToExclude(String normalPath, String tumorPath, int threads) throws Exception{
        String[] paths = new String[2];
        paths[0] = normalPath;
        paths[1] = tumorPath;
        ExecutorService exec = Executors.newFixedThreadPool(threads);
        List<CompletableFuture<Set<String>>> futures = new ArrayList<>();
        retrieveFileHeaders(normalPath, true);
        retrieveFileHeaders(tumorPath, false);
        for(SimpleFeature partition : partitions){
            for(String path : paths){
                CompletableFuture<Set<String>> future = CompletableFuture.supplyAsync (() -> {
                    Set<String> answer;
                    try {
                        answer = queryReadsToExcludeInPartition(path, partition);
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

    private Set<String> queryReadsToExcludeInPartition(String filePath, SimpleFeature partition) throws IOException{
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
        IndexedFastaSequenceFile referenceWalker = new IndexedFastaSequenceFile(new File(refGenome));
        //SamReader normalSamReader = null;
        //SamReader tumoralSamReader = null;
        //try {
        //Open streams and anonymize
        //normalSamReader = factory.open(new File(normalPath));
        //tumoralSamReader = factory.open(new File(tumorPath));
        //normalFileHeader = normalSamReader.getFileHeader();
        //tumoralFileHeader = tumoralSamReader.getFileHeader();
        openOutputStreams(outputPrefix);
        anonymizeReadsInFile(new File(normalPath), true, referenceWalker);
        anonymizeReadsInFile(new File(tumorPath), false, referenceWalker);
        //}
        //finally {
        // Close every stream
        referenceWalker.close();
        //normalSamReader.close();
        //tumoralSamReader.close();
        //closeOutputStreams();
        //}
    }

    /**
     *
     * @param samFile
     * @param isNormalDataset
     * @param referenceGenome
     */
    private void anonymizeReadsInFile(File samFile, boolean isNormalDataset,
                                      IndexedFastaSequenceFile referenceGenome) throws IOException{
        String currentContig = "";
        byte[] referenceContigSequence = new byte[0];
        // Sort partitions because they are in disorder, could even be out of this class
        for (int i = 0; i < partitions.size(); i++){
            SimpleFeature partition = partitions.get(i);
            if(!currentContig.equals(partition.getContig())){
                currentContig = partition.getContig();
                referenceContigSequence = referenceGenome.getSequence(partition.getContig()).getBases();
            }
            try(SamReader samReader = factory.open(samFile)){
                SAMRecordIterator it = samReader.query(partition.getContig(), partition.getStart(),
                        partition.getEnd(), false);
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
