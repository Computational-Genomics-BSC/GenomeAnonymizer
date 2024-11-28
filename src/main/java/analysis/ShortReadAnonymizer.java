package analysis;

import genomicelements.*;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.logging.Logger;
import java.util.stream.Collectors;

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

    List<GenomicRegion> partitions;
    //Set that contaains all the reads that will be excluded from the result (e.g. Unmapped and MAPQ=0)
    Set<String> readsToExclude;
    // Written read alignments to avoid writing again
    Set<String> writtenReadAlignments;
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
        writtenReadAlignments = new HashSet<>();
        readGermlinesToAnonymize = new HashMap<>();
        factory = SamReaderFactory.makeDefault();
        removeUnmapped = true;
        writersAreOpen = false;
    }

    public void setPartitions(List<GenomicRegion> partitions) {
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
        //TODO: Check if it is possible to change partitions based on actual content
        // (Implement CoveredGenomicRegion), to improve runtime using parallelization
        for(GenomicRegion partition : partitions){
            for(String path : paths){
                CompletableFuture<Set<String>> future = CompletableFuture.supplyAsync (() -> {
                    Set<String> answer;
                    try {
                        answer = queryReadsToExcludeInPartition(path, partition);
                    }
                    catch (IOException e) {
                        LOGGER.severe("Exception in thread querying reads to exclude in region: "
                                + partition.getSequenceName()
                                + " " + partition.getStart() + " " + partition.getEnd());
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

    private Set<String> queryReadsToExcludeInPartition(String filePath, GenomicRegion partition) throws IOException{
        Set<String> readsToExcludeInPartition = new HashSet<>();
        try(SamReader reader = factory.open(new File(filePath))){
            SAMRecordIterator it = reader.query(partition.getSequenceName(), partition.getStart(), partition.getEnd(), false);
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
        try(IndexedFastaSequenceFile referenceWalker = new IndexedFastaSequenceFile(new File(refGenome))) {
            openOutputStreams(outputPrefix);
            anonymizeReadsInFile(new File(normalPath), true, referenceWalker);
            anonymizeReadsInFile(new File(tumorPath), false, referenceWalker);
        } finally {
            closeOutputStreams();
        }
    }

    private void openOutputStreams(String prefix){
        normalOutputFile = new File(getBAMOutputName(prefix, NORMAL_DATASET_IDX));
        tumoralOutputFile = new File(getBAMOutputName(prefix, TUMORAL_DATASET_IDX));
        SAMFileWriterFactory factory = new SAMFileWriterFactory();
        factory.setCreateIndex(true);
        normalWriter = factory.makeBAMWriter(normalFileHeader, true, normalOutputFile);
        tumoralWriter = factory.makeBAMWriter(tumoralFileHeader, true, tumoralOutputFile);
        writersAreOpen = true;
    }

    public String getBAMOutputName(String outputPrefix, int datasetIdx){
        String anonTag = ".anonymized";
        String datasetIdStr = datasetIdx == NORMAL_DATASET_IDX ? ".N" : ".T";
        //String extension = compressed ? ".BAM" : ".SAM";
        String extension = BAM_FILE;
        return outputPrefix + anonTag + datasetIdStr + extension;
    }

    /**
     * @param samFile
     * @param isNormalDataset
     * @param referenceGenome
     */
    private void anonymizeReadsInFile(File samFile, boolean isNormalDataset,
                                      IndexedFastaSequenceFile referenceGenome) throws IOException{
        String currentContig = "";
        byte[] referenceContigSequence = new byte[0];
        partitions = partitions.stream()
                .sorted(Comparator.comparing(GenomicRegion::getSequenceIdx)
                        .thenComparing(GenomicRegion::getStart))
                .collect(Collectors.toList());
        for(int i = 0; i < partitions.size(); i++){
            GenomicRegion partition = partitions.get(i);
            if(!currentContig.equals(partition.getSequenceName())){
                currentContig = partition.getSequenceName();
                referenceContigSequence = referenceGenome.getSequence(currentContig).getBases();
            }
            try(SamReader samReader = factory.open(samFile)){
                SAMRecordIterator it = samReader.query(partition.getSequenceName(), partition.getStart(),
                        partition.getEnd(), false);
                while (it.hasNext()) {
                    SAMRecord samRecord = it.next();
                    String readName = samRecord.getReadName();
                    String readAlnId = generateReadId(samRecord);
                    if(readsToExclude.contains(readName) ||
                            writtenReadAlignments.contains(readAlnId)) continue;
                    SAMRecord newSamRecord = samRecord;
                    if(readGermlinesToAnonymize.containsKey(readAlnId)){
                        ShortReadAlignment readAlignment = new ShortReadAlignment(samRecord);
                        ShortAnonymizedReadAlignment anonymizedReadAlignment = new ShortAnonymizedReadAlignment(readAlignment);
                        anonymizedReadAlignment.setReferenceContigSequence(referenceContigSequence);
                        anonymizedReadAlignment.setVariantsToAnonymize(readGermlinesToAnonymize.get(readAlnId));
                        anonymizedReadAlignment.anonymizeVariants();
                        newSamRecord = anonymizedReadAlignment.getAnonymizedSamRecord();
                    }
                    writeReadAlignment(newSamRecord, isNormalDataset, readAlnId);
                }
            }
        }
    }

    private void writeReadAlignment(SAMRecord readAln, boolean isNormalDataset, String readAlnId){
        if(isNormalDataset) normalWriter.addAlignment(readAln);
        else tumoralWriter.addAlignment(readAln);
        writtenReadAlignments.add(readAlnId);
    }

    private void closeOutputStreams(){
        normalWriter.close();
        tumoralWriter.close();
        writersAreOpen = false;
    }
}
