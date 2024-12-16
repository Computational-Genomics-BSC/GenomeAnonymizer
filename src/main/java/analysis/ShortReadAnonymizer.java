package analysis;

import genomicelements.*;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import io.GenomicRegionBedReader;

import java.io.File;
import java.io.IOException;
import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.util.*;
import java.util.concurrent.Future;
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
    // Map containing all potential germlines (value: List), per pair (nested key, 0 or 1), per read (key)
    Map<String, List<CalledVariation>> readGermlinesToAnonymize;
    SAMFileHeader normalFileHeader;
    SAMFileHeader tumoralFileHeader;
    File canvasNormal;
    File canvasTumoral;
    List<GenomicRegion> genomicRegions;
    SamReaderFactory factory;
    SAMFileWriter normalWriter;
    SAMFileWriter tumoralWriter;
    boolean removeUnmapped;

    public ShortReadAnonymizer() {
        readsToExclude = new HashSet<>();
        readGermlinesToAnonymize = new HashMap<>();
        factory = SamReaderFactory.makeDefault();
        factory.setUseAsyncIo(true);
        factory.setDefaultValidationStringency(ValidationStringency.SILENT);
        removeUnmapped = true;
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
        // TODO: Make sure the reference genome is the same for this.canvasTumoral and this.tumoralWriter and for this.canvasNormal and this.normalWriter
        try {
            openOutputStreams(outputPrefix);
            ExecutorService exec = Executors.newFixedThreadPool(2);
            Future<?> tumoralFuture = exec.submit(() -> {
                try {
                    mergeAnonymizedReads(new File(tumorPath), this.canvasTumoral, this.tumoralWriter, new File(refGenome));
                } catch (Exception e) {
                    LOGGER.severe("Exception in thread merging anonymized reads in tumoral dataset");
                    LOGGER.severe(e.getMessage());
                    throw new RuntimeException(e);
                }
            });
            Future<?> normalFuture = exec.submit(() -> {
                try {
                    mergeAnonymizedReads(new File(normalPath), this.canvasNormal, this.normalWriter, new File(refGenome));
                } catch (Exception e) {
                    LOGGER.severe("Exception in thread merging anonymized reads in normal dataset");
                    LOGGER.severe(e.getMessage());
                    throw new RuntimeException(e);
                }
            });
            exec.shutdown();
            normalFuture.get();
            tumoralFuture.get();
        } catch (Exception e) {
            LOGGER.severe("Exception when merging anonymized reads halting execution prematurely");
            LOGGER.severe(e.getMessage());
            throw new IOException(e);
        } finally {
            closeOutputStreams();
        }
    }

    @Override
    public void setRegions(String bedFilePath) throws IOException {
        this.genomicRegions = GenomicRegionBedReader.readGenomicRegionBED(bedFilePath);
    }

    @Override
    public void setCanvasFiles(String normalCanvasFileName, String tumoralCanvasFileName) {
        this.canvasNormal = new File(normalCanvasFileName);
        this.canvasTumoral = new File(tumoralCanvasFileName);
    }

    private void openOutputStreams(String prefix) throws IOException {
        File normalOutputFile = new File(getBAMOutputName(prefix, NORMAL_DATASET_IDX));
        File tumoralOutputFile = new File(getBAMOutputName(prefix, TUMORAL_DATASET_IDX));
        SAMFileWriterFactory factory = new SAMFileWriterFactory();
        factory.setCreateIndex(true);
        factory.setCompressionLevel(1);
        factory.setUseAsyncIo(true);
        normalWriter = factory.makeBAMWriter(normalFileHeader, true, normalOutputFile);
        tumoralWriter = factory.makeBAMWriter(tumoralFileHeader, true, tumoralOutputFile);
        normalWriter.setSortOrderChecking(false);
        tumoralWriter.setSortOrderChecking(false);
    }

    public String getBAMOutputName(String outputPrefix, int datasetIdx){
        String anonTag = ".anonymized";
        String datasetIdStr = datasetIdx == NORMAL_DATASET_IDX ? ".N" : ".T";
        //String extension = compressed ? ".BAM" : ".SAM";
        String extension = BAM_FILE;
        return outputPrefix + anonTag + datasetIdStr + extension;
    }

    private void mergeAnonymizedReads(File samFile, File canvasSamFile, SAMFileWriter samWriter, File referenceGenomeFile) throws IOException {
        IndexedFastaSequenceFile referenceGenome = new IndexedFastaSequenceFile(referenceGenomeFile);
        SamReader canvasSamReader = null;
        Set<String> canvasToExclude = new HashSet<>();
        if (canvasSamFile != null) {
            canvasSamReader = factory.referenceSequence(referenceGenomeFile).open(canvasSamFile);
            canvasToExclude = extractCanvasReadsToExclude(canvasSamReader, this.genomicRegions);
        }
        Iterator<SAMRecord> canvasReads = canvasSamReader != null ? canvasSamReader.iterator()
                : Collections.emptyIterator();
        SamReader anonymizedSamReader = factory.open(samFile);
        Iterator<SAMRecord> anonymizedReads = getAnonymizedReadsInFile(anonymizedSamReader, referenceGenome);
        SAMRecord lastCanvasRead = null;
        // Write the reads
        while (anonymizedReads.hasNext()) {
            SAMRecord anonymizedRead = anonymizedReads.next();
            int contigIdx = anonymizedRead.getReferenceIndex();
            // If the read contig is "*", set it to the maximum value to write the canvas reads first
            if (contigIdx == -1) {
                contigIdx = Integer.MAX_VALUE;
            }
            // Write the canvas reads that are before the current anonymized read
            while (canvasReads.hasNext() || lastCanvasRead != null) {
                SAMRecord canvasRead = lastCanvasRead != null ? lastCanvasRead : canvasReads.next();
                lastCanvasRead = null;
                if (canvasRead.getReferenceIndex() < contigIdx || (canvasRead.getAlignmentStart() < anonymizedRead.getAlignmentStart() && canvasRead.getReferenceIndex() == contigIdx)) {
                    if (canvasToExclude.contains(canvasRead.getReadName())) {
                        continue;
                    }
                    samWriter.addAlignment(canvasRead);
                } else {
                    lastCanvasRead = canvasRead;
                    break;
                }
            }
            samWriter.addAlignment(anonymizedRead);
        }
        // Iterate through the remaining canvas reads
        while (canvasReads.hasNext() || lastCanvasRead != null) {
            SAMRecord canvasRead = lastCanvasRead != null ? lastCanvasRead : canvasReads.next();
            lastCanvasRead = null;
            if (canvasToExclude.contains(canvasRead.getReadName())) {
                continue;
            }
            samWriter.addAlignment(canvasRead);
        }        
        anonymizedSamReader.close();
        if (canvasSamReader != null) {
            canvasSamReader.close();
        }
        referenceGenome.close();
    }

    private Set<String> extractCanvasReadsToExclude(SamReader canvasSamReader, Iterable<GenomicRegion> regions) {
        Set<String> canvasToExclude = new HashSet<>();
        for (GenomicRegion region : regions) {
            SAMRecordIterator it = canvasSamReader.query(region.getSequenceName(), region.getStart(), region.getEnd(),
                    true);
            while (it.hasNext()) {
                SAMRecord samRecord = it.next();
                canvasToExclude.add(samRecord.getReadName());
            }
            it.close();
        }
        return canvasToExclude;
    }

    private Iterator<SAMRecord> getAnonymizedReadsInFile(SamReader samReader, IndexedFastaSequenceFile referenceGenome) {
        return new Iterator<SAMRecord>() {
            String currentContig = "";
            byte[] referenceContigSequence = new byte[0];
            Iterator<GenomicRegion> partitionIterator = partitions.stream()
                    .sorted(Comparator.comparing(GenomicRegion::getSequenceIdx)
                            .thenComparing(GenomicRegion::getStart))
                    .iterator();
            SAMRecordIterator it = null;
            // Written read alignments to avoid writing again
            Set<String> returnedReads = new HashSet<>();
            SAMRecord samRecordToReturn = getNextSAMRecord();

            @Override
            public boolean hasNext() {
                return samRecordToReturn != null;
            }

            @Override
            public SAMRecord next() {
                if (!hasNext()) {
                    throw new NoSuchElementException();
                }
                SAMRecord samRecord = samRecordToReturn;
                samRecordToReturn = getNextSAMRecord();
                return samRecord;
            }

            private SAMRecord getNextSAMRecord() {
                SAMRecord samRecord = null;
                String readAlnId = "";
                while (samRecord == null) {
                    while (it == null || !it.hasNext()) {
                        if (it != null) {
                            it.close();
                        }
                        if (partitionIterator.hasNext()) {
                            GenomicRegion partition = partitionIterator.next();
                            if (!currentContig.equals(partition.getSequenceName())) {
                                currentContig = partition.getSequenceName();
                                referenceContigSequence = referenceGenome.getSequence(currentContig).getBases();
                            }
                            it = samReader.query(partition.getSequenceName(), partition.getStart(),
                                    partition.getEnd(), false);
                        } else {
                            // No more partitions to process
                            return null;
                        }
                    }
                    samRecord = it.next();
                    String readName = samRecord.getReadName();
                    readAlnId = generateReadId(samRecord);
                    if (readsToExclude.contains(readName) || returnedReads.contains(readAlnId)) {
                        samRecord = null;
                    }
                }
                SAMRecord newSamRecord = samRecord;
                if (readGermlinesToAnonymize.containsKey(readAlnId)) {
                    ShortReadAlignment readAlignment = new ShortReadAlignment(samRecord);
                    ShortAnonymizedReadAlignment anonymizedReadAlignment = new ShortAnonymizedReadAlignment(readAlignment);
                    anonymizedReadAlignment.setReferenceContigSequence(referenceContigSequence);
                    anonymizedReadAlignment.setVariantsToAnonymize(readGermlinesToAnonymize.get(readAlnId));
                    anonymizedReadAlignment.anonymizeVariants();
                    newSamRecord = anonymizedReadAlignment.getAnonymizedSamRecord();
                }
                returnedReads.add(readAlnId);
                return newSamRecord;
            }
        };
    }

    private void closeOutputStreams(){
        normalWriter.close();
        tumoralWriter.close();
    }
}
