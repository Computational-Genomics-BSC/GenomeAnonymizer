package analysis;

import genomicelements.*;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import utils.GlobalRandom;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Future;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.UUID;

import static analysis.GenomeAnonymizer.BAM_FILE;
import static io.SamplePairReadAlignmentReader.DEFAULT_MINIMUM_MAPPING_QUALITY;

/**
 * AnonymizerAlgorithm implementation for short read data
 * @author Nicolas Gaitan
 * @author Rodrigo Martin
 */
public class ShortReadAnonymizer implements AnonymizerAlgorithm {

    private static final Logger LOGGER = Logger.getLogger(ShortReadAnonymizer.class.getName());

    public static final int NORMAL_DATASET_IDX = 0;
    public static final int TUMORAL_DATASET_IDX = 1;

    private List<GenomicRegion> genomicPartitions;
    // Set that contains all the reads that will be excluded from the result (e.g. Unmapped and MAPQ=0)
    private Set<String> readsToExclude;
    // Map containing all potential germlines (value: List), per pair (nested key, 0 or 1), per read (key)
    private Map<String, List<Signal>> readGermlinesToAnonymize;
    private File canvasNormal;
    private File canvasTumoral;
    private List<GenomicRegion> queryRegions;
    private SamReaderFactory factory;
    private SAMFileWriter normalWriter;
    private SAMFileWriter tumoralWriter;
    private int hashSalt;

    public ShortReadAnonymizer() {
        readsToExclude = new HashSet<>();
        readGermlinesToAnonymize = new HashMap<>();
        factory = SamReaderFactory.makeDefault();
        factory.setUseAsyncIo(true);
        factory.validationStringency(ValidationStringency.SILENT);
        queryRegions = new ArrayList<>();
        // Initialize hash_salt with a random value
        hashSalt = GlobalRandom.getInstance().nextInt();
    }

    public void setGenomicPartitions(List<GenomicRegion> genomicPartitions) {
        this.genomicPartitions = genomicPartitions;
    }

    public void setReadGermlinesToAnonymize(Map<String, List<Signal>> readGermlinesToAnonymize) {
        this.readGermlinesToAnonymize = readGermlinesToAnonymize;
    }

    private SAMFileHeader buildFileHeader(String bamFile, String sampleSuffix) throws IOException {
        try (SamReader samReader = factory.open(new File(bamFile))) {
            // Retrieve the SAMFileHeader
            SAMFileHeader header = samReader.getFileHeader();
            // Remove all read groups from the header
            // Set the read group as the hash of the file name + salt
            String readGroupId = Integer.toHexString((bamFile + hashSalt).hashCode());
            SAMReadGroupRecord readGroup = new SAMReadGroupRecord(readGroupId);
            readGroup.setSample(readGroupId + sampleSuffix);
            header.setReadGroups(Collections.singletonList(readGroup));
            return header;
        }
    }

    public void queryReadsToExclude(String normalPath, String tumorPath, int threads) throws Exception{
        String[] paths = new String[2];
        paths[0] = normalPath;
        paths[1] = tumorPath;
        ExecutorService exec = Executors.newFixedThreadPool(threads);
        List<CompletableFuture<Set<String>>> futures = new ArrayList<>();
        //TODO: Check if it is possible to change partitions based on actual content
        // (Implement CoveredGenomicRegion), to improve runtime using parallelization
        for(GenomicRegion partition : genomicPartitions){
            for(String path : paths){
                CompletableFuture<Set<String>> future = CompletableFuture.supplyAsync (() -> {
                    Set<String> answer;
                    try {
                        answer = queryReadsToExcludeInPartition(path, partition);
                    }
                    catch (IOException e) {
                        LOGGER.log(Level.SEVERE,
                                "Exception in thread querying reads to exclude in region: "
                                + partition.getSequenceName()
                                + " " + partition.getStart() + " " + partition.getEnd(),
                                e);
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
            LOGGER.log(Level.SEVERE,
                    "Exception when retrieving excluded reads from thread, halting execution prematurely",
                    e);
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
                        (samRecord.getMappingQuality() < DEFAULT_MINIMUM_MAPPING_QUALITY && !samRecord.isSecondaryOrSupplementary())){
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
            openOutputStreams(outputPrefix, normalPath, tumorPath);
            ExecutorService exec = Executors.newFixedThreadPool(2);
            Future<?> tumoralFuture = exec.submit(() -> {
                try {
                    mergeAnonymizedReads(new File(tumorPath), this.canvasTumoral, this.tumoralWriter, new File(refGenome),
                            false);
                } catch (Exception e) {
                    LOGGER.log(Level.SEVERE,
                            "Exception in thread merging anonymized reads in tumoral dataset", e);
                    throw new RuntimeException(e);
                }
            });
            Future<?> normalFuture = exec.submit(() -> {
                try {
                    mergeAnonymizedReads(new File(normalPath), this.canvasNormal, this.normalWriter, new File(refGenome),
                            true);
                } catch (Exception e) {
                    LOGGER.log(Level.SEVERE,
                            "Exception in thread merging anonymized reads in normal dataset", e);
                    throw new RuntimeException(e);
                }
            });
            exec.shutdown();
            normalFuture.get();
            tumoralFuture.get();
        } catch (Exception e) {
            LOGGER.log(Level.SEVERE,
                    "IOException in reading/writing anonymized reads, halting execution prematurely", e);
            throw new IOException(e);
        } finally {
            closeOutputStreams();
        }
    }

    @Override
    public void setQueryRegions(List<GenomicRegion> regions){
        this.queryRegions = regions;
    }

    @Override
    public void setCanvasFiles(String normalCanvasFileName, String tumoralCanvasFileName) {
        this.canvasNormal = new File(normalCanvasFileName);
        this.canvasTumoral = new File(tumoralCanvasFileName);
    }

    private void openOutputStreams(String prefix, String normalPath, String tumorPath) throws IOException {
        File normalOutputFile = new File(getBAMOutputName(prefix, NORMAL_DATASET_IDX));
        File tumoralOutputFile = new File(getBAMOutputName(prefix, TUMORAL_DATASET_IDX));
        SAMFileWriterFactory factory = new SAMFileWriterFactory();
        //factory.setCreateIndex(true);
        factory.setCompressionLevel(1);
        factory.setUseAsyncIo(true);
        SAMFileHeader normalFileHeader = buildFileHeader(normalPath, "_N");
        SAMFileHeader tumoralFileHeader = buildFileHeader(tumorPath, "_T");
        //normalWriter = factory.makeBAMWriter(normalFileHeader, false, normalOutputFile);
        //tumoralWriter = factory.makeBAMWriter(tumoralFileHeader, false, tumoralOutputFile);
        normalWriter = factory.makeBAMWriter(normalFileHeader, true, normalOutputFile);
        tumoralWriter = factory.makeBAMWriter(tumoralFileHeader, true, tumoralOutputFile);
        normalWriter.setSortOrderChecking(false);
        tumoralWriter.setSortOrderChecking(false);
        LOGGER.info("Beginning writing Anonymized reads to\tNormal: " + normalOutputFile + "\tTumoral: " +
                tumoralOutputFile);
    }

    public String getBAMOutputName(String outputPrefix, int datasetIdx){
        String anonTag = ".anonymized";
        String datasetIdStr = datasetIdx == NORMAL_DATASET_IDX ? ".N" : ".T";
        String extension = BAM_FILE;
        return outputPrefix + anonTag + datasetIdStr + extension;
    }

    private void mergeAnonymizedReads(File samFile, File canvasSamFile, SAMFileWriter samWriter, File referenceGenomeFile,
                                      boolean isNormalDataset) throws IOException {
        IndexedFastaSequenceFile referenceGenome = new IndexedFastaSequenceFile(referenceGenomeFile);
        SamReader canvasSamReader = null;
        Set<String> canvasToExclude = new HashSet<>();
        int writtenReads = 0;
        String datasetName = isNormalDataset ? "Normal" : "Tumoral";
        if (canvasSamFile != null) {
            canvasSamReader = factory.referenceSequence(referenceGenomeFile).open(canvasSamFile);
            canvasToExclude = extractCanvasReadsToExclude(canvasSamReader, this.queryRegions);
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
                    writeRead(samWriter, canvasRead);
                    writtenReads++;
                    logWrittenReads(writtenReads, datasetName);
                } else {
                    lastCanvasRead = canvasRead;
                    break;
                }
            }
            writeRead(samWriter, anonymizedRead);
            writtenReads++;
            logWrittenReads(writtenReads, datasetName);
        }
        // Iterate through the remaining canvas reads
        while (canvasReads.hasNext() || lastCanvasRead != null) {
            SAMRecord canvasRead = lastCanvasRead != null ? lastCanvasRead : canvasReads.next();
            lastCanvasRead = null;
            if (canvasToExclude.contains(canvasRead.getReadName())) {
                continue;
            }
            writeRead(samWriter, canvasRead);
            writtenReads++;
            logWrittenReads(writtenReads, datasetName);
        }        
        anonymizedSamReader.close();
        if (canvasSamReader != null) {
            canvasSamReader.close();
        }
        referenceGenome.close();
    }

    private void logWrittenReads(int writtenReads, String datasetName) {
        if (writtenReads % 10_000_000 == 0) LOGGER.info("Written " + writtenReads + " to " + datasetName);
    }

    private Set<String> extractCanvasReadsToExclude(SamReader canvasSamReader, Iterable<GenomicRegion> regions) {
        Set<String> canvasToExclude = new HashSet<>();
        for (GenomicRegion region : regions) {
            SAMRecordIterator it = canvasSamReader.query(region.getSequenceName(), region.getStart(), region.getEnd(),
                    false);
            while (it.hasNext()) {
                SAMRecord samRecord = it.next();
                canvasToExclude.add(samRecord.getReadName());
            }
            it.close();
        }
        return canvasToExclude;
    }

    private void writeRead(SAMFileWriter samWriter, SAMRecord read) {
        // Set the read name to the hash of the read name + salt
        read.setReadName(UUID.nameUUIDFromBytes((read.getReadName() + hashSalt).getBytes()).toString());
        // Set the read group to the same value as the read group of the header
        read.setAttribute("RG", samWriter.getFileHeader().getReadGroups().get(0).getId());
        samWriter.addAlignment(read);
    }

    private Iterator<SAMRecord> getAnonymizedReadsInFile(SamReader samReader, IndexedFastaSequenceFile referenceGenome) {
        return new Iterator<SAMRecord>() {
            String currentContig = "";
            byte[] referenceContigSequence = new byte[0];
            Iterator<GenomicRegion> partitionIterator = genomicPartitions.stream()
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
                    readAlnId = ShortReadAlignment.generateReadAlnId(samRecord);
                    if (readsToExclude.contains(readName) || returnedReads.contains(readAlnId)) {
                        samRecord = null;
                    }
                }
                SAMRecord newSamRecord = samRecord;
                if (readGermlinesToAnonymize.containsKey(readAlnId)) {
                    ShortReadAlignment readAlignment = new ShortReadAlignment(samRecord);
                    ShortAnonymizedReadAlignment anonymizedReadAlignment = new ShortAnonymizedReadAlignment(readAlignment);
                    anonymizedReadAlignment.setReferenceContigSequence(referenceContigSequence);
                    anonymizedReadAlignment.setSignalsToAnonymize(readGermlinesToAnonymize.get(readAlnId));
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
