package analysis;

import genomicelements.*;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import utils.GlobalRandom;
import utils.Tuple;

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
import static genomicelements.ShortReadAlignment.generateReadId;

/**
 * AnonymizerAlgorithm implementation for short read data
 * @author Nicolas Gaitan
 * @author Rodrigo Martin
 */
public class ShortReadAnonymizer implements AnonymizerAlgorithm {

    private static final Logger LOGGER = Logger.getLogger(ShortReadAnonymizer.class.getName());

    public static final int NORMAL_DATASET_IDX = 0;
    public static final int TUMORAL_DATASET_IDX = 1;
    private static final int INSERT_SIZE_BIN_SIZE = 50;
    private static final int INSERT_SIZE_BIN_COUNT = 5000 / INSERT_SIZE_BIN_SIZE + 1;

    private List<GenomicRegion> genomicPartitions;
    // Set that contains all the reads that will be excluded from the result (e.g. Unmapped and MAPQ=0)
    private Set<String> readsToExclude;
    // Thresholds for the insert sizes
    private int insertSizeMinThreshold = Integer.MIN_VALUE;
    private int insertSizeMaxThreshold = Integer.MAX_VALUE;
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
        List<CompletableFuture<Tuple<Set<String>, List<Integer>>>> futures = new ArrayList<>();
        //TODO: Check if it is possible to change partitions based on actual content
        // (Implement CoveredGenomicRegion), to improve runtime using parallelization
        for(GenomicRegion partition : genomicPartitions){
            for(String path : paths){
                CompletableFuture<Tuple<Set<String>, List<Integer>>> future = CompletableFuture.supplyAsync (() -> {
                    Tuple<Set<String>, List<Integer>> answer;
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
        List<Integer> insertSizesBins = new ArrayList<>(Collections.nCopies(INSERT_SIZE_BIN_COUNT, 0));
        CompletableFuture<Void> allFutures = CompletableFuture.allOf(futures.toArray(new CompletableFuture[0]));
        try{
            for (CompletableFuture<Tuple<Set<String>, List<Integer>>> future : futures) {
                Set<String> answer = future.get().getFirst();
                readsToExclude.addAll(answer);
                List<Integer> insertSizesBinsInPartition = future.get().getSecond();
                // Sum each element of the list with the corresponding element of the other list
                for (int i = 0; i < insertSizesBins.size(); i++) {
                    insertSizesBins.set(i, insertSizesBins.get(i) + insertSizesBinsInPartition.get(i));
                }
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
        // Calculate the bottom 5% and top 95% quantiles of the insert sizes aproximated by the corresponding bins
        int totalReads = insertSizesBins.stream().mapToInt(Integer::intValue).sum();
        int bottom5Percentile = (int) Math.ceil(totalReads * 0.05);
        int top95Percentile = (int) Math.ceil(totalReads * 0.95);
        int currentReads = 0;
        for (int i = 0; i < insertSizesBins.size(); i++) {
            currentReads += insertSizesBins.get(i);
            if (currentReads >= bottom5Percentile && insertSizeMinThreshold == Integer.MIN_VALUE) {
            int previousReads = currentReads - insertSizesBins.get(i);
            double fraction = (bottom5Percentile - previousReads) / (double) insertSizesBins.get(i);
            insertSizeMinThreshold = (int) ((i - 1 + fraction) * INSERT_SIZE_BIN_SIZE);
            }
            if (currentReads >= top95Percentile && insertSizeMaxThreshold == Integer.MAX_VALUE) {
            int previousReads = currentReads - insertSizesBins.get(i);
            double fraction = (top95Percentile - previousReads) / (double) insertSizesBins.get(i);
            insertSizeMaxThreshold = (int) ((i - 1 + fraction) * INSERT_SIZE_BIN_SIZE);
            break;
            }
        }
    }

    private Tuple<Set<String>, List<Integer>> queryReadsToExcludeInPartition(String filePath, GenomicRegion partition) throws IOException {
        Set<String> readsToExcludeInPartition = new HashSet<>();
        List<Integer> insertSizesBinsInPartition = new ArrayList<>(Collections.nCopies(INSERT_SIZE_BIN_COUNT, 0));
        try (SamReader reader = factory.open(new File(filePath))) {
            SAMRecordIterator it = reader.query(partition.getSequenceName(), partition.getStart(), partition.getEnd(), false);
            while (it.hasNext()) {
                SAMRecord samRecord = it.next();
                if (samRecord.getReadUnmappedFlag() ||
                        (samRecord.getMappingQuality() == 0 && !samRecord.isSecondaryOrSupplementary())) {
                    readsToExcludeInPartition.add(samRecord.getReadName());
                    continue;
                }
                // Only use the first read
                if (!samRecord.getReadPairedFlag() || samRecord.getMateUnmappedFlag() || !samRecord.getFirstOfPairFlag()) {
                    continue;
                }
                // Get the insert sizes
                int insertSize = samRecord.getInferredInsertSize();
                if (insertSize <= 0 || samRecord.getReadNegativeStrandFlag() || !samRecord.getMateNegativeStrandFlag()) {
                    continue;
                }
                // Add to the corresponding bucket
                int binIdx = Math.min(insertSize / INSERT_SIZE_BIN_SIZE, INSERT_SIZE_BIN_COUNT - 1);
                insertSizesBinsInPartition.set(binIdx, insertSizesBinsInPartition.get(binIdx) + 1);
            }
        }
        return new Tuple<>(readsToExcludeInPartition, insertSizesBinsInPartition);
    }

    @Override
    public int getInsertSizeMinThreshold() {
        return insertSizeMinThreshold;
    }

    @Override
    public int getInsertSizeMaxThreshold() {
        return insertSizeMaxThreshold;
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
                    mergeAnonymizedReads(new File(tumorPath), this.canvasTumoral, this.tumoralWriter, new File(refGenome));
                } catch (Exception e) {
                    LOGGER.log(Level.SEVERE,
                            "Exception in thread merging anonymized reads in tumoral dataset", e);
                    throw new RuntimeException(e);
                }
            });
            Future<?> normalFuture = exec.submit(() -> {
                try {
                    mergeAnonymizedReads(new File(normalPath), this.canvasNormal, this.normalWriter, new File(refGenome));
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
    }

    public String getBAMOutputName(String outputPrefix, int datasetIdx){
        String anonTag = ".anonymized";
        String datasetIdStr = datasetIdx == NORMAL_DATASET_IDX ? ".N" : ".T";
        String extension = BAM_FILE;
        return outputPrefix + anonTag + datasetIdStr + extension;
    }

    private void mergeAnonymizedReads(File samFile, File canvasSamFile, SAMFileWriter samWriter, File referenceGenomeFile) throws IOException {
        IndexedFastaSequenceFile referenceGenome = new IndexedFastaSequenceFile(referenceGenomeFile);
        SamReader canvasSamReader = null;
        Set<String> canvasToExclude = new HashSet<>();
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
                } else {
                    lastCanvasRead = canvasRead;
                    break;
                }
            }
            writeRead(samWriter, anonymizedRead);
        }
        // Iterate through the remaining canvas reads
        while (canvasReads.hasNext() || lastCanvasRead != null) {
            SAMRecord canvasRead = lastCanvasRead != null ? lastCanvasRead : canvasReads.next();
            lastCanvasRead = null;
            if (canvasToExclude.contains(canvasRead.getReadName())) {
                continue;
            }
            writeRead(samWriter, canvasRead);
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
