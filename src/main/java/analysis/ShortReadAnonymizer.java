package analysis;

import genomicelements.*;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.FastaSequenceIndexEntry;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.IOUtil;
import utils.GlobalRandom;
import utils.Tuple;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import static analysis.GenomeAnonymizer.*;
import static genomicelements.ShortAnonymizedReadAlignment.*;

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
    private static final float DEFAULT_INSERT_SIZE_THRESHOLD_FRACTION = 0.05f;

    private String inputNormalPath;
    private String inputTumorPath;
    private String refGenomePath;
    private String outputPrefix;

    private List<GenomicRegion> genomicPartitions = new ArrayList<>();
    private Map<String, byte[]> referenceSequences = new HashMap<>();

    // Set that contains all the reads that will be excluded from the result (e.g. Unmapped and MAPQ < filter)
    private Set<String> readsToExclude;
    private Map<String, Integer> updatedMatePositions;
    private Map<String, Integer> readsToCorrectOrientation;

    // Thresholds for the insert sizes
    private int insertSizeMinThreshold = Integer.MIN_VALUE;
    private int insertSizeMedian = Integer.MIN_VALUE;
    private int insertSizeMaxThreshold = Integer.MAX_VALUE;
    private float insertSizeThresholdFraction = DEFAULT_INSERT_SIZE_THRESHOLD_FRACTION;
    private File tmpDir;
    private SamReaderFactory factory;

    private int minimumMappingQuality = DEFAULT_MIN_MAPPING_QUALITY;
    private int hashSalt;
    private int threads = 1;
    private int maxReadsInRam = DEFAULT_MAX_READS_IN_RAM;

    //TODO: Save partitioned files names sorted by region coordinates

    public ShortReadAnonymizer(String inputNormalPath, String inputTumorPath, String refGenomePath, String outputPrefix) {
        this.inputNormalPath = inputNormalPath;
        this.inputTumorPath = inputTumorPath;
        this.refGenomePath = refGenomePath;
        this.outputPrefix = outputPrefix;
        this.readsToExclude = new HashSet<>();
        this.updatedMatePositions = new HashMap<>();
        this.readsToCorrectOrientation = new HashMap<>();
        this.factory = SamReaderFactory.makeDefault();
        this.factory.setUseAsyncIo(true);
        this.factory.validationStringency(ValidationStringency.SILENT);
        // Initialize hash_salt with a random value
        this.hashSalt = GlobalRandom.getInstance().nextInt();
    }

    @Override
    public void setMinimumMappingQuality(int minimumMappingQuality) {
        this.minimumMappingQuality = minimumMappingQuality;
    }

    public void setTmpDir(File tmpDir) {
        // This is necessary to avoid all reads being kept in memory
        if (!tmpDir.exists()) tmpDir.mkdirs();
        tmpDir.setReadable(true, false);
        tmpDir.setWritable(true, false);
        System.setProperty("java.io.tmpdir", tmpDir.getAbsolutePath());
        this.tmpDir = tmpDir;
    }

    public void setThreadNumber(int threads) {
        this.threads = threads;
    }

    public void setMaxReadsInRam(int maxReadsInRam) {
        this.maxReadsInRam = maxReadsInRam;
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

    public void setPartitions() throws IOException {
        IndexedFastaSequenceFile reference = new IndexedFastaSequenceFile(new File(refGenomePath));
        FastaSequenceIndex refIndexes = reference.getIndex();
        long genomeSize = 0;
        List<FastaSequenceIndexEntry> sequences = new ArrayList<>();
        Map<String, Integer> refSequenceOrder = new HashMap<>();
        int idx = 0;
        for (FastaSequenceIndexEntry refEntry : refIndexes){
            String contig = refEntry.getContig();
            referenceSequences.put(contig, reference.getSequence(contig).getBases());
            sequences.add(refEntry);
            genomeSize += refEntry.getSize();
            refSequenceOrder.put(contig, idx);
            idx++;
        }
        reference.close();
        // Compute genome partitions
        sequences.sort(Comparator.comparing(FastaSequenceIndexEntry::getSize));
        Collections.reverse(sequences);
        int[] partitionsPerSequence = new int[sequences.size()];
        Arrays.fill(partitionsPerSequence, 1);
        long basesPerThread = genomeSize / (threads * 5L);
        // Estimate threads to be assigned to each contig
        for(int i = 0; i < sequences.size(); i++){
            partitionsPerSequence[i] += (int) (sequences.get(i).getSize() / basesPerThread);
        }
        // Make propper contig partitions into the genomicRegions, based on n assigned threads
        for(int i = 0; i < sequences.size(); i++){
            FastaSequenceIndexEntry currentContig = sequences.get(i);
            String contig = currentContig.getContig();
            int contigIdx = refSequenceOrder.get(contig);
            int contigLength = (int) currentContig.getSize();
            int partitionSize = contigLength / partitionsPerSequence[i];
            int currentFirst = 1;
            for(int j = 0; j < partitionsPerSequence[i]; j++){
                //Be careful with very large chromosomes, with humans there should not be a problem
                int currentLast = j == partitionsPerSequence[i]-1 ? (int) currentContig.getSize() : currentFirst + partitionSize;
                GenomicRegion region = new GenomicRegionBaseImpl(contig, contigIdx, currentFirst, currentLast);
                genomicPartitions.add(region);
                currentFirst += partitionSize + 1;
            }
        }
    }

    public void queryReadsToExclude() {
        String[] paths = new String[2];
        paths[0] = inputNormalPath;
        paths[1] = inputTumorPath;
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
        CompletableFuture<Void> allFutures = CompletableFuture.allOf(futures.toArray(new CompletableFuture[0]));
        List<Integer> insertSizesBins = new ArrayList<>(Collections.nCopies(INSERT_SIZE_BIN_COUNT, 0));
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
        // Calculate the bottom % and top % quantiles of the insert sizes aproximated by the corresponding bins
        int totalReads = insertSizesBins.stream().mapToInt(Integer::intValue).sum();
        int bottomPercentile = (int) Math.ceil(totalReads * insertSizeThresholdFraction);
        int medianPercentile = (int) Math.ceil(totalReads * 0.5);
        int topPercentile = (int) Math.ceil(totalReads * (1 - insertSizeThresholdFraction));
        int currentReads = 0;
        for (int i = 0; i < insertSizesBins.size(); i++) {
            currentReads += insertSizesBins.get(i);
            if (currentReads >= bottomPercentile && insertSizeMinThreshold == Integer.MIN_VALUE) {
                int previousReads = currentReads - insertSizesBins.get(i);
                double fraction = (bottomPercentile - previousReads) / (double) insertSizesBins.get(i);
                insertSizeMinThreshold = (int) ((i - 1 + fraction) * INSERT_SIZE_BIN_SIZE);
            }
            if(currentReads >= medianPercentile && insertSizeMedian == Integer.MIN_VALUE) {
                int previousReads = currentReads - insertSizesBins.get(i);
                double fraction = (medianPercentile - previousReads) / (double) insertSizesBins.get(i);
                insertSizeMedian = (int) ((i - 1 + fraction) * INSERT_SIZE_BIN_SIZE);
            }
            if (currentReads >= topPercentile && insertSizeMaxThreshold == Integer.MAX_VALUE) {
                int previousReads = currentReads - insertSizesBins.get(i);
                double fraction = (topPercentile - previousReads) / (double) insertSizesBins.get(i);
                insertSizeMaxThreshold = (int) ((i - 1 + fraction) * INSERT_SIZE_BIN_SIZE);
                break;
            }
        }
    }

    private Tuple<Set<String>, List<Integer>> queryReadsToExcludeInPartition(String filePath, GenomicRegion partition) throws IOException {
        Set<String> readsToExcludeInPartition = new HashSet<>();
        List<Integer> insertSizesBinsInPartition = new ArrayList<>(Collections.nCopies(INSERT_SIZE_BIN_COUNT, 0));
        try(SamReader reader = factory.open(new File(filePath))){
            SAMRecordIterator it = reader.query(partition.getSequenceName(), partition.getStart(), partition.getEnd(), false);
            while (it.hasNext()) {
                SAMRecord samRecord = it.next();
                if((samRecord.getReadUnmappedFlag() || samRecord.getMateUnmappedFlag()) ||
                        (samRecord.getMappingQuality() < minimumMappingQuality && !samRecord.isSecondaryOrSupplementary())){
                    readsToExcludeInPartition.add(samRecord.getReadName());
                    continue;
                }
                // Only use the first read
                if (!samRecord.getReadPairedFlag() || !samRecord.getFirstOfPairFlag()) {
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
    public void anonymizeReads() {
        // If this.tmpDir is not set, set it to the default temporary directory
        if (this.tmpDir == null) {
            setTmpDir(IOUtil.getDefaultTmpDir());
        }
        List<String> normalPaths = new ArrayList<>();
        List<String> tumorPaths = new ArrayList<>();
        PartitionProviderRunner[] partitionRunnables = new PartitionProviderRunner[genomicPartitions.size()];
        for(int i = 0; i < partitionRunnables.length; i++){
            GenomicRegion genomicPartition = genomicPartitions.get(i);
            String suffix = "_" + genomicPartition.toString();
            String normalOutputPath = getBAMOutputName(outputPrefix+suffix, NORMAL_DATASET_IDX);
            String tumorOutputPath = getBAMOutputName(outputPrefix+suffix, TUMORAL_DATASET_IDX);
            partitionRunnables[i] = new PartitionProviderRunner(genomicPartition, normalOutputPath, tumorOutputPath);
            normalPaths.add(normalOutputPath);
            tumorPaths.add(tumorOutputPath);
        }
        ExecutorService executorService = Executors.newFixedThreadPool(threads);
        List<CompletableFuture<Void>> futures = new ArrayList<>();
        for(PartitionProviderRunner runnable : partitionRunnables){
            CompletableFuture<Void> future = CompletableFuture.runAsync(runnable, executorService);
            futures.add(future);
        }
        CompletableFuture.allOf(futures.toArray(new CompletableFuture[0])).join();
        executorService.shutdown();
        // Collect the answers from all the threads
        for(PartitionProviderRunner runnable : partitionRunnables){
            PartitionAnswer answer = runnable.getAnswer();
            // Get all pairs to update with info from the mate (indexed by the  name of the pair to be updated)
            updatedMatePositions.putAll(answer.partitionPairsToUpdate());
            // Get all reads to correct orientation
            for (String readName : answer.readsToCorrectOrientation()) {
                readsToCorrectOrientation.putIfAbsent(readName, 0);
            }
        }
        // Free memory cleaning unused attributes
        readsToExclude = new HashSet<>();
        // Merge the anonymized reads from all partition files
        mergeAnonymizedReads(normalPaths, tumorPaths);
        // Delete the temporary files
        deleteTempFiles(normalPaths);
        deleteTempFiles(tumorPaths);
    }

    private void deleteTempFiles(List<String> paths) {
        for (String path : paths) {
            File file = new File(path);
            if (file.exists()) {
                try {
                    file.delete();
                } catch (SecurityException e) {
                    LOGGER.log(Level.WARNING, "Unable to delete file: " + file.getAbsolutePath(), e);
                }
            }
        }
    }

    public class PartitionProviderRunner implements Runnable {
        private final GenomicRegion genomicPartition;
        private final String normalOutputPath;
        private final String tumorOutputPath;

        private Set<String> partitionReadsToCorrectOrientation = new HashSet<>();
        private Map<String, Integer> partitionPairsToUpdate = new HashMap<>();

        private PartitionAnswer answer;

        public PartitionProviderRunner(GenomicRegion genomicPartition, String normalOutputPath, String tumorOutputPath) {
            this.genomicPartition = genomicPartition;
            this.normalOutputPath = normalOutputPath;
            this.tumorOutputPath = tumorOutputPath;
        }

        public PartitionAnswer getAnswer() {
            return answer;
        }

        @Override
        public void run() {
            try (SAMFileWriter partitionNormalWriter = openSingleOutputStream(inputNormalPath, normalOutputPath, true);
                 SAMFileWriter partitionTumoralWriter = openSingleOutputStream(inputTumorPath, tumorOutputPath, false);
                 AnonymizedReadAlignmentProvider anonymizedReadProvider = new AnonymizedReadAlignmentProvider()) {
                anonymizedReadProvider.setReadsToExclude(readsToExclude);
                anonymizedReadProvider.setInsertSizeMinThreshold(insertSizeMinThreshold);
                anonymizedReadProvider.setInsertSizeMaxThreshold(insertSizeMaxThreshold);
                anonymizedReadProvider.setRefSequence(referenceSequences.get(genomicPartition.getSequenceName()));
                anonymizedReadProvider.init(inputNormalPath, inputTumorPath, refGenomePath, genomicPartition);
                long startcallVariation = System.currentTimeMillis();
                for (AnonymizedRead anonymizedRead : anonymizedReadProvider) {
                    ShortAnonymizedReadAlignment shortAnonymizedRead = (ShortAnonymizedReadAlignment) anonymizedRead;
                    SAMFileWriter writer = shortAnonymizedRead.isFromNormalDataset() ? partitionNormalWriter : partitionTumoralWriter;
                    int pairIdx = shortAnonymizedRead.getPairIdx();
                    String readName = shortAnonymizedRead.getReadName();
                    String thisPairName = getPairedReadName(readName, pairIdx);
                    String otherPairName = getPairedReadName(readName, 1 - pairIdx);
                    //If the first pair contained a signal mandating orientation to be fixed, fix the second pair orientation also
                    if(partitionReadsToCorrectOrientation.contains(readName) && !shortAnonymizedRead.isSupplementary()) {
                        shortAnonymizedRead.setFixOrientation();
                        partitionReadsToCorrectOrientation.remove(readName);
                    }
                    //Anonymize read if it has germline signals
                    if(!shortAnonymizedRead.isAnonymized()) shortAnonymizedRead.anonymizeRead();
                    //If the first pair contains a signal mandating orientation to be fixed, save the read name to fix the second pair orientation
                    if(shortAnonymizedRead.fixOrientation()) partitionReadsToCorrectOrientation.add(readName);
                    // Update info from the other mate into this pair
                    if (partitionPairsToUpdate.containsKey(thisPairName)){
                        int updatedMatePos = partitionPairsToUpdate.get(thisPairName);
                        shortAnonymizedRead.getAnonymizedSamRecord().setMateAlignmentStart(updatedMatePos);
                        partitionPairsToUpdate.remove(thisPairName);
                    }
                    // Update info from this pair into the other mate
                    if (shortAnonymizedRead.updateInfoForMate() && !shortAnonymizedRead.isSupplementary()){
                        int thisMateUpdatedPos = shortAnonymizedRead.getPairUpdatedPos();
                        partitionPairsToUpdate.put(otherPairName, thisMateUpdatedPos);
                    }
                    if(shortAnonymizedRead.hasDestructiveSignal() && !shortAnonymizedRead.isSupplementary()){
                        partitionPairsToUpdate.remove(thisPairName);
                        if(shortAnonymizedRead.mateOriginalPosIsEqual()) partitionPairsToUpdate.put(otherPairName, -1);
                        else partitionPairsToUpdate.remove(otherPairName);
                        SAMRecord newPairAlnRecord = shortAnonymizedRead.getNewPair(insertSizeMedian);
                        writeRead(writer, newPairAlnRecord);
                    }
                    // With all updated information, get this anonymized sam record
                    SAMRecord alnRecord = shortAnonymizedRead.getAnonymizedSamRecord();
                    long startWriteRead = System.currentTimeMillis();
                    writeRead(writer, alnRecord);
                    long endWriteRead = System.currentTimeMillis();
                    long writeReadTime = endWriteRead - startWriteRead;
                    anonymizedReadProvider.METHOD_TIME_MAP.compute("writeReadTime",  (k,v) -> v == null ?
                            writeReadTime :
                            v + writeReadTime);
                }
                answer = new PartitionAnswer(partitionPairsToUpdate, partitionReadsToCorrectOrientation);
                long endcallVariation = System.currentTimeMillis();
                //TIME DEBUG
                anonymizedReadProvider.METHOD_TIME_MAP.put("callVariation", endcallVariation - startcallVariation);
                LOGGER.info("Finished variation analysis of genomic region: SEQ=" + genomicPartition.getSequenceName()
                        + " POS=" + genomicPartition.getStart() + " END=" + genomicPartition.getEnd());
                StringBuilder msg = new StringBuilder();
                for (Map.Entry<String, Long> entry : anonymizedReadProvider.METHOD_TIME_MAP.entrySet()) {
                    msg.append(entry.getKey()).append(": ");
                    long timeInSeconds = entry.getValue();
                    msg.append(timeInSeconds).append("\n");
                }
                LOGGER.info("PARTITION TIME TABLE: " + "\n" + msg);
                //TIME DEBUG
            }
            catch (Exception e) {
                LOGGER.log(Level.SEVERE,
                        "Exception in thread Anonymizing reads in region: "
                                + genomicPartition.getSequenceName()
                                + " " + genomicPartition.getStart() + " " + genomicPartition.getEnd(),
                        e);
                System.exit(1);
//                    throw new RuntimeException(e);
            }
        }
    }

    public record PartitionAnswer(Map<String, Integer> partitionPairsToUpdate, Set<String> readsToCorrectOrientation) { }

    private SAMFileWriter openSingleOutputStream(String path, String outputPath, boolean isNormalDataset) throws IOException {
        SAMFileWriterFactory factory = new SAMFileWriterFactory();
        factory.setCompressionLevel(1);
        factory.setMaxRecordsInRam(maxReadsInRam / threads);
        SAMFileHeader fileHeader = buildFileHeader(path, isNormalDataset ? "_N" : "_T");
        fileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return factory.makeWriter(fileHeader, false, new File(outputPath), new File(refGenomePath));
    }

    public String getBAMOutputName(String outputPrefix, int datasetIdx){
        String anonTag = ".anonymized";
        String datasetIdStr = datasetIdx == NORMAL_DATASET_IDX ? ".N" : ".T";
        return outputPrefix + anonTag + datasetIdStr + BAM_FILE;
    }

    private void writeRead(SAMFileWriter samWriter, SAMRecord read) {
        // Set the read group to the same value as the read group of the header
        read.setAttribute("RG", samWriter.getFileHeader().getReadGroups().get(0).getId());
        samWriter.addAlignment(read);
    }

    public String getPairedReadName(String readName, int pairIdx) {
        final StringBuilder builder = new StringBuilder(64);
        builder.append(readName);
        if (pairIdx == 0) {
            builder.append(" 1/2");
        }
        else {
            builder.append(" 2/2");
        }
        return builder.toString();
    }

    private void logWrittenReads(int writtenReads, String datasetName) {
        if (writtenReads % 10_000_000 == 0) LOGGER.info("Written " + writtenReads + " to " + datasetName);
    }

    public void mergeAnonymizedReads(List<String> normalPaths, List<String> tumorPaths) {
        String outputNormal = outputPrefix + ".anonymized.N" + BAM_FILE;
        String outputTumor = outputPrefix + ".anonymized.T" + BAM_FILE;
        LOGGER.log(Level.INFO, "Merging normal to " + outputNormal);
        LOGGER.log(Level.INFO, "Merging tumoral to " + outputTumor);
        try {
            ExecutorService exec = Executors.newFixedThreadPool(2);
            Future<?> tumoralFuture = exec.submit(() -> {
                try {
                    mergeReads(outputTumor, tumorPaths);
                } catch (Exception e) {
                    LOGGER.log(Level.SEVERE,
                            "Exception in thread merging anonymized reads in tumoral dataset", e);
                    throw new RuntimeException(e);
                }
            });
            Future<?> normalFuture = exec.submit(() -> {
                try {
                    mergeReads(outputNormal, normalPaths);
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
            throw new RuntimeException(e);
        }
    }

    private void mergeReads(String outputPath, List<String> partitionPaths) throws IOException {
        // Open the files for reading and writing
        final List<SamReader> readers = new ArrayList<>();
        final List<SAMFileHeader> headers = new ArrayList<>();
        for (final String partitionPath : partitionPaths) {
            final SamReader reader = factory.open(new File(partitionPath));
            readers.add(reader);
            headers.add(reader.getFileHeader());
        }
        // Merge the headers
        final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SAMFileHeader.SortOrder.coordinate, headers, false);
        SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
        writerFactory.setCompressionLevel(1);
        writerFactory.setMaxRecordsInRam(maxReadsInRam / 2);
        writerFactory.setUseAsyncIo(true);
        writerFactory.setCreateIndex(true);
        final SAMFileWriter writer = writerFactory.makeWriter(headerMerger.getMergedHeader(), true,
                new File(outputPath), new File(refGenomePath));
        // Merge the records
        final MergingSamRecordIterator iterator = new MergingSamRecordIterator(headerMerger, readers, false);
        Map<String, Integer> generatedReadPairs = new HashMap<>();
        while (iterator.hasNext()) {
            final SAMRecord record = iterator.next();
            boolean keepPair = keepPair(record, generatedReadPairs);
            if(keepPair) {
                // Set the read name to the hash of the read name + salt
//            record.setReadName(UUID.nameUUIDFromBytes((record.getReadName() + hashSalt).getBytes()).toString());
                if(readsToCorrectOrientation.containsKey(record.getReadName())){
                    int pairCount = readsToCorrectOrientation.get(record.getReadName());
                    if (pairCount == 0) {
                        ShortAnonymizedReadAlignment.correctOrientation(record);
                        readsToCorrectOrientation.put(record.getReadName(), 1);
                    }
                    else if (pairCount == 1) {
                        ShortAnonymizedReadAlignment.correctOrientation(record);
                        readsToCorrectOrientation.remove(record.getReadName());
                    }
                }
                record.setAttribute(ORIGIN_PAIR_TAG, null);
                record.setAttribute("MC", null);
                writer.addAlignment(record);
            }
        }
        // Close the files
        writer.close();
        for (final SamReader reader : readers) {
            reader.close();
        }
    }

    /**
     * Check if the read pair should be kept or not.
     * @param record The SAMRecord object representing the read alignment.
     * @return true if the read pair should be kept, false otherwise.
     */

    public boolean keepPair(SAMRecord record, Map<String, Integer> generatedReadPairs) {
        boolean hasOriginPairTag = record.getAttribute(ORIGIN_PAIR_TAG) != null;
        if(generatedReadPairs.containsKey(record.getReadName())){
            if(hasOriginPairTag){
                int originPair = record.getIntegerAttribute(ORIGIN_PAIR_TAG);
                int keptOriginPair = generatedReadPairs.get(record.getReadName());
                return keptOriginPair == originPair;
            }
            else return false;
        }
        if(hasOriginPairTag) {
            generatedReadPairs.put(record.getReadName(), record.getIntegerAttribute(ORIGIN_PAIR_TAG));
        }
        else{
            String pairName = record.getPairedReadName();
            if(updatedMatePositions.containsKey(pairName)){
                int updatedMateInfo = updatedMatePositions.get(pairName);
                record.setMateAlignmentStart(updatedMateInfo);
            }
            if(record.getMateAlignmentStart() == -1) return false;
        }
        return true;
    }
}
