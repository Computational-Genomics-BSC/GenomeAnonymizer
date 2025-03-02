package analysis;

import genomicelements.*;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.FastaSequenceIndexEntry;
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
import static analysis.GenomeAnonymizer.DEFAULT_MIN_MAPPING_QUALITY;

/**
 * AnonymizerAlgorithm implementation for short read data
 * @author Nicolas Gaitan
 * @author Rodrigo Martin
 */
public class ShortReadAnonymizer implements AnonymizerAlgorithm {

    private static final Logger LOGGER = Logger.getLogger(ShortReadAnonymizer.class.getName());

    public static final int NORMAL_DATASET_IDX = 0;
    public static final int TUMORAL_DATASET_IDX = 1;

    private String normalPath;
    private String tumorPath;
    private String refGenomePath;
    private String outputPrefix;

    private List<GenomicRegion> genomicPartitions = new ArrayList<>();
    private Map<String, byte[]> referenceSequences = new HashMap<>();

    // Set that contains all the reads that will be excluded from the result (e.g. Unmapped and MAPQ < filter)
    private Set<String> readsToExclude;
    private File canvasNormal;
    private File canvasTumoral;
    private SamReaderFactory factory;
    private SAMFileWriter normalWriter;
    private SAMFileWriter tumoralWriter;

    private int minimumMappingQuality = DEFAULT_MIN_MAPPING_QUALITY;
    private int hashSalt;

    private int threads = 1;
    //Query regions are only used to merge with a synthetic genome
    private List<GenomicRegion> queryRegions;

    //TODO: Save partitioned files names sorted by region coordinates

    private Map<String, Map<Integer, PairCalledVariation>> somaticVariantsToKeep = new HashMap<>();

    public ShortReadAnonymizer(String normalPath, String tumorPath, String refGenomePath, String outputPrefix) {
        this.normalPath = normalPath;
        this.tumorPath = tumorPath;
        this.refGenomePath = refGenomePath;
        this.outputPrefix = outputPrefix;
        readsToExclude = new HashSet<>();
        factory = SamReaderFactory.makeDefault();
        factory.setUseAsyncIo(true);
        factory.validationStringency(ValidationStringency.SILENT);
        queryRegions = new ArrayList<>();
        // Initialize hash_salt with a random value
        hashSalt = GlobalRandom.getInstance().nextInt();
    }

    @Override
    public void setMinimumMappingQuality(int minimumMappingQuality) {
        this.minimumMappingQuality = minimumMappingQuality;
    }

    public void setThreadNumber(int threads) {
        this.threads = threads;
    }

    public void setVCFVariantsToKeep(Map<String, Map<Integer, PairCalledVariation>> variantsToKeep){
        this.somaticVariantsToKeep = variantsToKeep;
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
        // Sort queryRegions with the updated information on sequence order from the reference genome
        queryRegions.forEach(v -> v.setSequenceIdx(refSequenceOrder.get(v.getSequenceName())));
        queryRegions.sort(Comparator.comparingInt(GenomicRegion::getSequenceIdx)
                .thenComparing(GenomicRegion::getStart));
        // Compute genome partitions
        sequences.sort(Comparator.comparing(FastaSequenceIndexEntry::getSize));
        Collections.reverse(sequences);
        int[] partitionsPerSequence = new int[sequences.size()];
        Arrays.fill(partitionsPerSequence, 1);
//        long basesPerThread = genomeSize / (nThreads * 10L);
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
                        (samRecord.getMappingQuality() < minimumMappingQuality && !samRecord.isSecondaryOrSupplementary())){
                    readsToExcludeInPartition.add(samRecord.getReadName());
                }
            }
        }
        return readsToExcludeInPartition;
    }

    @Override
    public void anonymizeReads() {
        ExecutorService executorService = Executors.newFixedThreadPool(threads);
        List<CompletableFuture<Void>> futures = new ArrayList<>();
        for(GenomicRegion genomicPartition : genomicPartitions){
            CompletableFuture<Void> future = CompletableFuture.runAsync (() -> {
                new PartitionProviderRunner(genomicPartition).run();
            }, executorService);
            futures.add(future);
        }
        CompletableFuture.allOf(futures.toArray(new CompletableFuture[0])).join();
        executorService.shutdown();
    }

    class PartitionProviderRunner implements Runnable {
        private final GenomicRegion genomicPartition;

        public PartitionProviderRunner(GenomicRegion genomicPartition) {
            this.genomicPartition = genomicPartition;
        }

        @Override
        public void run() {
            try (SAMFileWriter partitionNormalWriter = openSingleOutputStream(normalPath, genomicPartition, true);
                 SAMFileWriter partitionTumoralWriter = openSingleOutputStream(tumorPath, genomicPartition, false);
                 AnonymizedReadAlignmentProvider anonymizedReadProvider = new AnonymizedReadAlignmentProvider();) {
                anonymizedReadProvider.setReadsToExclude(readsToExclude);
                anonymizedReadProvider.setRefSequence(referenceSequences.get(genomicPartition.getSequenceName()));
                anonymizedReadProvider.setMinMappingQuality(minimumMappingQuality);
                anonymizedReadProvider.setVCFVariantsToKeep(somaticVariantsToKeep);
                anonymizedReadProvider.init(normalPath, tumorPath, refGenomePath, genomicPartition);
                long startcallVariation = System.currentTimeMillis();
                for (AnonymizedRead anonymizedRead : anonymizedReadProvider) {
                    boolean isNormalDataset = anonymizedRead.isFromNormalDataset();
                    long startWriteRead = System.currentTimeMillis();
                    writeRead(isNormalDataset ? partitionNormalWriter : partitionTumoralWriter,
                            anonymizedRead.getAnonymizedSamRecord());
                    long endWriteRead = System.currentTimeMillis();
                    long writeReadTime = endWriteRead - startWriteRead;
                    anonymizedReadProvider.METHOD_TIME_MAP.compute("writeReadTime",  (k,v) -> v == null ?
                            writeReadTime :
                            v + writeReadTime);
                }
                long endcallVariation = System.currentTimeMillis();
                //DEBUG
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
                //DEBUG
            }
            catch (Exception e) {
                LOGGER.log(Level.SEVERE,
                        "Exception in thread Anonymizing reads to exclude in region: "
                                + genomicPartition.getSequenceName()
                                + " " + genomicPartition.getStart() + " " + genomicPartition.getEnd(),
                        e);
                System.exit(1);
//                    throw new RuntimeException(e);
            }
        }
    }

    private SAMFileWriter openSingleOutputStream(String path, GenomicRegion region, boolean isNormalDataset) throws IOException {
        int idx = isNormalDataset ? NORMAL_DATASET_IDX : TUMORAL_DATASET_IDX;
        String suffix = "_" + region.toString() ;
        File outputFile = new File(getBAMOutputName(outputPrefix+suffix, idx));
        SAMFileWriterFactory factory = new SAMFileWriterFactory();
        //factory.setCreateIndex(true);
        factory.setCompressionLevel(1);
        factory.setUseAsyncIo(true);
        SAMFileHeader fileHeader = buildFileHeader(path, isNormalDataset ? "_N" : "_T");
        SAMFileWriter writer = factory.makeBAMWriter(fileHeader, true, outputFile);
        writer.setSortOrderChecking(false);
        return writer;
    }

    public String getBAMOutputName(String outputPrefix, int datasetIdx){
        String anonTag = ".anonymized";
        String datasetIdStr = datasetIdx == NORMAL_DATASET_IDX ? ".N" : ".T";
        String extension = BAM_FILE;
        return outputPrefix + anonTag + datasetIdStr + extension;
    }

    private void writeRead(SAMFileWriter samWriter, SAMRecord read) {
        // Set the read name to the hash of the read name + salt
        read.setReadName(UUID.nameUUIDFromBytes((read.getReadName() + hashSalt).getBytes()).toString());
        // Set the read group to the same value as the read group of the header
        read.setAttribute("RG", samWriter.getFileHeader().getReadGroups().get(0).getId());
        samWriter.addAlignment(read);
    }

    private void logWrittenReads(int writtenReads, String datasetName) {
        if (writtenReads % 10_000_000 == 0) LOGGER.info("Written " + writtenReads + " to " + datasetName);
    }

    //TODO: @Rodrigo Martin Change this to be the merger code
    public void mergeAnonymizedReads() throws IOException {
        // TODO: Make sure the reference genome is the same for this.canvasTumoral and this.tumoralWriter and for this.canvasNormal and this.normalWriter
        try {
            openFinalOutputStreams();
            ExecutorService exec = Executors.newFixedThreadPool(2);
            Future<?> tumoralFuture = exec.submit(() -> {
                try {
                    mergeAnonymizedReads(new File(tumorPath), this.canvasTumoral, this.tumoralWriter, new File(refGenomePath),
                            false);
                } catch (Exception e) {
                    LOGGER.log(Level.SEVERE,
                            "Exception in thread merging anonymized reads in tumoral dataset", e);
                    throw new RuntimeException(e);
                }
            });
            Future<?> normalFuture = exec.submit(() -> {
                try {
                    mergeAnonymizedReads(new File(normalPath), this.canvasNormal, this.normalWriter, new File(refGenomePath),
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

    private void openFinalOutputStreams() throws IOException {
        File normalOutputFile = new File(getBAMOutputName(outputPrefix, NORMAL_DATASET_IDX));
        File tumoralOutputFile = new File(getBAMOutputName(outputPrefix, TUMORAL_DATASET_IDX));
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
                    readAlnId = ShortAnonymizedReadAlignment.generateReadAlnId(samRecord);
                    if (readsToExclude.contains(readName) || returnedReads.contains(readAlnId)) {
                        samRecord = null;
                    }
                }
                SAMRecord newSamRecord = samRecord;
//                if (readGermlinesToAnonymize.containsKey(readAlnId)) {
//                    ShortAnonymizedReadAlignment anonymizedReadAlignment = new ShortAnonymizedReadAlignment(samRecord);
//                    anonymizedReadAlignment.setReferenceContigSequence(referenceContigSequence);
//                    anonymizedReadAlignment.setSignalsToAnonymize(readGermlinesToAnonymize.get(readAlnId));
////                    anonymizedReadAlignment.anonymizeVariants();
//                    newSamRecord = anonymizedReadAlignment.getAnonymizedSamRecord();
//                }
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
