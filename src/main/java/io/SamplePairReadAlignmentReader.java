package io;

import genomicelements.*;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.FastaSequenceIndexEntry;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import static utils.Operations.compare;

/**
* Read mapping reader to get paired pileups from tumor normal samples,
* only in positions covered in both alignment files
* @author Nicolas Gaitan
 */
public class SamplePairReadAlignmentReader implements Iterable<PairedPileup>, Closeable {

    public static final int DEFAULT_MINIMUM_MAPPING_QUALITY = 1;

    private SAMFileHeader normalSamHeader;
    private SAMFileHeader tumorSamHeader;
    private GenomicRegion region;
    private String platform;
    private SamReader normalSamReader;
    private SamReader tumorSamReader;
    private boolean returnNormal = false;
    private int minimumMappingQuality = DEFAULT_MINIMUM_MAPPING_QUALITY;
    private boolean includeDuplicates = true;
    private Set<String> readsToExclude;
    private int pileupMode;
    private byte[] referenceSequence = new byte[0];

    public SamplePairReadAlignmentReader(String normalFilePath, String tumorFilePath, String referenceGenomeFile, GenomicRegion region) throws IOException {
        this.region = region;
        init(new File(normalFilePath), new File(tumorFilePath), new File(referenceGenomeFile));
    }

    public SamplePairReadAlignmentReader(String normalFilePath, String tumorFilePath, String referenceGenomeFile,
                                         byte[] referenceSequence, GenomicRegion region) throws IOException {
        this.region = region;
        this.referenceSequence = referenceSequence;
        init(new File(normalFilePath), new File(tumorFilePath), new File(referenceGenomeFile));
    }

    /**
     * Initializes the iterators for the normal and tumor files. If intervals is provided, the iteration is restricted to that genomic region
     * @param normalFile
     * @param tumorFile
     * @param referenceGenomeFile
     */
    private void init(File normalFile, File tumorFile, File referenceGenomeFile) throws IOException {
        SamReaderFactory normalSamFactory = SamReaderFactory.makeDefault();
        SamReaderFactory tumorSamFactory = SamReaderFactory.makeDefault();
        normalSamFactory.referenceSequence(referenceGenomeFile);
        tumorSamFactory.referenceSequence(referenceGenomeFile);
        normalSamReader = normalSamFactory.open(normalFile);
        tumorSamReader = tumorSamFactory.open(tumorFile);
        normalSamHeader = normalSamReader.getFileHeader();
        tumorSamHeader = tumorSamReader.getFileHeader();
        // This is temporary, as it relies on all read groups having the same information
        platform = normalSamHeader.getReadGroups().get(0).getPlatform();
        readsToExclude = new HashSet<>();
    }

    public void setMinimumMappingQuality(int minimumMappingQuality) {
        this.minimumMappingQuality = minimumMappingQuality;
    }

    public void setIncludeDuplicates(boolean includeDuplicates) {
        this.includeDuplicates = includeDuplicates;
    }

    public void setReturnNormal(boolean returnNormal){
        this.returnNormal = returnNormal;
    }

    public void setReadsToExclude(Set<String> readsToExclude){
        this.readsToExclude = readsToExclude;
    }

    public void setPileupMode(int pileupMode){
        this.pileupMode = pileupMode;
    }

    /**
     * Iterator that allows sorted access to pileups coming from two different SAM files, from tumor-normal pairs,
     * as Double Records for simultaneous processing
     * @return iterator
     */
    @Override
    public Iterator<PairedPileup> iterator() {
        if (normalSamReader == null) {
            throw new IllegalStateException("Normal file reader is null");
        }
        if (tumorSamReader == null) {
            throw new IllegalStateException("Tumoral file reader is null");
        }
        return new PairedPileupIterator();
    }

    @Override
    public void forEach(Consumer<? super PairedPileup> action) {
        Iterable.super.forEach(action);
    }

    @Override
    public void close() throws IOException {
        normalSamReader.close();
        tumorSamReader.close();
    }

    public SAMFileHeader getNormalSamHeader() {
        return normalSamHeader;
    }

    public SAMFileHeader getTumorSamHeader() {
        return tumorSamHeader;
    }

    private class PairedPileupIterator implements Iterator<PairedPileup> {

        private PairedPileup nextPileup;
        private LocusPileUp nextNormalLocus;
        private LocusPileUp nextTumorLocus;
        Iterator<LocusPileUp> normalPileupIter;
        Iterator<LocusPileUp> tumorPileupIter;

        public PairedPileupIterator(){
            normalPileupIter = getLocusPileupIterator(normalSamReader);
            nextNormalLocus = normalPileupIter.hasNext() ? normalPileupIter.next() : null;
            tumorPileupIter = getLocusPileupIterator(tumorSamReader);
            nextTumorLocus = tumorPileupIter.hasNext() ? tumorPileupIter.next() : null;
            nextPileup = getNext();
        }

        private Iterator<LocusPileUp> getLocusPileupIterator(SamReader reader){
            LocusPileupIterator pileupClass = new LocusPileupIterator(reader, region.getSequenceName(), region.getStart(), region.getEnd(),
                    LocusPileupIterator.CLAIM_NEW_AND_DIFFERING_READS_PILEUP_MODE, referenceSequence);
            pileupClass.setReadsToExclude(readsToExclude);
            pileupClass.setIncludeDuplicates(includeDuplicates);
            pileupClass.setMinimumMappingQuality(minimumMappingQuality);
            return pileupClass.iterator();
        }

        @Override
        public PairedPileup next() {
            if(nextPileup == null) throw new NoSuchElementException();
            PairedPileup currentPileup = nextPileup;
            nextPileup = getNext();
            return currentPileup;
        }

        private PairedPileup getNext(){
            if(returnNormal) return getPairsAndNormals();
            return getPairs();
        }

        /**
         *
         * @return currentPileup
         */
        private PairedPileup getPairsAndNormals() {
            PairedPileup currentPileup;
            while(true){
                if (nextNormalLocus != null && nextTumorLocus != null){
                    int cmp = compare(nextNormalLocus.getSequenceIdx(), nextNormalLocus.getLocation(), nextNormalLocus.getEnd(),
                            nextTumorLocus.getSequenceIdx(), nextTumorLocus.getLocation(), nextTumorLocus.getEnd());
                    if (cmp < -1){
                        currentPileup = new PairedPileup(nextNormalLocus);
                        nextNormalLocus =  nextOrNull(normalPileupIter);
                        break;
                    }
                    else if(cmp > 1){
                        nextTumorLocus = nextOrNull(tumorPileupIter);
                    }
                    else{
                        currentPileup = new PairedPileup(nextNormalLocus, nextTumorLocus);
                        nextNormalLocus = nextOrNull(normalPileupIter);
                        nextTumorLocus = nextOrNull(tumorPileupIter);
                        break;
                    }
                }
                else if(nextNormalLocus == null && nextTumorLocus == null){
                    currentPileup = null;
                    break;
                }
                else{
                    if(nextTumorLocus == null){
                        currentPileup = new PairedPileup(nextNormalLocus);
                        nextNormalLocus = nextOrNull(normalPileupIter);
                    }
                    else{
                        currentPileup = null;
                    }
                    break;
                }
            }
            return currentPileup;
        }

        private LocusPileUp nextOrNull(Iterator<LocusPileUp> iterator){
            return iterator.hasNext() ? iterator.next() : null;
        }

        /**
         *
         * @return currentPileup
         */
        private PairedPileup getPairs() {
            PairedPileup currentPileup;
            while(true){
                if(nextNormalLocus == null || nextTumorLocus == null){
                    currentPileup = null;
                    break;
                }
                else{
                    int cmp = compare(nextNormalLocus.getSequenceIdx(), nextNormalLocus.getLocation(), nextNormalLocus.getEnd(),
                            nextTumorLocus.getSequenceIdx(), nextTumorLocus.getLocation(), nextTumorLocus.getEnd());
                    if (cmp < -1){
                        nextNormalLocus = nextOrNull(normalPileupIter);
                    }
                    else if(cmp > 1){
                        nextTumorLocus = nextOrNull(tumorPileupIter);
                    }
                    else{
                        currentPileup = new PairedPileup(nextNormalLocus, nextTumorLocus);
                        nextNormalLocus = nextOrNull(normalPileupIter);
                        nextTumorLocus = nextOrNull(tumorPileupIter);
                        break;
                    }
                }
            }
            return currentPileup;
        }

        @Override
        public boolean hasNext(){
            return nextPileup != null;
        }
    }

    //TESTS

    //Temp test for new pileup iterator
    public static void main(String[] args) throws Exception{
        String normalFilePath = args[0];
        String tumorFilePath =  args[1];
        String referenceGenome =  args[2];
//        GenomicRegion region1 = new GenomicRegionBaseImpl("1", 196256968, 196329019);
//        GenomicRegion region1 = new GenomicRegionBaseImpl("1", 1, 249250621);
        int nTests = 10;
        boolean passedAll = true;
        boolean passed;
        for(int i = 0; i < nTests; i++){
            passed = verifyPileupMethod(tumorFilePath, referenceGenome);
            if(!passed){
                System.out.println("Did not pass test #" + i);
                passedAll = false;
            }
        }
        if(passedAll) System.out.println("Passed all tests");
//        try(SamplePairReadAlignmentReader pairPileupReader = new SamplePairReadAlignmentReader(normalFilePath, tumorFilePath, referenceGenome, region1);
//            IndexedFastaSequenceFile referenceWalker = new IndexedFastaSequenceFile(new File(referenceGenome))){
//            pairPileupReader.setReturnNormal(true);
//            Iterator<PairedPileup> pileupIterator = pairPileupReader.iterator();
//            while (pileupIterator.hasNext()) {
//                PairedPileup pileup = pileupIterator.next();
//                int pos = pileup.getReferencePos();
//                System.out.println("POS=" + pos);
//
//            }
//        }
    }
    public static boolean verifyPileupMethod(String tumorFilePath, String refGenome) throws IOException {
        boolean passed = true;
        int testWindow = 1_000_000;
        Map<Integer, Set<String>> truthPileupReadNames = new HashMap<>();
        Map<Integer, Set<String>> myPileupReadNames = new HashMap<>();
        long myPileupExecTime = 0;
        long samLocusPileupExecTime = 0;
        //List<String> truthReads = new ArrayList<>();
        Random rand = new Random();
        IndexedFastaSequenceFile reference = new IndexedFastaSequenceFile(new File(refGenome));
        FastaSequenceIndex refIndexes = reference.getIndex();
        Map<String, Integer> refSequencesLength = new HashMap<>();
        for (FastaSequenceIndexEntry refEntry : refIndexes){
            refSequencesLength.put(refEntry.getContig(), (int) refEntry.getSize());
        }
        reference.close();
        String chr = Integer.toString(rand.nextInt(1, 22));
        int start = rand.nextInt(1, refSequencesLength.get(chr) - testWindow);
        int end = start + testWindow;
//        String chr = "1";
//        int start = 180513307;
//        int end = 181513307;
        System.out.println("Testing window=" + chr + " " + start + " " + end);
        System.out.println("Testing SamLocusIterator");
        int samLocusPositions = 0;
        try(SamReader tumorSamReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(refGenome))
                .open(new File(tumorFilePath))){
            IntervalList intervalList = new IntervalList(tumorSamReader.getFileHeader());
            Interval intervalRegion = new Interval(chr, start, end);
            intervalList.add(intervalRegion);
            SamLocusIterator samIter = new SamLocusIterator(tumorSamReader, intervalList, true);
            //SAMRecordIterator samIter = tumorSamReader.query(region1.getSequenceName(), region1.getStart(), region1.getEnd(), false);
            samIter.setSamFilters(null);
            samIter.setEmitUncoveredLoci(false);
            samIter.setMappingQualityScoreCutoff(1);
            samIter.setIncludeNonPfReads(true);
            long samLocusPileupBegin = System.currentTimeMillis();
            while (samIter.hasNext()) {
                SamLocusIterator.LocusInfo rec = samIter.next();
//                if(rec.getStart()<=196256969 && rec.getEnd()>=196256969){
//                    truthReads.add(rec.getReadName());
//                }

//                Set<String> readNames = rec.getRecordAndOffsets().stream()
//                        //.filter(v -> v.getRecord().getMappingQuality() > 0)
//                        .map(SamLocusIterator.RecordAndOffset::getReadName)
//                        //.map(SAMRecord::getReadName)
//                        //.distinct()
//                        .collect(Collectors.toSet());
//                truthPileupReadNames.put(rec.getPosition(), readNames);

                samLocusPositions++;
            }
            long samLocusPileupEnd = System.currentTimeMillis();
            samLocusPileupExecTime = samLocusPileupEnd - samLocusPileupBegin;
        }
//        List<String> testReads = new ArrayList<>();
        System.out.println("Testing MyPileup");
        int myPileupPositions = 0;
        try(SamReader tumorSamReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(refGenome))
                .open(new File(tumorFilePath))){
//            SAMRecordIterator samIter = tumorSamReader.iterator();
            SAMRecordIterator samIter = tumorSamReader.query(chr, start, end, false);
            LocusPileupIterator pileupIterator = new LocusPileupIterator(samIter, start, end);
            pileupIterator.setMinimumMappingQuality(1);
            pileupIterator.setIncludeDuplicates(true);
            Iterator<LocusPileUp> it = pileupIterator.iterator();
            long myPileupBegin = System.currentTimeMillis();
            while (it.hasNext()) {
                LocusPileUp pileup = it.next();
//                int pos = pileup.getLocation();
//                for(SAMRecord rec : pileup.getRecords()){
//                   if(rec.getReadName().equals("HWI-ST805:328:C1LRFACXX:3:2314:14375:100552")) {
//                       System.out.println(rec.toString());
//                       System.exit(0);
//                   }
//                }
//                Set<String> readNames = pileup.getRecords().stream()

//                Set<String> readNames = pileup.getRecords().stream()
//                        //.filter(v -> !v.isSecondaryOrSupplementary())
//                        .map(SAMRecord::getReadName)
//                        //.distinct()
//                        .collect(Collectors.toSet());
//                myPileupReadNames.put(pos, readNames);

//                System.out.println("Reading pileup at pos=" + pos + " with " + pileup.size() + " reads");
//                if(pos==196256969){
//                    readNames.forEach(System.out::println);
//                    testReads.addAll(readNames);
//                    System.exit(0);
//                }
//                System.out.println("POS=" + pos);
                myPileupPositions++;
            }
            long myPileupEnd = System.currentTimeMillis();
            myPileupExecTime = myPileupEnd - myPileupBegin;
        }
        System.out.print("SamLocusPileup executed in : " + samLocusPileupExecTime + "ms ");
        System.out.println("across " + samLocusPositions + " positions");
        System.out.print("MyPileup executed in : " + myPileupExecTime + "ms ");
        System.out.println("across " + myPileupPositions + " positions");
        Set<String> allTruthReads = new HashSet<>();
        for(int i = start; i <= end; i++){
            Set<String> posTruthReads = truthPileupReadNames.getOrDefault(i, new HashSet<>());
            allTruthReads.addAll(posTruthReads);
            Set<String> posMyPileupReads = myPileupReadNames.getOrDefault(i, new HashSet<>());
            Set<String> commonReads = new HashSet<>(posTruthReads);
            commonReads.retainAll(posMyPileupReads);
            boolean areEqualAt = (commonReads.size() == posTruthReads.size()) && (commonReads.size() == posMyPileupReads.size());
            if(!areEqualAt) {
                passed = false;
                Set<String> posMyPileupReadsDiff = new HashSet<>(posMyPileupReads);
                posMyPileupReadsDiff.removeAll(commonReads);
                Set<String> posTruthReadsDiff = new HashSet<>(posTruthReads);
                posTruthReadsDiff.removeAll(commonReads);
                System.out.print("Reads are not equal at pos=" + i);
                System.out.println(" where common reads are: " + commonReads.size() + ", SamLocusIterator reads are: " + posTruthReads.size() +
                        ", myPileup reads are: " + posMyPileupReads.size());
                if(!posTruthReadsDiff.isEmpty()) {
                    System.out.println("Different reads in SamLocusIterator are: ");
                    posTruthReadsDiff.forEach(System.out::println);
                }
                if(!posMyPileupReadsDiff.isEmpty()) {
                    System.out.println("Different reads in myPileup are: ");
                    posMyPileupReadsDiff.forEach(System.out::println);
                }
            }
        }
        myPileupReadNames.clear();
        System.out.println("Testing new on-demand read retrieval functionality");
        Set<String> onDemandTestReads = new HashSet<>();
        int onDemandPileupPositions = 0;
        long onDemandPileupExecTime;
        try(SamReader tumorSamReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(refGenome))
                .open(new File(tumorFilePath))){
            SAMRecordIterator samIter = tumorSamReader.query(chr, start, end, false);
            LocusPileupIterator pileupIterator = new LocusPileupIterator(samIter, start, end);
            pileupIterator.setMinimumMappingQuality(1);
            pileupIterator.setIncludeDuplicates(true);
            Iterator<LocusPileUp> it = pileupIterator.iterator();
            long onDemandPileupBegin = System.currentTimeMillis();
            while (it.hasNext()) {
                LocusPileUp pileup = it.next();
//                Set<String> readNames = pileup.claimNewReadsOnPileup().stream()
//                        .map(SAMRecord::getReadName)
//                        .collect(Collectors.toSet());
//                onDemandTestReads.addAll(readNames);
                onDemandPileupPositions++;
            }
            long onDemandPileupEnd = System.currentTimeMillis();
            onDemandPileupExecTime = onDemandPileupEnd - onDemandPileupBegin;
        }
        System.out.println("onDemandPileup executed across: " + onDemandPileupPositions + " positions" );
        System.out.println("onDemandPileup executed in : " + onDemandPileupExecTime + "ms while " +
                "SamLocusPileup executed in : " + samLocusPileupExecTime + "ms ,and " +
                "MyPileup executed in : " + myPileupExecTime + "ms ");
        if (!onDemandTestReads.equals(allTruthReads)) {
            Set<String> uniqueToTestReads = new HashSet<>(onDemandTestReads);
            uniqueToTestReads.removeAll(allTruthReads);
            Set<String> uniqueToAllTruthReads = new HashSet<>(allTruthReads);
            uniqueToAllTruthReads.removeAll(onDemandTestReads);
            if(!uniqueToTestReads.isEmpty()) {
                System.out.println("Unique elements in onDemandPileupReads: " + uniqueToTestReads.size() + " from " + onDemandTestReads.size() + " total" );
                System.out.println("Showing sample of unique elements in onDemandPileup: ");
                for(int i = 0; i < 10; i++){
                    String readName = uniqueToTestReads.stream().skip(new Random().nextInt(uniqueToTestReads.size())).findFirst().orElse(null);
                    if(readName == null) break;
                    System.out.println(readName);
                }
            }
            if(!uniqueToAllTruthReads.isEmpty()) {
                System.out.println("Unique elements in allTruthReads: " + uniqueToAllTruthReads.size() + " from " + allTruthReads.size() + " total");
                System.out.println("Showing sample of unique elements in allTruthReads: ");
                for(int i = 0; i < 10; i++){
                    String readName = uniqueToAllTruthReads.stream().skip(new Random().nextInt(uniqueToAllTruthReads.size())).findFirst().orElse(null);
                    if(readName == null) break;
                    System.out.println(readName);
                }
            }
            passed = false;
        } else {
            System.out.println("onDemandPileupReads and allTruthReads have the same elements.");
        }
        return passed;
    }
}
