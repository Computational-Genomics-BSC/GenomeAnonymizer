package io;

import genomicelements.PairedPileup;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import htsjdk.tribble.SimpleFeature;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.function.Consumer;

import static utils.Operations.compare;

/**
* Read mapping reader to get paired pileups from tumor normal samples,
* only in positions covered in both alignment files
* @author Nicolas Gaitan
 */
public class SamplePairReadAlignmentReader implements Iterable<PairedPileup>, Closeable {

    public static final int MINIMUM_MAPPING_QUALITY = 0;
    public static final boolean INCLUDE_NON_PF_READS = true;

    private SAMFileHeader normalSamHeader;
    private SAMFileHeader tumorSamHeader;
    private String platform;
    private SamLocusIterator iterTumor;
    private SamLocusIterator iterNormal;
    // private Iterator<PairedPileup> ;

    public SamplePairReadAlignmentReader(String normalFilePath, String tumorFilePath, String referenceGenome) throws IOException {
        init(new File(normalFilePath), new File(tumorFilePath), new File(referenceGenome), null);
    }

    public SamplePairReadAlignmentReader(String normalFilePath, String tumorFilePath, String referenceGenome, SimpleFeature region) throws IOException {
        init(new File(normalFilePath), new File(tumorFilePath), new File(referenceGenome), region);
    }

    /**
     * Initializes the iterators for the normal and tumor files. If intervals is provided, the iteration is restricted to that genomic region
     * @param normalFile
     * @param tumorFile
     * @param referenceGenome
     * @param region Can be any range of valid genomic regions, but is intended to be a list of one region for multithreading
     */
    private void init(File normalFile, File tumorFile, File referenceGenome, SimpleFeature region)throws IOException {
        SamReaderFactory normalSamFactory = SamReaderFactory.makeDefault();
        SamReaderFactory tumorSamFactory = SamReaderFactory.makeDefault();
        normalSamFactory.referenceSequence(referenceGenome);
        tumorSamFactory.referenceSequence(referenceGenome);
        SamReader normalSamReader = tumorSamFactory.open(normalFile);
        SamReader tumorSamReader = tumorSamFactory.open(tumorFile);
        normalSamHeader = normalSamReader.getFileHeader();
        tumorSamHeader = tumorSamReader.getFileHeader();
        // This is temporary, as it relies on all read groups having the same information
        platform = normalSamHeader.getReadGroups().get(0).getPlatform();
        if (region != null){
            assert (getNormalSamHeader().getSequenceDictionary().isSameDictionary(getTumorSamHeader().getSequenceDictionary())):
                    "Headers have different sequence dictionaries";
            IntervalList intervalList = new IntervalList(this.getNormalSamHeader());
            Interval intervalRegion = new Interval(region.getContig(), region.getStart(), region.getEnd());
            intervalList.add(intervalRegion);
            //DEBUG
            //intervalList.forEach(v -> System.out.println(v.getContig() + " " + v.getStart() + " " + v.getEnd()));
            //System.exit(0);
            //DEBUG
            iterTumor = new SamLocusIterator(tumorSamReader, intervalList, true);
            iterNormal = new SamLocusIterator(normalSamReader, intervalList, true);
        }
        else{
            iterTumor = new SamLocusIterator(tumorSamReader);
            iterNormal = new SamLocusIterator(normalSamReader);
        }
        // DEBUG
        //System.out.println("hola");
//        if (iterTumor==null){
//            System.out.println("iterTumor is null");
//        }
        // DEBUG
        // DEBUG
//        if (iterNormal==null){
//            System.out.println("iterNormal is null");
//        }
        // DEBUG
        setFilteringForSAMIterators(iterNormal);
        setFilteringForSAMIterators(iterTumor);
    }

    private void setFilteringForSAMIterators(SamLocusIterator iterator){
        iterator.setSamFilters(null);
        iterator.setEmitUncoveredLoci(false);
        iterator.setMappingQualityScoreCutoff(MINIMUM_MAPPING_QUALITY);
        iterator.setIncludeNonPfReads(INCLUDE_NON_PF_READS);
    }

    /**
     * Iterator that allows sorted access to pileups coming from two different SAM files, from tumor-normal pairs,
     * as Double Records for simultaneous processing
     * @return iterator
     */
    @Override
    public Iterator<PairedPileup> iterator() {
        if (iterNormal == null) {
            throw new IllegalStateException("Normal file reader is closed");
        }
        if (iterTumor == null) {
            throw new IllegalStateException("Tumoral file reader is closed");
        }
        return new PairedPileupIterator();
    }

    @Override
    public void forEach(Consumer<? super PairedPileup> action) {
        Iterable.super.forEach(action);
    }

    @Override
    public void close() throws IOException {
        iterNormal.close();
        iterTumor.close();
    }

    public SAMFileHeader getNormalSamHeader() {
        return normalSamHeader;
    }

    public SAMFileHeader getTumorSamHeader() {
        return tumorSamHeader;
    }

    private class PairedPileupIterator implements Iterator<PairedPileup> {

        private PairedPileup nextPileup;
        private LocusInfo nextNormalLocus;
        private LocusInfo nextTumorLocus;

        public PairedPileupIterator(){
            nextNormalLocus = iterNormal.next();
            nextTumorLocus = iterTumor.next();
            nextPileup = getNextOrAdvance();
        }

        @Override
        public PairedPileup next() {
            if(nextPileup == null) throw new NoSuchElementException();
            PairedPileup currentPileup = nextPileup;
            nextPileup = getNextOrAdvance();
            return currentPileup;
        }

//        /**
//         * <pre> Pileups must pertain to positions in the same sequence </pre>
//         * @return currentPileup
//         */
//        private PairedPileup getNextOrAdvance() {
//            PairedPileup currentPileup;
//            while(true){
//                if(nextNormalLocus == null || nextTumorLocus == null){
//                    currentPileup = null;
//                    break;
//                }
//                else{
//                    if (nextNormalLocus.getPosition() < nextTumorLocus.getPosition()){
//                        nextNormalLocus = iterNormal.next();
//                    }
//                    else if(nextNormalLocus.getPosition() > nextTumorLocus.getPosition()){
//                        nextTumorLocus = iterTumor.next();
//                    }
//                    else{
//                        currentPileup = new PairedPileup(nextNormalLocus, nextTumorLocus);
//                        nextNormalLocus = iterNormal.next();
//                        nextTumorLocus = iterTumor.next();
//                        break;
//                    }
//                }
//            }
//            return currentPileup;
//        }

        /**
         *
         * @return currentPileup
         */
        private PairedPileup getNextOrAdvance() {
            PairedPileup currentPileup;
            while(true){
                if(nextNormalLocus == null || nextTumorLocus == null){
                    currentPileup = null;
                    break;
                }
                else{
                    int cmp = compare(nextNormalLocus.getSequenceIndex(), nextNormalLocus.getPosition(), nextNormalLocus.getEnd(),
                            nextTumorLocus.getSequenceIndex(), nextTumorLocus.getPosition(), nextTumorLocus.getEnd());
                    if (cmp < -1){
                        nextNormalLocus = iterNormal.next();
                    }
                    else if(cmp > 1){
                        nextTumorLocus = iterTumor.next();
                    }
                    else{
                        currentPileup = new PairedPileup(nextNormalLocus, nextTumorLocus);
                        nextNormalLocus = iterNormal.next();
                        nextTumorLocus = iterTumor.next();
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
}
