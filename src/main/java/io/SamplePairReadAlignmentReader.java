package io;

import GenomicElements.PairedPileup;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.function.Consumer;

/**
* Simple read mapping reader to get paired pileups from tumor normal samples,
* only from positions covered in both alignment files
* @author Nicolas Gaitan
 */
public class SamplePairReadAlignmentReader implements Iterable<PairedPileup>, Closeable {

    public static final int MINIMUM_MAPPING_QUALITY = 0;
    public static final boolean INCLUDE_NON_PF_READS = true;

    private SAMFileHeader normalSamHeader;
    private SAMFileHeader tumorSamHeader;
    private SamLocusIterator iterTumor;
    private SamLocusIterator iterNormal;
    // private Iterator<PairedPileup> ;

    public SamplePairReadAlignmentReader(String normalFilePath, String tumorFilePath, String referenceGenome) throws IOException {
        init(new File(normalFilePath), new File(tumorFilePath), new File(referenceGenome), null);
    }

    public SamplePairReadAlignmentReader(String normalFilePath, String tumorFilePath, String referenceGenome, IntervalList intervals) throws IOException {
        init(new File(normalFilePath), new File(tumorFilePath), new File(referenceGenome), intervals);
    }

    private void init(File normalFile, File tumorFile, File referenceGenome, IntervalList intervals) {
        SamReaderFactory normalSamFactory = SamReaderFactory.makeDefault();
        SamReaderFactory tumorSamFactory = SamReaderFactory.makeDefault();
        normalSamFactory.referenceSequence(referenceGenome);
        tumorSamFactory.referenceSequence(referenceGenome);
        SamReader normalSamReader = tumorSamFactory.open(normalFile);
        SamReader tumorSamReader = tumorSamFactory.open(tumorFile);
            normalSamHeader = normalSamReader.getFileHeader();
            tumorSamHeader = tumorSamReader.getFileHeader();
            if (intervals != null){
                iterTumor = new SamLocusIterator(tumorSamReader, intervals);
                iterNormal = new SamLocusIterator(normalSamReader, intervals);
            }
            else{
                iterTumor = new SamLocusIterator(tumorSamReader);
                iterNormal = new SamLocusIterator(normalSamReader);
            }
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
            throw new IllegalStateException("File reader for normal file is null");
        }
        if (iterTumor == null) {
            throw new IllegalStateException("File reader for tumor file is null");
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
            if (!iterNormal.hasNext()){
                throw new IllegalStateException("Normal set alignment file is empty");
            }
            if (!iterTumor.hasNext()){
                throw new IllegalStateException("Tumor set alignment file is empty");
            }
            nextNormalLocus = iterNormal.next();
            nextTumorLocus = iterTumor.next();
            nextPileup = getNextOrAdvance();
            // PairedPileup nextPileup = new PairedPileup(iterNormal, iterTumor);
        }

        @Override
        public PairedPileup next() {
            if(nextPileup == null) throw new NoSuchElementException();
            PairedPileup currentPileup = nextPileup;
            nextPileup = getNextOrAdvance();
            return currentPileup;
        }

        /**
         * <pre> Pileups must pertain to positions in the same sequence </pre>
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
                    if (nextNormalLocus.getPosition() < nextTumorLocus.getPosition()){
                        nextNormalLocus = iterNormal.next();
                    }
                    else if(nextNormalLocus.getPosition() > nextTumorLocus.getPosition()){
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
