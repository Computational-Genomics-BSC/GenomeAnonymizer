package io;

import genomicelements.*;
import htsjdk.samtools.*;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Consumer;

import static analysis.GenomeAnonymizer.DEFAULT_MIN_MAPPING_QUALITY;
import static utils.Operations.compare;

/**
* Read mapping reader to get paired pileups from tumor normal samples,
* only in positions covered in both alignment files
* @author Nicolas Gaitan
 */
public class SamplePairReadAlignmentReader implements Iterable<PairedPileup>, Closeable {

    private SAMFileHeader normalSamHeader;
    private SAMFileHeader tumorSamHeader;
    private GenomicRegion region;
    private String platform;
    private SamReader normalSamReader;
    private SamReader tumorSamReader;

    private int minimumMappingQuality = DEFAULT_MIN_MAPPING_QUALITY;
    private boolean includeDuplicates = false;

    private Set<String> readsToExclude;
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

    public SAMFileHeader getNormalSamHeader() {
        return normalSamHeader;
    }

    public SAMFileHeader getTumorSamHeader() {
        return tumorSamHeader;
    }

    public void setIncludeDuplicates(boolean includeDuplicates) {
        this.includeDuplicates = includeDuplicates;
    }

    public void setReadsToExclude(Set<String> readsToExclude){
        this.readsToExclude = readsToExclude;
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
            //Creates a LocusPileupIterator for the given SAM file reader in CLAIM_NEW_AND_DIFFERING_READS_ONLY_PILEUP_MODE
            LocusPileupIterator pileupClass = new LocusPileupIterator(reader, region.getSequenceName(), region.getStart(), region.getEnd(), referenceSequence);
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

        /**
         *
         * @return currentPileup
         */
        private PairedPileup getNext() {
            PairedPileup currentPileup;
            if (nextNormalLocus != null && nextTumorLocus != null){
                int cmp = compare(nextNormalLocus.getSequenceIdx(), nextNormalLocus.getLocation(), nextNormalLocus.getEnd(),
                        nextTumorLocus.getSequenceIdx(), nextTumorLocus.getLocation(), nextTumorLocus.getEnd());
                if (cmp < -1){
                    currentPileup = new PairedPileup(nextNormalLocus, true);
                    nextNormalLocus = nextOrNull(normalPileupIter);
                }
                else if(cmp > 1){
                    currentPileup = new PairedPileup(nextTumorLocus, false);
                    nextTumorLocus = nextOrNull(tumorPileupIter);
                }
                else{
                    currentPileup = new PairedPileup(nextNormalLocus, nextTumorLocus);
                    nextNormalLocus = nextOrNull(normalPileupIter);
                    nextTumorLocus = nextOrNull(tumorPileupIter);
                }
            }
            else if(nextNormalLocus == null && nextTumorLocus == null){
                return null;
            }
            else{
                if(nextTumorLocus == null){
                    currentPileup = new PairedPileup(nextNormalLocus, true);
                    nextNormalLocus = nextOrNull(normalPileupIter);
                }
                else{
                    currentPileup = new PairedPileup(nextTumorLocus, false);
                    nextTumorLocus = nextOrNull(tumorPileupIter);
                }
            }
            return currentPileup;
        }

        private LocusPileUp nextOrNull(Iterator<LocusPileUp> iterator){
            return iterator.hasNext() ? iterator.next() : null;
        }

        @Override
        public boolean hasNext(){
            return nextPileup != null;
        }
    }
}
