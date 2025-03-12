package genomicelements;

import htsjdk.samtools.*;
import htsjdk.samtools.SAMRecord;
import utils.MapCacheFIFO;

import java.util.*;
import java.util.function.Consumer;

import genomicelements.PileupRead.PileupReadStatus;

import static utils.Operations.overlap;

/**
 * @author Nicolas Gaitan
 * Iterable interface to traverse locus pileups over a stream of alignments, as SAMRecord objects
 */
public class LocusPileupIterator implements Iterable<LocusPileUp> {

    public static final Set<Character> ALPHABET = new HashSet<>(
            Arrays.asList(
                    'A', 'T', 'C', 'G', 'N'
            )
    );

    public static final int CLAIM_ALL_READS_PILEUP_MODE = 0;
    public static final int CLAIM_NEW_READS_ONLY_PILEUP_MODE = 1;
    public static final int CLAIM_NEW_AND_DIFFERING_READS_ONLY_PILEUP_MODE = 2;

    public static final int DEFAULT_PILEUP_READ_LIMIT = 100_000;
    public static final int DEFAULT_READ_MINIMUM_MAPQ = 0;

    private SamReader samReaderStream;
    private MapCacheFIFO<Integer, LocusPileUp> cache;
    private int start;
    private int end;
    private int nextReferencePosition;
    private OnPileupQueue readQueue;
    private Set<String> readsToExclude;

    //Pileup behaviour modifiers
    private int minimumMappingQuality = DEFAULT_READ_MINIMUM_MAPQ;
    private boolean includeDuplicates = false;
    private int pileupMode;

    //This reference sequence is 0-based, given it comes in a byte array
    private byte[] refSequence;
    private String sequenceName;

    //TODO: Handle pileup when traversing different chromosomes or specify it cannot be used for that

    public LocusPileupIterator(SamReader samReaderStream, String sequenceName, int start, int end){
        this(samReaderStream, sequenceName, start, end, CLAIM_ALL_READS_PILEUP_MODE);
    }

    public LocusPileupIterator(SamReader samReaderStream, String sequenceName, int start, int end, byte[] refSequence) throws IllegalStateException{
        this(samReaderStream, sequenceName, start, end, CLAIM_NEW_AND_DIFFERING_READS_ONLY_PILEUP_MODE);
        if(this.pileupMode == CLAIM_NEW_AND_DIFFERING_READS_ONLY_PILEUP_MODE && refSequence == null) {
            throw new IllegalStateException("Reference sequence must be provided when using the new and different reads pileup mode.");
        }
        setRefSequence(refSequence);
    }

    public LocusPileupIterator(SamReader samReaderStream, String sequenceName, int start, int end, int pileupMode) throws IllegalStateException{
        setPileupMode(pileupMode);
        if(this.pileupMode != CLAIM_NEW_AND_DIFFERING_READS_ONLY_PILEUP_MODE && pileupMode != CLAIM_ALL_READS_PILEUP_MODE) {
            throw new IllegalArgumentException("Invalid pileup mode: " + pileupMode + ". Valid modes are: " +
                    "0 - Claim all reads pileup mode, 1 - Claim new and different reads pileup mode.");
        }
        init(samReaderStream, sequenceName, start, end);
    }

    public void init(SamReader samReaderStream, String sequenceName, int start, int end) throws IllegalStateException{
        this.samReaderStream = samReaderStream;
        this.sequenceName = sequenceName;
        // 1-based start and end of the region
        this.start = start;
        this.end = end;
        this.nextReferencePosition = this.start;
        this.cache = new MapCacheFIFO<>();
        this.readQueue = new OnPileupQueue();
        this.readsToExclude = new HashSet<>();
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public int getMinimumMappingQuality() {
        return minimumMappingQuality;
    }

    public boolean includesDuplicates() {
        return includeDuplicates;
    }

    public void setMinimumMappingQuality(int minimumMappingQuality) {
        this.minimumMappingQuality = minimumMappingQuality;
    }

    public void setIncludeDuplicates(boolean includeDuplicates) {
        this.includeDuplicates = includeDuplicates;
    }

    public void setReadsToExclude(Set<String> readsToExclude) {
        this.readsToExclude = readsToExclude;
    }

    public void setPileupMode(int pileupMode){
        this.pileupMode = pileupMode;
    }

    public void setRefSequence(byte[] refSequence) {
        this.refSequence = refSequence;
    }

    private void createOrUpdatePileupsFromRead(SAMRecord read) throws IllegalArgumentException{
        byte[] bases = read.getReadBases();
        int refPos = read.getAlignmentStart();
        int readPos = 0;
        boolean isReadsFirstPileup = true;
        int alignedCount = 0;
        String sequenceName = read.getContig();
        List<CigarElement> cigarElementList = read.getCigar().getCigarElements();
        for (CigarElement cigarElement : cigarElementList) {
            CigarOperator op = cigarElement.getOperator();
            int length = cigarElement.getLength();
            if(CigarOperator.M == op || CigarOperator.EQ == op || CigarOperator.X == op){
                for (int i = 0; i < length; i++) {
                    int pileupPos = refPos + i;
                    int pileUpReadPos = readPos + i;
                    char readBaseUpper = Character.toUpperCase((char) bases[pileUpReadPos]);
                    if (!ALPHABET.contains(readBaseUpper)) {
                        throw new IllegalArgumentException("Invalid nucleotide detected: " + readBaseUpper + ". Expected one of: A, G, C, T, N.");
                    }
                    LocusPileUp pileup = cache.get(pileupPos);
                    if(pileup == null) {
                        if(pileupMode == CLAIM_ALL_READS_PILEUP_MODE) pileup = new LocusPileUp(sequenceName, pileupPos);
                        else pileup = new LocusPileUp(sequenceName, pileupPos, readQueue);
                        cache.putEntry(pileupPos, pileup);
                    }
                    PileupRead pileupRead = new PileupRead(read, pileUpReadPos, readBaseUpper, pileupPos);
                    if(isReadsFirstPileup){
                        isReadsFirstPileup = false;
                    }
                    else pileupRead.setStatus(PileupReadStatus.PILEUP_READ_STATUS_REPEATED);
                    if(pileupMode == CLAIM_NEW_AND_DIFFERING_READS_ONLY_PILEUP_MODE){
                        char referenceBaseUpper = Character.toUpperCase((char) refSequence[pileupPos-1]);
                        if(ALPHABET.contains(referenceBaseUpper)) pileupRead.setReferenceBase(referenceBaseUpper);
                    }
                    pileup.addPileupRead(pileupRead);
                    alignedCount++;
                }
            }
            if(op.consumesReferenceBases()){
                refPos += length;
            }
            if(op.consumesReadBases()){
                readPos += length;
            }
        }
    }

    @Override
    public Iterator<LocusPileUp> iterator() {
        return new Iterator<LocusPileUp>() {

            final SAMRecordIterator readIterator = samReaderStream.query(sequenceName, start, end, false);

            SAMRecord nextRead = readIterator.hasNext() ? readIterator.next() : null;
            LocusPileUp next = getNext();

            @Override
            public boolean hasNext() {
                return next != null;
            }

            @Override
            public LocusPileUp next() {
                if(next == null) throw new NoSuchElementException();
                LocusPileUp current = next;
                next = getNext();
                return current;
            }

            private LocusPileUp getNext() {
                LocusPileUp currentPileup;
                // Process pileups from reads and retrieve them when they are not covered by the latest read
                while (nextRead != null){
                    if(nextReferencePosition < nextRead.getStart()) {
                        currentPileup = cache.get(nextReferencePosition);
                        nextReferencePosition++;
                        if(currentPileup == null) continue;
                        cache.remove(nextReferencePosition-1);
                        readQueue.flushUnclaimedReads();
                        return currentPileup;
                    }
                    if(passesFilters(nextRead)){
                        createOrUpdatePileupsFromRead(nextRead);
                    }
                    if(readIterator.hasNext()){
                        nextRead = readIterator.next();
                    }
                    else nextRead = null;
                }
                // Retrieve remaining pileups
                while(nextReferencePosition <= end){
                    currentPileup = cache.get(nextReferencePosition);
                    nextReferencePosition++;
                    if(currentPileup != null) return currentPileup;
                }
                return null;
            }
        };
    }

    public boolean passesFilters(SAMRecord read){
        if(readsToExclude.contains(read.getReadName())) return false;
        if(read.getMappingQuality() < minimumMappingQuality) return false;
        if(read.getReadUnmappedFlag()) return false;
        if(!includeDuplicates && read.getDuplicateReadFlag()) return false;
        return true;
    }

    @Override
    public void forEach(Consumer<? super LocusPileUp> action) {
        Iterable.super.forEach(action);
    }

    /**
     * OnPileupQueue is a specialized implementation that extends a LinkedList of SAMRecord
     * and implements the Queue interface. This class is used to manage a queue of
     * SAMRecords for pileup generation while providing additional functionalities
     * for preventing duplicate reads and handling reads in a specific genomic context.
     */
    public class OnPileupQueue extends LinkedList<PileupRead> implements Queue<PileupRead> {

        private static final int DEFAULT_MAX_SIZE_LIMIT = 10_000_000;

        public OnPileupQueue(){
            super();
        }

        public boolean offerNew(PileupRead read){
            //Avoid memory issues when PileupReads from older pileups have not been claimed
            if(this.size() > DEFAULT_MAX_SIZE_LIMIT){
                flushUnclaimedReads();
            }
            PileupReadStatus status = updatePileupReadStatus(read);
            if(status != null){
                read.setStatus(status);
                return super.offer(read);
            }
            return false;
        }

        /**
         * Identifies and retrieves a list of SAMRecord objects from the queue that overlap a specified reference position.
         * Reads are considered overlapping if their start position is ≤ the pileup position and their end position is ≥ the pileup reference position.
         * The overlapping reads are removed from the queue and returned as part of the resulting list.
         * @return A list of SAMRecord objects that overlap the specified reference position. Each SAMRecord will no longer be in the queue.
         */
        public List<PileupRead> claimReadsOnPileup(LocusPileUp pileup){
            List<PileupRead> answer = new ArrayList<>();
            PileupRead nextRead = this.peek();
            while (!this.isEmpty() && overlap(pileup, nextRead)){
                answer.add(this.poll());
                nextRead = this.peek();
            }
            return answer;
        }

        /**
         * Removes and discards reads from the queue that have a start position less than the specified pileup position.
         * This is used to ensure that stale or unclaimed reads, which are no longer relevant
         * based on the current pileup position, are cleared from the queue.
         * @implicitParam nextReferencePosition the pileup position used as a threshold. Reads with a start position less than this
         *                       value will be removed from the queue if previously unclaimed.
         */
        public void flushUnclaimedReads(){
            this.removeIf(nextRead -> nextRead.getEnd() < nextReferencePosition-1);
        }
    }

    /**
     * Determines the pileup status of a PileupRead based on its previous status and whether it differs from the reference.
     * @param read    The PileupRead to evaluate.
     * @return The corresponding PileupReadStatus for the read, or null if no valid status applies.
     */
    private PileupReadStatus updatePileupReadStatus(PileupRead read) {
        boolean isNew = read.getStatus() == PileupReadStatus.PILEUP_READ_STATUS_UNKNOWN;
        if (pileupMode == CLAIM_NEW_READS_ONLY_PILEUP_MODE) {
            return isNew ?  PileupReadStatus.PILEUP_READ_STATUS_NEW
                    : null; // No status for repeated reads in this mode, will not be returned
        }
        else if (pileupMode == CLAIM_NEW_AND_DIFFERING_READS_ONLY_PILEUP_MODE) {
            if (isNew) {
                return read.differsFromReferenceAtPileup()
                        ? PileupReadStatus.PILEUP_READ_STATUS_NEW_AND_DIFFERING_BASE
                        : PileupReadStatus.PILEUP_READ_STATUS_NEW;
            }
            else {
                return read.differsFromReferenceAtPileup()
                        ? PileupReadStatus.PILEUP_READ_STATUS_REPEATED_AND_DIFFERING_BASE
                        : null; // No status for repeated reads without differences in this mode, will not be returned
            }
        }
        return null;
    }
}