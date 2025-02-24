package genomicelements;

import htsjdk.samtools.*;
import htsjdk.samtools.SAMRecord;
import utils.MapCacheFIFO;

import java.util.*;
import java.util.function.Consumer;

import genomicelements.PileupRead.PileupReadStatus;

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
    public static final int CLAIM_NEW_AND_DIFFERING_READS_PILEUP_MODE = 1;

    private SAMRecordIterator readIterator;
    private MapCacheFIFO<Integer, LocusPileUp> cache;
    private int start;
    private int end;
    private int nextReferencePosition;
    private int minimumMappingQuality = 0;
    private boolean includeDuplicates = false;
    private OnPileupQueue readQueue;
    private Set<String> readsToExclude;
    private int pileupMode = 0;
    //This reference sequence is 0-based, given it comes in a byte array
    private byte[] refSequence;

    //TODO: Handle pileup when traversing different chromosomes or specify it cannot be used for that

    public LocusPileupIterator(SamReader samReaderStream, String sequenceName, int start, int end) throws IllegalStateException{
        setPileupMode(CLAIM_ALL_READS_PILEUP_MODE);
        init(samReaderStream.query(sequenceName, start, end, false), start, end);
    }

    public LocusPileupIterator(SamReader samReaderStream, String sequenceName, int start, int end, int pileupMode, byte[] refSequence) throws IllegalStateException{
        init(samReaderStream.query(sequenceName, start, end, false), start, end);
        setPileupMode(pileupMode);
        setRefSequence(refSequence);
        if(this.pileupMode != CLAIM_NEW_AND_DIFFERING_READS_PILEUP_MODE && pileupMode != CLAIM_ALL_READS_PILEUP_MODE) {
            throw new IllegalArgumentException("Invalid pileup mode: " + pileupMode + ". Valid modes are: " +
                    "0 - Claim all reads pileup mode, 1 - Claim new and different reads pileup mode.");
        }
        if(this.pileupMode == CLAIM_NEW_AND_DIFFERING_READS_PILEUP_MODE && refSequence == null) {
            throw new IllegalStateException("Reference sequence must be provided when using the new and different reads pileup mode.");
        }
    }

    public LocusPileupIterator(SAMRecordIterator readIterator, int start, int end) throws IllegalStateException{
        init(readIterator, start, end);
    }

    public void init(SAMRecordIterator readIterator, int start, int end) throws IllegalStateException{
        this.readIterator = readIterator;
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

    //TODO: Implement on demand read return, if it hasn´t been seen, or if at the pileup it has a change in base

    private void createOrUpdatePileupsFromRead(SAMRecord read) throws IllegalArgumentException{
        byte[] bases = read.getReadBases();
        int refPos = read.getAlignmentStart();
        int readPos = 0;
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
                        pileup = new LocusPileUp(sequenceName, pileupPos, readQueue);
                        cache.putEntry(pileupPos, pileup);
                    }
                    PileupRead pileupRead = new PileupRead(read, pileUpReadPos, readBaseUpper, pileupPos);
                    if(pileupMode == CLAIM_NEW_AND_DIFFERING_READS_PILEUP_MODE){
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
        if(!readIterator.hasNext()) return Collections.emptyIterator();
        return new Iterator<LocusPileUp>() {

            // To make sure the readIterator is not null, it is checked in the Iterable class
            SAMRecord nextRead = readIterator.next();
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
        if (readsToExclude.contains(read.getReadName())) return false;
        if(read.getReadUnmappedFlag()) return false;
        if(read.getMappingQuality() < minimumMappingQuality) return false;
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

        private Set<String> addedReadIds;

        public OnPileupQueue(){
            super();
            addedReadIds = new HashSet<>();
        }

        public boolean offerNew(PileupRead read){
            //Avoid memory issues when PileupReads from older pileups have not been claimed
            if(this.size() > DEFAULT_MAX_SIZE_LIMIT){
                flushUnclaimedReads();
            }
            String readAlnId = read.getReadAlignmentId();
            PileupReadStatus status = determinePileupReadStatus(read, addedReadIds.contains(readAlnId));
            if(status != null){
                read.setStatus(status);
                if(status == PileupReadStatus.PILEUP_READ_STATUS_NEW || status == PileupReadStatus.PILEUP_READ_STATUS_NEW_AND_DIFFERING_BASE) {
                    addedReadIds.add(readAlnId);
                }
                return super.offer(read);
            }
            return false;
        }
        
        
        
        /**
         * Identifies and retrieves a list of SAMRecord objects from the queue that overlap a specified reference position.
         * Reads are considered overlapping if their start position is ≤ the pileup position and their end position is ≥ the pileup reference position.
         * The overlapping reads are removed from the queue and returned as part of the resulting list.
         * @implicitParam nextReferencePosition the pileup position used as a threshold.
         * @return A list of SAMRecord objects that overlap the specified reference position. Each SAMRecord will no longer be in the queue.
         */
        public List<PileupRead> claimNewReadsOnPileup(LocusPileUp pileup){
            List<PileupRead> answer = new ArrayList<>();
            PileupRead nextRead = this.peek();
            while (!this.isEmpty() &&
                    (nextRead.getStart() <= pileup.getLocation() && nextRead.getEnd() >= pileup.getLocation())){
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
            PileupRead nextRead = this.peek();
            while (!this.isEmpty() && nextRead.getEnd() < nextReferencePosition-1){
                this.poll();
                addedReadIds.remove(nextRead.getReadAlignmentId());
                nextRead = this.peek();
            }
        }
    }

    /**
     * Determines the pileup status of a PileupRead based on its presence in the set and whether it differs from the reference.
     *
     * @param read    The PileupRead to evaluate.
     * @param isInSet A boolean indicating whether the read is already in the addedReadIds set.
     * @return The corresponding PileupReadStatus for the read, or null if no valid status applies.
     */
    private PileupReadStatus determinePileupReadStatus(PileupRead read, boolean isInSet) {
        if (pileupMode == CLAIM_ALL_READS_PILEUP_MODE) {
            return isInSet ? PileupReadStatus.PILEUP_READ_STATUS_REPEATED : PileupReadStatus.PILEUP_READ_STATUS_NEW;
        }
        else if (pileupMode == CLAIM_NEW_AND_DIFFERING_READS_PILEUP_MODE) {
            if (!isInSet) {
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
