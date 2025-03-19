package genomicelements;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;

import java.util.Comparator;

import static genomicelements.ShortAnonymizedReadAlignment.generateReadAlnId;
import static utils.Operations.compare;

public class PileupRead implements GenomicRegion{
    private final int location;
    private SAMRecord read;
    private int readPosition;
    private char baseAtPileup;
    private char referenceBase;
    private PileupReadStatus status;
    private String readAlignmentId;

    private int sequenceIdx = 0;

    public PileupRead(SAMRecord read, int readPosition, char baseAtPileup, int location) {
        this.read = read;
        this.readPosition = readPosition;
        this.baseAtPileup = baseAtPileup;
        this.readAlignmentId = generateReadAlnId(read);
        this.location = location;
        setStatus(PileupReadStatus.PILEUP_READ_STATUS_UNKNOWN);
        this.referenceBase = 0;
    }

    public SAMRecord getRead() {
        return read;
    }

    public int getLocation() {
        return location;
    }

    public String getSequenceName(){
        return read.getReferenceName();
    }

    @Override
    public int getSequenceIdx() {
        return sequenceIdx;
    }

    public int getStart(){
        return read.getAlignmentStart();
    }

    public int getEnd(){
        return read.getAlignmentEnd();
    }

    @Override
    public void setSequenceIdx(int sequenceIdx) {
        this.sequenceIdx = sequenceIdx;
    }

    public String getReadName(){
        return read.getReadName();
    }

    public String getReadAlignmentId(){
        return readAlignmentId;
    }

    public byte[] getReadBases(){
        return read.getReadBases();
    }

    /**
     * Retrieves the 0-based read position corresponding to the current pileup.
     * @return the read position as an integer.
     */
    public int getReadPosition() {
        return readPosition;
    }

    public char getBaseAtPileup() {
        return baseAtPileup;
    }

    public Cigar getCigar() {
        return read.getCigar();
    }

    public char getReferenceBase() {
        return referenceBase;
    }

    public void setReferenceBase(char referenceBase) {
        this.referenceBase = referenceBase;
    }

    public boolean isNew(){
        return status == PileupRead.PileupReadStatus.PILEUP_READ_STATUS_NEW ||
                status == PileupRead.PileupReadStatus.PILEUP_READ_STATUS_NEW_AND_DIFFERING_BASE;
    }

    /**
     * Checks if the base at the pileup position differs from the corresponding reference base.
     * Does not inform about deletions or insertions, neither breakpoints
     * @return true if the base at the pileup position does not match the reference base;
     * otherwise or if it is not a certainty (e.g. N bases) false.
     */
    public boolean differsFromReferenceAtPileup() {
        if(baseAtPileup=='N') return false;
        if(referenceBase == 0) return false;
        return referenceBase != baseAtPileup;
    }

    /**
     * Retrieves the current status of the read in the context of the pileup.
     * @return The status of the pileup read, represented as a {@code PileupReadStatus} enum.
     * PILEUP_READ_STATUS_UNKNOWN: The pileup read status is unknown, only for cases where all PileupReads are returned
     * PILEUP_READ_STATUS_NEW: Pileup read that is claimed for the first time
     * PILEUP_READ_STATUS_REPEATED: Pileup read that has been claimed in a previous pileup
     * PILEUP_READ_STATUS_NEW_AND_DIFFERING_BASE: Pileup read that is new and differs from reference at this position
     * PILEUP_READ_STATUS_REPEATED_AND_DIFFERING_BASE: Previously claimed Pileup read which differs from reference at this position
     */
    public PileupReadStatus getStatus() {
        return status;
    }

    public void setStatus(PileupReadStatus status) {
        this.status = status;
    }

    @Override
    public int compareTo(GenomicRegion genomicRegion) {
        int cmp = compare(this, genomicRegion);
        if (cmp < -1 || cmp > 1) {
            return cmp;
        }
        if (genomicRegion instanceof PileupRead pileupRead2){
            return Comparator.comparingInt(PileupRead::getLocation)
                    .compare(this, pileupRead2);
        }
        return cmp;
    }

    /**
     * Enum representing the status of a read within a pileup.
     * This is used to classify the current state of a read in the context
     * of an alignment pileup at a specific genomic location.
     */
    public enum PileupReadStatus {
        PILEUP_READ_STATUS_UNKNOWN,
        PILEUP_READ_STATUS_NEW,
        PILEUP_READ_STATUS_REPEATED,
        PILEUP_READ_STATUS_NEW_AND_DIFFERING_BASE,
        PILEUP_READ_STATUS_REPEATED_AND_DIFFERING_BASE;
    }

    @Override
    public String toString() {
        return "PileupRead [read=" + read + ", readPosition=" + readPosition + " pileupLocus=" + location + ", baseAtPileup=" + baseAtPileup + " status=" + status + "]";
    }
}
