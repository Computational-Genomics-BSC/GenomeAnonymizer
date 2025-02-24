package genomicelements;

import genomicelements.LocusPileupIterator.OnPileupQueue;

import java.util.List;

/**
 * The LocusPileUp class represents a specific genomic location and maintains a collection of alignments
 * (SAMRecord objects) that overlap at the given locus. It is designed to store and manage piled-up reads,
 * their corresponding bases, and positions in the reads. It also interacts with an OnPileupQueue to claim
 * new alignments that overlap the current pileup location.
 * @author Nicolas Gaitan
 */
public class LocusPileUp implements GenomicRegion{
    //TODO: Implement on demand read return, if it hasn´t been seen, or if at the pileup it has a change in base
    
    private String sequenceName;
    private int start;
    private int sequenceIdx;
    private int size;
    //Alignments, their bases and position in read are provided in lists with corresponding indexes
    private OnPileupQueue readQueue;

    public LocusPileUp(String sequenceName, int start, OnPileupQueue readQueue) {
        this.sequenceName = sequenceName;
        this.start = start;
        this.readQueue = readQueue;
    }

    public void addPileupRead(PileupRead read){
        size++;
        readQueue.offerNew(read);
    }

    /**
     * Retrieves and claims only new SAMRecord objects from the associated OnPileupQueue that overlap the current pileup position.
     * Overlapping reads are removed from the queue and returned.
     * @return A list of SAMRecord objects representing the newly claimed reads overlapping the pileup position.
     */
    public List<PileupRead> claimNewReadsOnPileup(){
        return readQueue.claimNewReadsOnPileup(this);
    }

    @Override
    public String getSequenceName() {
        return sequenceName;
    }

    @Override
    public int getSequenceIdx() {
        return sequenceIdx;
    }

    public int getLocation(){
        return start;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getEnd() {
        return start;
    }

    public int size() {
        return size;
    }

    @Override
    public void setSequenceIdx(int sequenceIdx) {
        this.sequenceIdx = sequenceIdx;
    }
}
