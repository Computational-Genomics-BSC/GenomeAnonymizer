package genomicelements;

import genomicelements.LocusPileupIterator.OnPileupQueue;

import java.util.ArrayList;
import java.util.List;

import static genomicelements.LocusPileupIterator.DEFAULT_PILEUP_READ_LIMIT;

/**
 * The LocusPileUp class represents a specific genomic location and maintains a collection of alignments
 * (SAMRecord objects) that overlap at the given locus. It is designed to store and manage piled-up reads,
 * their corresponding bases, and positions in the reads. It also interacts with an OnPileupQueue to claim
 * new and/or different alignments at this pileup location.
 * @author Nicolas Gaitan
 */
public class LocusPileUp implements GenomicRegion{

    private String sequenceName;
    private int start;
    private int sequenceIdx;
    private int size = 0;
    private List<PileupRead> pileupReads;
    //Limit to the number of reads that can be piled up at a given locus
    private int pileupReadLimit = DEFAULT_PILEUP_READ_LIMIT;
    //Alignments, their bases and position in read are provided in lists with corresponding indexes
    private OnPileupQueue readQueue = null;

    public LocusPileUp(String sequenceName, int start) {
        this.sequenceName = sequenceName;
        this.start = start;
        this.pileupReads = new ArrayList<>();
    }

    public LocusPileUp(String sequenceName, int start, OnPileupQueue readQueue) {
        this.sequenceName = sequenceName;
        this.start = start;
        this.readQueue = readQueue;
    }

    public void addPileupRead(PileupRead read){
        size++;
        //If a valid OnPileupQueue is provided when instantiated, the PileupRead is added to the queue
        // (modes: CLAIM_NEW_READS_ONLY_PILEUP_MODE, CLAIM_NEW_AND_DIFFERING_READS_ONLY_PILEUP_MODE)
        if(readQueue != null){
            readQueue.offerNew(read);
        }
        else{
            //If the pileupReadLimit is reached, the read is not added to the pileupReads list
            // (mode: CLAIM_ALL_READS_PILEUP_MODE)
            if(pileupReads.size() <= pileupReadLimit){
                pileupReads.add(read);
            }
        }
    }

    public List<PileupRead> getPileupReads(){
        return pileupReads;
    }

    /**
     * Retrieves and claims only new SAMRecord objects from the associated OnPileupQueue that overlap the current pileup position,
     *  or those that differ at this pileup position.
     * @return A list of SAMRecord objects representing the newly claimed reads overlapping the pileup position.
     */
    public List<PileupRead> claimReadsOnPileup(){
        return readQueue.claimReadsOnPileup(this);
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

    public int getPileupReadLimit() {
        return pileupReadLimit;
    }

    public int size() {
        return size;
    }

    public void setPileupReadLimit(int pileupReadLimit) {
        this.pileupReadLimit = pileupReadLimit;
    }

    @Override
    public void setSequenceIdx(int sequenceIdx) {
        this.sequenceIdx = sequenceIdx;
    }
}