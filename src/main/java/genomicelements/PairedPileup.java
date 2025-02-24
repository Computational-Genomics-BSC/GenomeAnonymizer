package genomicelements;

import genomicelements.LocusPileupIterator.OnPileupQueue;

/**
 * Representation of a tuple of read alignment pileups, with one coming from a normal sample and the other from the tumor pair
 * allowing deep comparison between both mappings
 * @author Nicolas Gaitan
 */
public class PairedPileup {

    private LocusPileUp normalLocus;
    private LocusPileUp tumorLocus;
    private String refenceSequenceName;
    private int referenceSequenceIdx;
    private int referencePos;

    public PairedPileup(LocusPileUp normalLocus, LocusPileUp tumorLocus) {
        this.normalLocus = normalLocus;
        this.tumorLocus = tumorLocus;
        refenceSequenceName = normalLocus.getSequenceName();
        referenceSequenceIdx = normalLocus.getSequenceIdx();
        referencePos = normalLocus.getLocation();
    }

    public PairedPileup(LocusPileUp normalLocus) {
        this.normalLocus = normalLocus;
        this.tumorLocus = null;
        refenceSequenceName = normalLocus.getSequenceName();
        referenceSequenceIdx = normalLocus.getSequenceIdx();
        referencePos = normalLocus.getLocation();
    }

    public LocusPileUp getNormalPileup(){
        return normalLocus;
    }

    public LocusPileUp getTumorPileup(){
        return tumorLocus == null ? null : tumorLocus;
    }

    public int totalSize(){
        return normalLocus.size() + tumorLocus.size();
    }

    public String getRefenceSequenceName() {
        return refenceSequenceName;
    }

    public int getReferenceSequenceIdx() {
        return referenceSequenceIdx;
    }

    public int getReferencePos() {
        return referencePos;
    }

    public boolean isOnlyNormal(){
        return tumorLocus == null;
    }
}
