package genomicelements;

/**
 * Representation of a tuple of read alignment pileups, with one coming from a normal sample and the other from the tumor pair
 * allowing deep comparison between both mappings
 * @author Nicolas Gaitan
 */
public class PairedPileup {

    private LocusPileUp normalLocus = null;
    private LocusPileUp tumorLocus = null;
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

    public PairedPileup(LocusPileUp locus, boolean isNormal) {
        if(isNormal) this.normalLocus = locus;
        else this.tumorLocus = locus;
        refenceSequenceName = locus.getSequenceName();
        referenceSequenceIdx = locus.getSequenceIdx();
        referencePos = locus.getLocation();
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

    public int getLocation() {
        return referencePos;
    }

    public boolean isOnlyNormal(){
        return tumorLocus == null;
    }
}
