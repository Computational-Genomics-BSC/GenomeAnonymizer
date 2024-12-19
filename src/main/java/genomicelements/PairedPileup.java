package genomicelements;

import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;

import java.util.List;

/**
 * Representation of a tuple of read alignment pileups, with one coming from a normal sample and the other from the tumor pair
 * allowing deep comparison between both mappings
 * @author Nicolas Gaitan
 */
public class PairedPileup {

    private LocusInfo normalLocus;
    private LocusInfo tumorLocus;
    private String refenceSequenceName;
    private int referenceSequenceIdx;
    private int referencePos;

    public PairedPileup(LocusInfo normalLocus, LocusInfo tumorLocus) {
        // DEBUG
        assert normalLocus.getPosition() == tumorLocus.getPosition(): "Normal and tumor locuses are not in the same position " +
                normalLocus.getPosition() + tumorLocus.getPosition();
        // DEBUG
        this.normalLocus = normalLocus;
        this.tumorLocus = tumorLocus;
        refenceSequenceName = normalLocus.getSequenceName();
        referenceSequenceIdx = normalLocus.getSequenceIndex();
        referencePos = normalLocus.getPosition();
    }

    public PairedPileup(LocusInfo normalLocus) {
        this.normalLocus = normalLocus;
        this.tumorLocus = null;
        refenceSequenceName = normalLocus.getSequenceName();
        referenceSequenceIdx = normalLocus.getSequenceIndex();
        referencePos = normalLocus.getPosition();
    }

    public List<RecordAndOffset> getNormalPileup(){
        return normalLocus.getRecordAndOffsets();
    }

    public List<RecordAndOffset> getTumorPileup(){
        return tumorLocus == null ? null : tumorLocus.getRecordAndOffsets();
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
