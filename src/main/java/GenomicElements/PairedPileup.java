package GenomicElements;

import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;

import java.util.List;

public class PairedPileup {

    public static final String TUMORAL_DATASET = "TUMORAL_DATASET";
    public static final String NORMAL_DATASET = "NORMAL_DATASET";

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

    public List<RecordAndOffset> getNormalPileup(){
        return normalLocus.getRecordAndOffsets();
    }

    public List<RecordAndOffset> getTumorPileup(){
        return tumorLocus.getRecordAndOffsets();
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

}
