package genomicelements;

/**
 * Generic representation of a simple genomic region defined by coordinates on a reference
 * @author Nicolas Gaitan
 */
public class GenomicRegionBaseImpl implements GenomicRegion{
    public String sequenceName;
    //Represents the order of a sequence relative to the reference genome
    public int sequenceIdx;
    public int start;
    public int end;


    public GenomicRegionBaseImpl(String sequenceName, int start, int end) {
        this.sequenceName = sequenceName;
        this.start = start;
        this.end = end;
    }

    public GenomicRegionBaseImpl(String sequenceName, int sequenceIdx, int start, int end) {
        this.sequenceName = sequenceName;
        this.sequenceIdx = sequenceIdx;
        this.start = start;
        this.end = end;
    }

    public String getSequenceName() {
        return sequenceName;
    }

    public void setSequenceName(String sequenceName) {
        this.sequenceName = sequenceName;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public int getSequenceIdx() {
        return sequenceIdx;
    }

    public void setSequenceIdx(int sequenceIdx) {
        this.sequenceIdx = sequenceIdx;
    }

    public String toString(){
        return sequenceName + ":" + start + "-" + end;
    }
}
