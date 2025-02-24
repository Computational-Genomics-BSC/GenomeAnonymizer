package genomicelements;

public interface GenomicRegion {
    public String getSequenceName();
    public int getSequenceIdx();
    public int getStart();
    public int getEnd();
    public void setSequenceIdx(int sequenceIdx);
    public String toString();
}
