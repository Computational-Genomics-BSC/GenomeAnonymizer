package genomicelements;

/**
 * Represents a signal from a variation at a genomic location in a given read.
 * @author Nicolas Gaitan
 */
public class Signal {

    private String sequenceName;
    private int location;
    private int sequenceIdx;
    private String readAlnName;
    private int inReadPosition;
    private int length;
    private byte[] sequenceBytes;
    private Source source;
    private boolean isFromNormalDataset;
    private PairCalledVariation pairCalledVariation = null;


    public Signal(String sequenceName, int location, String readAlnName, int readPosition, int length, Source source) {
        this.sequenceName = sequenceName;
        this.location = location;
        this.readAlnName = readAlnName;
        this.inReadPosition = readPosition;
        this.length = length;
        this.source = source;
        //this.isFromNormalDataset = isFromNormalDataset;
    }

    /**
     * Generate a simple signal from a PairCalledVariation, this will not have the information about the dataset it comes from
     * @param readAlnName
     * @param pairCalledVariation
     */
    public Signal(String readAlnName, PairCalledVariation pairCalledVariation) {
        this(pairCalledVariation.getSeqName(), pairCalledVariation.getPos(), readAlnName, pairCalledVariation.getInReadPosition(readAlnName),
                pairCalledVariation.getLength(), Source.SIMPLE_VARIATION);
        this.pairCalledVariation = pairCalledVariation;
    }

    public String getSequenceName() {
        return sequenceName;
    }

    /**
     * @return Returns the reference coordinate of the closest aligned base pair to the signal
     *  e.g. For soft clips, if it is at the start of a read, returns the start position of the mapped read,
     *  if it is at the end, returns the position of the alignment end
     */
    public int getLocation() {
        return location;
    }

    public int getSequenceIdx() {
        return sequenceIdx;
    }

    public int getLength() {
        return length;
    }

    public byte[] getSequenceBytes(){
        return sequenceBytes;
    }

    public Source getSource() {
        return source;
    }

    public String getReadAlnName() {
        return readAlnName;
    }

    public int getInReadPosition() {
        return inReadPosition;
    }

    public boolean isFromNormalDataset() {
        return isFromNormalDataset;
    }

    public boolean isFromTumoralDataset() {
        return !isFromNormalDataset;
    }

    public PairCalledVariation getCalledVariation() {
        return pairCalledVariation;
    }

    public void setSequenceIdx(int sequenceIdx) {
        this.sequenceIdx = sequenceIdx;
    }

    public void setSequenceBytes(byte[] seqBytes){
        this.sequenceBytes = seqBytes;
    }

    public void setIsFromNormalDataset(boolean isNormal){
        this.isFromNormalDataset = isNormal;
    }

    public boolean locatedAtReadStart(){
        return inReadPosition == 0;
    }

    public boolean locatedAtReadEnd(){
        return !locatedAtReadStart();
    }

    @Override
    public String toString() {
        return "Signal{" +
                "readAlnName='" + readAlnName + '\'' +
                ", length=" + length +
                ", source=" + source +
                ", isFromNormalDataset=" + isFromNormalDataset +
                ", sequenceName='" + sequenceName + '\'' +
                ", location=" + location +
                ", inReadPosition=" + inReadPosition +
                '}';
    }

    public enum Source{
        SIMPLE_VARIATION(1, "DEFINABLE_VARIATION"),
        SOFT_CLIP(2, "SOFT_CLIP"),
        //HARD_CLIP(),
        INSERT_SIZE(3, "INSERT_SIZE");

        private final int value;
        private final String name;

        Source(int value, String name){
            this.name = name;
            this.value = value;
        }

        public int getValue() {
            return value;
        }

        public String getName() {
            return name;
        }
    }
}
