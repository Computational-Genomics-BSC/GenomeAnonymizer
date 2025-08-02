package genomicelements;

import static utils.Operations.compare;

/**
 * Represents a signal from a variation at a genomic location in a given read.
 * @author Nicolas Gaitan
 */
public class Signal implements GenomicRegion, ReadSignal, Cloneable{

    private String sequenceName;
    private int location;
    private int sequenceIdx;
    private String readAlnName;
    private int inReadPosition;
    private int length;
    private Source source;
    private boolean isFromNormalDataset;
    private boolean isGermline = false;
    private boolean isHandled = false; // Indicates if the signal has been tested for classification
    private boolean classifiedByDistance = false;
    private byte[] altAllele;
    private byte[] refAllele;
    private IndelSignalType indelSignalType = IndelSignalType.NOT_INDEL;


    public Signal(String sequenceName, int location, String readAlnName, int readPosition, int length, Source source) {
        this.sequenceName = sequenceName;
        this.location = location;
        this.readAlnName = readAlnName;
        this.inReadPosition = readPosition;
        this.length = length;
        this.source = source;
        //this.isFromNormalDataset = isFromNormalDataset;
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

    public int getStart() {
        return location;
    }

    public int getEnd() {
        return location;
    }

    public int getSequenceIdx() {
        return sequenceIdx;
    }

    public int getLength() {
        return length;
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

    public boolean isSNV() {
        return source == Source.SNV;
    }

    public boolean isIndel() {
        return source == Source.INDEL;
    }

    public boolean isSoftClip() {
        return source == Source.SOFT_CLIP;
    }

    public boolean isDestructiveSignal() {
        return source == Source.INSERT_SIZE || source == Source.CHROM_CHANGE;
    }

    public boolean isIntraAlignmentSignal() {
        return source == Source.SNV || source == Source.INDEL || source == Source.SOFT_CLIP;
    }

    public boolean isInterAlignmentSignal() {
        return source == Source.INSERT_SIZE || source == Source.STRAND_ORIENTATION
                || source == Source.CHROM_CHANGE;
    }

    public boolean isFromNormalDataset() {
        return isFromNormalDataset;
    }

    public boolean isFromTumoralDataset() {
        return !isFromNormalDataset;
    }

    public boolean isGermline() {
        return isGermline;
    }

    public boolean isClassifiedByDistance() {
        return classifiedByDistance;
    }

    public byte[] getAltAllele() {
        return altAllele;
    }

    public byte[] getRefAllele() {
        return refAllele;
    }

    public void setSequenceIdx(int sequenceIdx) {
        this.sequenceIdx = sequenceIdx;
    }

    public IndelSignalType getIndelSignalType() {
        return indelSignalType;
    }

    public void setIsFromNormalDataset(boolean isNormal){
        this.isFromNormalDataset = isNormal;
    }

    public void setIsGermline(boolean isGermline) {
        this.isGermline = isGermline;
    }

    public void setClassifiedByDistance(boolean classifiedByDistance) {
        this.classifiedByDistance = classifiedByDistance;
    }

    public void setAltAllele(byte[] altAllele) {
        this.altAllele = altAllele;
    }

    public void setRefAllele(byte[] refAllele) {
        this.refAllele = refAllele;
    }

    public void setIndelSignalType(IndelSignalType indelSignalType) {
        this.indelSignalType = indelSignalType;
    }

    public void setReadAlnName(String readAlnName) {
        this.readAlnName = readAlnName;
    }

    public void setInReadPosition(int inReadPosition) {
        this.inReadPosition = inReadPosition;
    }

    public boolean locatedAtReadStart(){
        return inReadPosition == 0;
    }

    public boolean locatedAtReadEnd(){
        return !locatedAtReadStart();
    }
    /**
     * Generates a unique key for the signal based on its properties.
     * This key can be used to identify the signal in a collection or database.
     * @return A string representing the unique key for the signal.
     */
    public String getSignalKey(){
        if (isInterAlignmentSignal()) return getReadAlnName();
        return this.getSequenceName() + ":" + this.getLocation() + ":" + this.getLength() +
                ":" + this.getSource() + ":" + this.getIndelSignalType();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Signal)) return false;

        Signal signal = (Signal) o;

        if (location != signal.location) return false;
        if (sequenceIdx != signal.sequenceIdx) return false;
        if (inReadPosition != signal.inReadPosition) return false;
        if (length != signal.length) return false;
        if (isFromNormalDataset != signal.isFromNormalDataset) return false;
        if (!sequenceName.equals(signal.sequenceName)) return false;
        if (!readAlnName.equals(signal.readAlnName)) return false;
        return source == signal.source;
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

    public String toAbridgedString() {
        return this.getLocation() + ":" + this.getLength() + ":" + this.getSource() + ":" +
                (this.isFromNormalDataset() ? "normal" : "tumor");
    }

    @Override
    public int compareTo(GenomicRegion genomicRegion) {
        return compare(this, genomicRegion);
    }

    /**
     * Creates a clone of this Signal with new readAlnName and inReadPosition values.
     * All other fields are shared with the original instance.
     *
     * @return A new Signal instance with the specified fields
     */
    @Override
    public Object clone() throws CloneNotSupportedException {
        Signal answer = (Signal) super.clone();
        answer.readAlnName = this.readAlnName; // This is a new reference, not shared
        answer.inReadPosition = this.inReadPosition;
        return answer;
    }

    public boolean isHandled() {
        return isHandled;
    }

    public void setHandled(boolean handled) {
        isHandled = handled;
    }

    public enum Source{
        SNV(0, "SNV"),
        INDEL(1, "INDEL"),
        SOFT_CLIP(2, "SOFT_CLIP"),
        //HARD_CLIP(),
        INSERT_SIZE(3, "INSERT_SIZE"),
        STRAND_ORIENTATION(4, "STRAND_ORIENTATION"),
        CHROM_CHANGE(5, "CHROM_CHANGE");

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

    public enum IndelSignalType {
        NOT_INDEL,
        INSERTION,
        DELETION,
    }
}
