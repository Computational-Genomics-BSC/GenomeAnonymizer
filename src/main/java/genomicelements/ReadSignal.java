package genomicelements;

public interface ReadSignal extends GenomicRegion {
    String getSequenceName();

    int getLocation();

    int getStart();

    int getEnd();

    int getSequenceIdx();

    int getLength();

    String getReadAlnName();

    int getInReadPosition();

    boolean isFromNormalDataset();

    boolean isGermline();

    boolean isHandled();

    void setHandled(boolean handled);

    boolean isClassifiedByDistance();

    byte[] getAltAllele();

    boolean isSNV();

    boolean isIndel();

    boolean isSoftClip();

    boolean isDestructiveSignal();

    boolean isFromTumoralDataset();

    void setIsGermline(boolean isGermline);

    void setIsFromNormalDataset(boolean isFromNormalDataset);

    void setClassifiedByDistance(boolean classifiedByDistance);

    String getSignalKey();

    Signal.IndelSignalType getIndelSignalType();

    Signal.Source getSource();

    boolean isIntraAlignmentSignal();

    boolean isInterAlignmentSignal();

}
