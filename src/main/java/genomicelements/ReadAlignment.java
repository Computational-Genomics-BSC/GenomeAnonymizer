package genomicelements;

import htsjdk.samtools.Cigar;

public interface ReadAlignment {

    String getReadAlignmentId();
    String getReadName();
    String getSequenceName();
    int getStart();
    int getEnd();
    int getLength();
    byte[] getSequenceArray();
    byte[] getQualitiesArray();
    String getSequence();
    String getQualities();
    int getMappingQuality();
    Cigar getCigar();
    boolean isSupplementary();
    boolean isSecondary();


}
