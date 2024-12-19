package genomicelements;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.List;

/**
 * Represents a short read alignment
 * @author Nicolas Gaitan
 */
public class ShortReadAlignment extends ReadAlignmentBaseImpl implements ReadAlignment{

    public static final int PAIR_1_IDX = 0;
    public static final int PAIR_2_IDX = 1;
    public static final String DEFAULT_ID_NAME_SEPARATOR = ";";

    private int pairIdx;
    private boolean isReverse;
    private List<SAMRecord.SAMTagAndValue> tags;

    public ShortReadAlignment(String readName, String sequenceName, int flags, int start, int end, boolean isSupplementary, int mapQ,
                              Cigar cigar, boolean isFirstOfPair, boolean isReverse) {
        super(readName, sequenceName, flags, start, end, isSupplementary, mapQ, cigar);
        this.setReadAlignmentId(generateReadAlignmentId(readName));
        this.pairIdx = isFirstOfPair ? PAIR_1_IDX : PAIR_2_IDX;
        this.isReverse = isReverse;
    }

    public ShortReadAlignment(SAMRecord samRecord){
        this(samRecord.getReadName(), samRecord.getContig(), samRecord.getFlags(), samRecord.getStart(), samRecord.getEnd(), samRecord.isSecondaryOrSupplementary(),
                samRecord.getMappingQuality(), samRecord.getCigar(), samRecord.getFirstOfPairFlag(), samRecord.getReadNegativeStrandFlag());
        this.setSequenceArray(samRecord.getReadBases());
        this.setQualitiesArray(samRecord.getBaseQualities());
        this.setReadAlignmentId(generateReadId(samRecord));
        this.setHeader(samRecord.getHeader());
        this.setTags(samRecord.getAttributes());

    }


    public boolean isPair1(){
        return pairIdx == PAIR_1_IDX;
    }

    public boolean isPair2(){
        return pairIdx == PAIR_2_IDX;
    }

    public int getPairIdx() {
        return pairIdx;
    }

    public boolean isReverse() {
        return isReverse;
    }

    public void setReverse(boolean reverse) {
        isReverse = reverse;
    }

    public void setPairIdx(int pairIdx) {
        this.pairIdx = pairIdx;
    }

    public List<SAMRecord.SAMTagAndValue> getTags(){
        return tags;
    }

    public void setTags(List<SAMRecord.SAMTagAndValue> tags){
        this.tags = tags;
    }

    @Override
    public String toString(){
        return "Read name=" + this.getReadName() + " ReadAlignmentID=" + this.getReadName() + " seq=" + this.getSequenceName() + " pos=" + this.getStart() + " end=" +
                this.getEnd() + " length=" + this.getLength() + " sequence=" +  new String(this.getSequenceArray(), StandardCharsets.UTF_8)
                + " qualities=" + Arrays.toString(this.getQualitiesArray()) + " pairIdx=" + this.getPairIdx() +
                " isSupplementaryOrSecondary=" + this.isSupplementary() + " isAnonymized=" + " isReverse=" + this.isReverse;
    }


    //Provides a unique ID for each read alignment
    public static String generateReadId(String readName, int pairIdx, int refPos){
        return readName + DEFAULT_ID_NAME_SEPARATOR + pairIdx + DEFAULT_ID_NAME_SEPARATOR + refPos;
    }

    //Provides a unique ID for each read alignment, directly from a SAMRecord
    public static String generateReadId(SAMRecord alignment){
        int pairIdx = alignment.getFirstOfPairFlag() ? PAIR_1_IDX : PAIR_2_IDX;
        String baseName = generateReadId(alignment.getReadName(), pairIdx, alignment.getAlignmentStart());
        String complement = generateAlignmentHash(alignment);
        return baseName + DEFAULT_ID_NAME_SEPARATOR + complement;
    }

    public static String generateAlignmentHash(SAMRecord samRec) {
        long answer = 17;
        answer = 37*answer + samRec.getCigar().toString().hashCode();
        answer = 37*answer + samRec.getBaseQualityString().hashCode();
        answer = 37*answer + samRec.getReadString().hashCode();
        return String.valueOf(answer);
    }

}
