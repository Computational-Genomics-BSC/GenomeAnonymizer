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

    public ShortReadAlignment(SAMRecord samRecord){
        super(samRecord);
        setReadAlignmentId(generateReadId(samRecord));
    }

    public boolean isPair1(){
        return record.getFirstOfPairFlag();
    }

    public boolean isPair2(){
        return record.getSecondOfPairFlag();
    }

    public int getPairIdx() {
        return isPair1() ? PAIR_1_IDX : PAIR_2_IDX;
    }

    public boolean isReverse() {
        return record.getReadNegativeStrandFlag();
    }

    public void setReverse(boolean reverse) {
        record.setReadNegativeStrandFlag(reverse);
    }

    public void setPairIdx(int pairIdx) {
        if(pairIdx == PAIR_1_IDX){
            record.setFirstOfPairFlag(true);
            record.setSecondOfPairFlag(false);
        }else if(pairIdx == PAIR_2_IDX){
            record.setFirstOfPairFlag(false);
            record.setSecondOfPairFlag(true);
        }
    }

    @Override
    public String toString(){
        return "Read name=" + this.getReadName() + " ReadAlignmentID=" + this.getReadName() + " seq=" + this.getSequenceName() + " pos=" + this.getStart() + " end=" +
                this.getEnd() + " length=" + this.getLength() + " sequence=" +  new String(this.getSequenceArray(), StandardCharsets.UTF_8)
                + " qualities=" + Arrays.toString(this.getQualitiesArray()) + " pairIdx=" + this.getPairIdx() +
                " isSupplementaryOrSecondary=" + this.isSupplementary() + " isReverse=" + this.isReverse();
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
