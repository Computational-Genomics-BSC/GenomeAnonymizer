package genomicelements;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.List;

/**
 * Base implementation of a generic read alignment, only to be used for inheritance
 * @author Nicolas Gaitan
 */
public class ReadAlignmentBaseImpl implements ReadAlignment{

    private String readAlignmentId;
    protected SAMRecord record;

    public ReadAlignmentBaseImpl(SAMRecord record) {
        this.record = record;
        this.readAlignmentId = record.getReadName();
    }

    public SAMRecord cloneRecord(){
        try {
            return (SAMRecord) record.clone();
        } catch (CloneNotSupportedException e) {
            e.printStackTrace();
            throw new RuntimeException("Error cloning record " + record.getReadName(), e);
        }
    }

    @Override
    public String getReadAlignmentId() {
        return readAlignmentId;
    }

    public void setReadAlignmentId(String readAlignmentId) {
        this.readAlignmentId = readAlignmentId;
    }

    public SAMFileHeader getHeader(){
        return record.getHeader();
    }

    @Override
    public String getReadName() {
        return record.getReadName();
    }

    public void setReadName(String readName) {
        record.setReadName(readName);
    }

    @Override
    public String getSequenceName() {
        return record.getReferenceName();
    }

    public void setSequenceName(String sequenceName) {
        record.setReferenceName(sequenceName);
    }

    public int getFlags() {
        return record.getFlags();
    }

    public void setFlags(int flags) {
        record.setFlags(flags);
    }

    public void setCigar(Cigar cigar) {
        record.setCigar(cigar);
    }

    @Override
    public int getStart() {
        return record.getAlignmentStart();
    }

    public void setStart(int start) {
        record.setAlignmentStart(start);
    }

    @Override
    public int getEnd() {
        return record.getAlignmentEnd();
    }

    @Override
    public int getLength() {
        return record.getReadLength();
    }

    @Override
    public byte[] getSequenceArray() {
        return record.getReadBases();
    }

    public void setSequenceArray(byte[] sequenceArray) {
        record.setReadBases(sequenceArray);
    }

    @Override
    public byte[] getQualitiesArray() {
        return record.getBaseQualities();
    }

    @Override
    public String getSequence() {
        return new String(getSequenceArray(), StandardCharsets.UTF_8);
    }

    @Override
    public String getQualities() {
        return Arrays.toString(getQualitiesArray());
    }

    @Override
    public int getMappingQuality() {
        return record.getMappingQuality();
    }

    @Override
    public Cigar getCigar() {
        return record.getCigar();
    }

    public List<CigarElement> getCigarElements(){
        return getCigar().getCigarElements();
    }

    public void setQualitiesArray(byte[] qualitiesArray) {
        record.setBaseQualities(qualitiesArray);
    }

    @Override
    public boolean isSupplementary() {
        return record.getSupplementaryAlignmentFlag();
    }

    public void setSupplementary(boolean supplementary) {
        record.setSupplementaryAlignmentFlag(supplementary);
    }

    @Override
    public boolean isSecondary() {
        return record.getNotPrimaryAlignmentFlag();
    }

    public void setSecondary(boolean secondary) {
        record.setNotPrimaryAlignmentFlag(secondary);
    }

    public List<SAMRecord.SAMTagAndValue> getTags(){
        return record.getAttributes();
    }
}
