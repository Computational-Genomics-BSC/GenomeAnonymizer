package genomicelements;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.List;

/**
 * Base implementation of a generic read alignment, only to be used for inheritance
 * @author Nicolas Gaitan
 */
public class ReadAlignmentBaseImpl implements ReadAlignment{

    private String readAlignmentId;
    private String readName;
    private int flags;
    private String sequenceName;
    private int start;
    private int end;
    private int length;
    private byte[] sequenceArray;
    private byte[] qualitiesArray;
    private int mapQ;
    private Cigar cigar;
    boolean isSupplementary;
    boolean isSecondary;
    private SAMFileHeader header;


    public ReadAlignmentBaseImpl(String readName, String sequenceName, int flags, int start, int end,
                                 boolean isSupplementary, int mapQ, Cigar cigar) {
        this.readAlignmentId = generateReadAlignmentId(readName);
        this.readName = readName;
        this.flags = flags;
        this.sequenceName = sequenceName;
        this.start = start;
        this.end = end;
        this.length = end-start+1;
        this.isSupplementary = isSupplementary;
        this.mapQ = mapQ;
        this.cigar = cigar;
    }

    @Override
    public String getReadAlignmentId() {
        return readAlignmentId;
    }

    public void setReadAlignmentId(String readAlignmentId) {
        this.readAlignmentId = readAlignmentId;
    }

    public void setHeader(SAMFileHeader header){
        this.header = header;
    }

    public SAMFileHeader getHeader(){
        return header;
    }

    @Override
    public String getReadName() {
        return readName;
    }

    public void setReadName(String readName) {
        this.readName = readName;
    }

    @Override
    public String getSequenceName() {
        return sequenceName;
    }

    public void setSequenceName(String sequenceName) {
        this.sequenceName = sequenceName;
    }

    public int getFlags() {
        return flags;
    }

    public void setFlags(int flags) {
        this.flags = flags;
    }

    public int getMapQ() {
        return mapQ;
    }

    public void setMapQ(int mapQ) {
        this.mapQ = mapQ;
    }

    public void setCigar(Cigar cigar) {
        this.cigar = cigar;
    }

    @Override
    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    @Override
    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    @Override
    public int getLength() {
        return length;
    }

    public void setLength(int length) {
        this.length = length;
    }

    @Override
    public byte[] getSequenceArray() {
        return sequenceArray;
    }

    public void setSequenceArray(byte[] sequenceArray) {
        this.sequenceArray = sequenceArray;
    }

    @Override
    public byte[] getQualitiesArray() {
        return qualitiesArray;
    }

    @Override
    public String getSequence() {
        return new String(this.sequenceArray, StandardCharsets.UTF_8);
    }

    @Override
    public String getQualities() {
        return Arrays.toString(this.qualitiesArray);
    }

    @Override
    public int getMappingQuality() {
        return mapQ;
    }

    @Override
    public Cigar getCigar() {
        return cigar;
    }

    public List<CigarElement> getCigarElements(){
        return cigar.getCigarElements();
    }

    public void setQualitiesArray(byte[] qualitiesArray) {
        this.qualitiesArray = qualitiesArray;
    }

    @Override
    public boolean isSupplementary() {
        return isSupplementary;
    }

    public void setSupplementary(boolean supplementary) {
        isSupplementary = supplementary;
    }

    @Override
    public boolean isSecondary() {
        return isSecondary;
    }

    public void setSecondary(boolean secondary) {
        isSecondary = secondary;
    }

    // Read Id depends on the type of read that is aligned (e.g. Sequencing platform)
    public static String generateReadAlignmentId(String readName){
        return readName;
    }
}
