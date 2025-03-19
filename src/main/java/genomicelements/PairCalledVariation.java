package genomicelements;

import java.nio.charset.StandardCharsets;
import java.util.*;

import static utils.Operations.compare;
import static utils.Operations.computeThreeDimEuclideanDistance;

/**
 * Class that represents a genomic variation with supporting evidence,
 * with the caveat that it is does not always represent a real genomic variant
 * @author Nicolas Gaitan
 */
public class PairCalledVariation implements GenomicRegion {
    public static final String GENERIC_TYPE_SNV = "SNV";
    public static final String GENERIC_TYPE_INDEL = "INDEL";
    public static final String GENERIC_TYPE_SV = "SV";

    private String seqName;
    private int pos;
    private int end;
    private VariantType variantType;
    private int length;
    private byte[] allele;
    private byte[] refAllele;
    private BreakendSVRecord breakendRecord;
    private ShorthandSVRecord shortHandRecord;
    private SomaticVariationType somaticVariationType;
    private boolean hasDiffused;
    private boolean isLinkedToAnotherGermline;
    private Map<String, Integer> supportingReadPositions;

    public PairCalledVariation(String seqName, int pos, int end, VariantType varType, int length, byte[] allele, byte[] refAllele) {
        this.seqName = seqName;
        this.pos = pos;
        this.end = end;
        this.variantType = varType;
        this.length = length;
        this.allele = allele;
        this.refAllele = refAllele;
        this.somaticVariationType = SomaticVariationType.UNCLASSIFIED;
        this.hasDiffused = false;
        this.isLinkedToAnotherGermline = false;
        //Holds the reads that supported this call as keys, and their 1-based position in-read for SNVs, or 0 based index in-CIGAR for indels and SVs
        this.supportingReadPositions = new HashMap<>();
        // OPTIONAL: Breakend and Shorthand record representations for SVs
        this.breakendRecord = null;
        this.shortHandRecord = null;
    }

    /*
    public static CalledGenomicVariant fromVariantRecord(VariantRecord variantRecord) {
        // VariantRecord coordinates are 1-based, while CalledGenomicVariant coordinates are 0-based
        return new CalledGenomicVariant(variantRecord.getContig(), variantRecord.getPos() - 1, variantRecord.getEnd() - 1,
                variantRecord.getVariantType(), variantRecord.getLength(), variantRecord.getAlt(), variantRecord.getRef());
    }
    */

    public void addSupportingRead(String readId, int varReadPos) {
        supportingReadPositions.put(readId, varReadPos);
    }

    public void setDiffusedStatus(boolean diffusedStatus) {
        this.hasDiffused = diffusedStatus;
    }

    public void setLinkToAnotherGermline(boolean isLinked) {
        this.isLinkedToAnotherGermline = isLinked;
    }

    public boolean isCandidateForDiffusion() {
        return !isLinkedToAnotherGermline;
    }

    public boolean hasDiffused() {
        return hasDiffused;
    }

    public double calculateDistanceToAnother(PairCalledVariation variant2) {
        return computeThreeDimEuclideanDistance(this.pos, this.end, this.length, variant2.pos, variant2.end, variant2.length);
    }

    public int getInReadPosition(AnonymizedRead anonRead){
        return getInReadPosition(anonRead.getReadAlignmentId());
    }

    public int getInReadPosition(String anonReadId){
        return supportingReadPositions.get(anonReadId);
    }

    @Override
    public String getSequenceName() {
        return seqName;
    }

    @Override
    public int getSequenceIdx() {
        return 0;
    }

    public int getStart() {
        return pos;
    }

    public int getEnd() {
        return end;
    }

    @Override
    public void setSequenceIdx(int sequenceIdx) {

    }

    public VariantType getVariantType() {
        return variantType;
    }

    public int getLength() {
        return length;
    }

    public byte[] getAllele() {
        return allele;
    }

    public byte[] getRefAllele() {
        return refAllele;
    }

    public Map<String, Integer> getSupportingReadPositions() {
        return supportingReadPositions;
    }

    public void setEnd(int end){
        this.end = end;
    }

    public void setLength(int length) {
        this.length = length;
    }

    public void setVariantType(VariantType variantType) {
        this.variantType = variantType;
    }

    public void setAllele(byte[] allele) {
        this.allele = allele;
    }

    public void setRefAllele(byte[] refAllele) {
        this.refAllele = refAllele;
    }

    public SomaticVariationType getSomaticVariationType() {
        return somaticVariationType;
    }

    public void setSomaticVariationType(SomaticVariationType somaticVariationType) {
        this.somaticVariationType = somaticVariationType;
    }

    public BreakendSVRecord getBreakendRecord() {
        return breakendRecord;
    }

    public void setBreakendRecord(BreakendSVRecord breakendRecord) {
        this.breakendRecord = breakendRecord;
    }

    public ShorthandSVRecord getShortHandRecord() {
        return shortHandRecord;
    }

    public void setShortHandRecord(ShorthandSVRecord shortHandRecord) {
        this.shortHandRecord = shortHandRecord;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj instanceof PairCalledVariation) {
            boolean answer;
            PairCalledVariation var2 = (PairCalledVariation) obj;
            answer =  this.seqName.equals(var2.seqName) &&
                    this.variantType.equals(var2.variantType) &&
                    this.pos == var2.pos &&
                    this.end == var2.end &&
                    this.length == var2.length &&
                    Arrays.equals(this.allele, var2.allele);
            return answer;
        }
        return false;
    }

    @Override
    public String toString() {
        return "seq_name: " + seqName + " pos: " + pos + " end: " + end + " var_type: " + variantType +
                " length: " + length + " alt_allele: " + new String(allele, StandardCharsets.UTF_8) + " ref_allele: " + new String(refAllele, StandardCharsets.UTF_8) +
                " somatic_variation_type: " + somaticVariationType;
    }

    @Override
    public int compareTo(GenomicRegion genomicRegion) {
        return compare(this, genomicRegion);
    }

    public record BreakendSVRecord(String prefix, String bracket, String contig, int pos, String suffix) {}

    public record ShorthandSVRecord(String type, String[] extraInfo){}

    public enum VariantType {
        SNV(0, "SNV"),
        DEL(1, "DEL"),
        INS(2, "INS"),
        LARGE_DEL(11, "LARGE_DEL"),
        LARGE_INS(12, "LARGE_INS"),
        DUP(3, "DUP"),
        INV(4, "INV"),
        CNV(5, "CNV"),
        TRA(6, "TRA"),
        SGL(7, "SGL"),
        //Generic type used as placeholder for now
        SV(13, "SV");;



        private final int value;
        private final String name;

        VariantType(int value, String name) {
            this.value = value;
            this.name = name;
        }

        public int getValue() {
            return value;
        }

        public String getName() {
            return name;
        }
    }

    public enum SomaticVariationType {
        NOT_SOMATIC(-1),
        UNCLASSIFIED(0),
        NORMAL_SINGLE_READ_VARIANT(1),
        TUMORAL_SINGLE_READ_VARIANT(2),
        NORMAL_ONLY_VARIANT(3),
        TUMORAL_ONLY_VARIANT(4),
        TUMORAL_NORMAL_VARIANT(5);

        private final int value;

        SomaticVariationType(int value) {
            this.value = value;
        }

        public int getValue() {
            return value;
        }
    }
}

