package GenomicElements;

import java.util.HashMap;
import java.util.Map;

import java.util.EnumSet;

import static Utils.Operations.estimateEuclideanDistance;

public class CalledVariation {
    private String seqName;
    private int pos;
    private int end;
    private VariantType variantType;
    private int length;
    private String allele;
    private String refAllele;
    private SomaticVariationType somaticVariationType;
    private boolean hasDiffused;
    private boolean isLinkedToAnotherGermline;
    private Map<String, Integer> supportingReads;

    public CalledVariation(String seqName, int pos, int end, VariantType varType, int length, String allele, String refAllele) {
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
        this.supportingReads = new HashMap<>();
    }

    /*
    public static CalledGenomicVariant fromVariantRecord(VariantRecord variantRecord) {
        // VariantRecord coordinates are 1-based, while CalledGenomicVariant coordinates are 0-based
        return new CalledGenomicVariant(variantRecord.getContig(), variantRecord.getPos() - 1, variantRecord.getEnd() - 1,
                variantRecord.getVariantType(), variantRecord.getLength(), variantRecord.getAlt(), variantRecord.getRef());
    }
    */

    public void addSupportingRead(String readId, int varReadPos) {
        supportingReads.put(readId, varReadPos);
    }

    public void setDiffusedStatus() {
        this.hasDiffused = true;
    }

    public void setLinkToAnotherGermline() {
        this.isLinkedToAnotherGermline = true;
    }

    public boolean isCandidateForDiffusion() {
        return !isLinkedToAnotherGermline;
    }

    public boolean hasDiffused() {
        return hasDiffused;
    }

    public double calculateDistanceToAnother(CalledVariation variant2) {
        return estimateEuclideanDistance(this.pos, this.end, this.length, variant2.pos, variant2.end, variant2.length);
    }

    public String getSeqName() {
        return seqName;
    }

    public int getPos() {
        return pos;
    }

    public int getEnd() {
        return end;
    }

    public VariantType getVariantType() {
        return variantType;
    }

    public void setVariantType(VariantType variantType) {
        this.variantType = variantType;
    }

    public int getLength() {
        return length;
    }

    public String getAllele() {
        return allele;
    }

    public void setAllele(String allele) {
        this.allele = allele;
    }

    public String getRefAllele() {
        return refAllele;
    }

    public void setRefAllele(String refAllele) {
        this.refAllele = refAllele;
    }

    public SomaticVariationType getSomaticVariationType() {
        return somaticVariationType;
    }

    public void setSomaticVariationType(SomaticVariationType somaticVariationType) {
        this.somaticVariationType = somaticVariationType;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (!(obj instanceof CalledVariation var2)) return false;
        return seqName.equals(var2.seqName) &&
                variantType.equals(var2.variantType) &&
                pos == var2.pos &&
                end == var2.end &&
                length == var2.length &&
                allele.equals(var2.allele);
    }

    @Override
    public String toString() {
        return "seq_name: " + seqName + " pos: " + pos + " end: " + end + " var_type: " + variantType +
                " length: " + length + " alt_allele: " + allele + " ref_allele: " + refAllele +
                " somatic_variation_type: " + somaticVariationType;
    }

    public enum VariantType {
        SNV(0),
        DEL(1),
        INS(2),
        LARGE_DEL(11),
        LARGE_INS(12),
        DUP(3),
        INV(4),
        CNV(5),
        TRA(6),
        SGL(7);

        private final int value;

        VariantType(int value) {
            this.value = value;
        }

        public int getValue() {
            return value;
        }
    }

    public enum SomaticVariationType {
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

