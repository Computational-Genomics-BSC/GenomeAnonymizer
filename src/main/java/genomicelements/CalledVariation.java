package genomicelements;

import htsjdk.variant.variantcontext.VariantContext;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import static utils.Operations.estimateEuclideanDistance;

public class CalledVariation {
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
    private SomaticVariationType somaticVariationType;
    private boolean hasDiffused;
    private boolean isLinkedToAnotherGermline;
    private Map<String, Integer> supportingReads;

    public CalledVariation(String seqName, int pos, int end, VariantType varType, int length, byte[] allele, byte[] refAllele) {
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

    /**
     * Constructor to create a CalledVarition object from the fields in a VariantContext from a VCF file
     * @param varContext
     */
    public CalledVariation(VariantContext varContext){
        this(varContext.getContig(), varContext.getStart(), varContext.getEnd(), VariantType.SNV, varContext.getLengthOnReference(),
                varContext.getAlternateAllele(0).getBases(), new byte[0]);
        this.setVariantType(getTypeFromVarContext(varContext));
        this.somaticVariationType = SomaticVariationType.UNCLASSIFIED;
        this.hasDiffused = false;
        this.isLinkedToAnotherGermline = false;
        this.supportingReads = new HashMap<>();
    }

    private VariantType getTypeFromVarContext(VariantContext context){
        VariantType varType = VariantType.SNV;
        if(context.isSNP()) varType = VariantType.SNV;
        if(context.isSimpleDeletion()) varType = VariantType.DEL;
        if(context.isSimpleInsertion()) varType = VariantType.INS;
        if(context.isSymbolicOrSV()) varType = VariantType.SV;
        return varType;
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

    public int getInReadPosition(AnonymizedRead anonRead){
        return supportingReads.get(anonRead.getUniqueReadName());
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

    public int getLength() {
        return length;
    }

    public byte[] getAllele() {
        return allele;
    }

    public byte[] getRefAllele() {
        return refAllele;
    }

    public Map<String, Integer> getSupportingReads() {
        return supportingReads;
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

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj instanceof CalledVariation) {
            CalledVariation var2 = (CalledVariation) obj;
            return this.seqName.equals(var2.seqName) &&
                    this.variantType.equals(var2.variantType) &&
                    this.pos == var2.pos &&
                    this.end == var2.end &&
                    //this.length == var2.length &&
                    Arrays.equals(this.allele, var2.allele);
        }
        return false;
    }

    @Override
    public String toString() {
        return "seq_name: " + seqName + " pos: " + pos + " end: " + end + " var_type: " + variantType +
                " length: " + length + " alt_allele: " + new String(allele, StandardCharsets.UTF_8) + " ref_allele: " + new String(refAllele, StandardCharsets.UTF_8) +
                " somatic_variation_type: " + somaticVariationType;
    }

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

