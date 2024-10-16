package genomicelements;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.Locatable;

import java.util.*;

import static analysis.VariationClassifier.getShortReadPairName;
import static genomicelements.ShortAnonymizedReadPair.*;

public class ShortAnonymizedRead implements AnonymizedRead{

    private String readName;
    private String contig;
    private int start;
    private int end;
    private int length;
    private int pair;
    private boolean isSupplementaryOrSecondary;
    private boolean isAnonymized;
    private byte[] sequenceArray;
    private byte[] qualitiesArray;
    private Map<String, List<CalledVariation>> variantsToAnonymize;

    public ShortAnonymizedRead(String readName, String contig, int start, int end) {
        this.readName = readName;
        this.contig = contig;
        this.start = start;
        this.end = end;
        this.isAnonymized = false;
        this.sequenceArray = new byte[0];
        this.qualitiesArray = new byte[0];
        variantsToAnonymize = new HashMap<>();
    }

    public ShortAnonymizedRead fromSAMRecord(SAMRecord samRec, boolean updateSequenceIfPossible) {
        ShortAnonymizedRead instance = new ShortAnonymizedRead(samRec.getReadName(), samRec.getContig(), samRec.getStart(), samRec.getEnd());
        int pairIdx = samRec.getFirstOfPairFlag() ? PAIR_1_IDX : PAIR_2_IDX;
        instance.setPair(pairIdx);
        isSupplementaryOrSecondary = samRec.isSecondaryOrSupplementary();
        if (updateSequenceIfPossible && !isSupplementaryOrSecondary) setSequenceArray(samRec);
        return instance;
    }

    public ShortAnonymizedRead fromSAMRecord(SAMRecord samRec) {
        return fromSAMRecord(samRec, false);
    }

    @Override
    public void anonymizeVariantsInRead() {
        // DEBUG
        assert !isSupplementaryOrSecondary: "Trying to mask variants in AnonymizedRead " + getUniqueReadName() + " without processing the primary mapping";
        // DEBUG
        if (variantsToAnonymize.isEmpty()){
            isAnonymized = true;
            return;
        }
        // DEBUG
        assert sequenceArray.length>0: "Sequence array is not set, therefore anonymization should not be attempted";
        // DEBUG
        List<CalledVariation> snvsToAnonymize = variantsToAnonymize.getOrDefault(CalledVariation.GENERIC_TYPE_SNV, new ArrayList<>());
        List<CalledVariation> indelsToAnonymize = variantsToAnonymize.getOrDefault(CalledVariation.GENERIC_TYPE_INDEL, new ArrayList<>());;
        if(!snvsToAnonymize.isEmpty()) {
            snvsToAnonymize.sort(Comparator.comparing(var -> var.getInReadPosition(this)));
            for (CalledVariation var : snvsToAnonymize) {
                int inReadPos = var.getInReadPosition(this);
                modifyBaseInRead(inReadPos, var.getRefAllele()[0]);
            }
        }
        if(!indelsToAnonymize.isEmpty()){
            int offset = 0;
            for (CalledVariation var : snvsToAnonymize) {
                int inReadPos = var.getInReadPosition(this);
                offset += modifyIndel(inReadPos + offset, var);
            }
        }
        isAnonymized = true;
    }

    public void modifyBaseInRead(int inReadPosition, byte asciiBase){
        int inArrayPosition = inReadPosition - 1;
        // DEBUG
        assert inArrayPosition <= sequenceArray.length: "In read position is bigger than the length of the read sequence: readpos=" + inReadPosition + " seq_length=" + sequenceArray.length;
        // DEBUG
        sequenceArray[inArrayPosition] = asciiBase;
    }

    private int modifyIndel(int inReadPosition, CalledVariation var) {
        int addedOffset = 0;

        return addedOffset;
    }

    public void modifyBaseAndQualityInRead(int inReadPosition, byte asciiBase, byte asciiBaseQuality){
        modifyBaseInRead(inReadPosition, asciiBase);
        qualitiesArray[inReadPosition] = asciiBaseQuality;
    }

    public boolean isPair1(){
        return pair == PAIR_1_IDX;
    }

    public boolean isPair2(){
        return pair == PAIR_2_IDX;
    }

    public String getReadName() {
        return readName;
    }

    public String getUniqueReadName() {
        return getShortReadPairName(this.readName, this.getPairIdx());
    }

    public int getLength() {
        return length;
    }

    @Override
    public byte[] getSequenceArray() {
        return sequenceArray;
    }

    public byte[] getQualitiesArray() {
        return qualitiesArray;
    }

    public Map<String, List<CalledVariation>> getVariantsToAnonymize() {
        return variantsToAnonymize;
    }

    public int getPairIdx() {
        return pair;
    }

    public boolean addVariantToAnonymize(CalledVariation variation){
        String varType = CalledVariation.VariantType.SNV.equals(variation.getVariantType()) ?
                CalledVariation.GENERIC_TYPE_SNV : CalledVariation.GENERIC_TYPE_INDEL;
        List<CalledVariation> variationPerType = variantsToAnonymize.computeIfAbsent(varType,
                v -> new ArrayList<>());
        return variationPerType.add(variation);
    }

    @Override
    public String getContig() {
        return contig;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getEnd() {
        return end;
    }

    public boolean isSupplementaryOrSecondary() {
        return isSupplementaryOrSecondary;
    }

    @Override
    public boolean isAnonymized() {
        return isAnonymized;
    }

    @Override
    public FastqRecord getFastqRecord() {

        return null;
    }

    /**
     * <pre> Sequence array must come from the primary mapping in any case </pre>
     * <pos> This anonymized read is no longer considered supplementary nor secondary </pos>
     * @param sequenceArray sequences as array of bytes
     * @param qualitiesArray base qualities as array of bytes
     */
    public void setSequenceArray(byte[] sequenceArray, byte[] qualitiesArray) {
        this.sequenceArray = sequenceArray;
        this.qualitiesArray = qualitiesArray;
        this.isSupplementaryOrSecondary = false;
    }

    /**
     * <pre> samRecord must come from the primary alignment </pre>
     * @param samRecord
     */
    public void setSequenceArray(SAMRecord samRecord){
        setSequenceArray(samRecord.getReadBases(), samRecord.getBaseQualities());
    }

    public void setPair(int pair){
        this.pair = pair;
    }

    // TODO: Fix these methods

    @Override
    public int getLengthOnReference() {
        return 0;
    }

    @Override
    public boolean overlaps(Locatable other) {
        return false;
    }

    @Override
    public boolean withinDistanceOf(Locatable other, int distance) {
        return false;
    }

    @Override
    public boolean contains(Locatable other) {
        return false;
    }

    @Override
    public boolean contigsMatch(Locatable other) {
        return false;
    }
}
