package genomicelements;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;

import java.nio.charset.StandardCharsets;
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
    private boolean isReverse;
    private boolean isSupplementaryOrSecondary;
    private boolean isAnonymized;
    private byte[] sequenceArray;
    private byte[] qualitiesArray;
    private Map<String, List<CalledVariation>> variantsToAnonymize;

    public ShortAnonymizedRead(String readName, String contig, int start, int end, boolean isReverse) {
        this.readName = readName;
        this.contig = contig;
        this.start = start;
        this.end = end;
        this.isReverse = isReverse;
        this.isAnonymized = false;
        this.sequenceArray = new byte[0];
        this.qualitiesArray = new byte[0];
        variantsToAnonymize = new HashMap<>();
    }

    public ShortAnonymizedRead(SAMRecord samRec) {
        this(samRec.getReadName(), samRec.getContig(), samRec.getStart(), samRec.getEnd(), samRec.getReadNegativeStrandFlag());
        int pairIdx = samRec.getFirstOfPairFlag() ? PAIR_1_IDX : PAIR_2_IDX;
        this.setPair(pairIdx);
        isSupplementaryOrSecondary = samRec.isSecondaryOrSupplementary();
        if (!isSupplementaryOrSecondary) setSequenceArray(samRec);
    }

//    public ShortAnonymizedRead fromSAMRecord(SAMRecord samRec) {
//        return new ShortAnonymizedRead(samRec, false);
//    }

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
            indelsToAnonymize.sort(Comparator.comparing(var -> var.getInReadPosition(this)));
            int offset = 0;
            for (CalledVariation var : indelsToAnonymize) {
                int inReadPos = var.getInReadPosition(this);
                offset += modifyIndel(inReadPos + offset, var);
            }
        }
        isAnonymized = true;
    }

    public void modifyBaseInRead(int inReadPosition, byte asciiBase){
        int inArrayPosition = inReadPosition - 1;
        // DEBUG
//        if (inArrayPosition==100 && sequenceArray.length==95){
//            System.out.println("# Read=" + this.readName + " seq=" + Arrays.toString(sequenceArray) + " inArrayPos=" + inArrayPosition + " base=" + asciiBase + " pair=" + pair);
//        }
        assert inArrayPosition < sequenceArray.length: "In read position is bigger than the length of the read sequence: readpos=" + inReadPosition + " seq_length=" + sequenceArray.length;
        // DEBUG
        sequenceArray[inArrayPosition] = asciiBase;
    }

    public void modifyBaseAndQualityInRead(int inReadPosition, byte asciiBase, byte asciiBaseQuality){
        modifyBaseInRead(inReadPosition, asciiBase);
        qualitiesArray[inReadPosition] = asciiBaseQuality;
    }

    private int modifyIndel(int inReadPosition, CalledVariation var) {
        int addedOffset = 0;
        int inArrayPosition = inReadPosition - 1;
        int varLength = var.getLength();
        byte [] newSequenceArray;
        byte[] newQualitiesArray;
        if (CalledVariation.VariantType.INS.equals(var.getVariantType())){
            // Deletes the insertion array
            newSequenceArray = removeInsertion(sequenceArray, inArrayPosition, varLength);
            newQualitiesArray = removeInsertion(qualitiesArray, inArrayPosition, varLength);
            addedOffset = -varLength;
        }
        else if (CalledVariation.VariantType.DEL.equals(var.getVariantType())){
            newSequenceArray = removeDeletion(sequenceArray, inArrayPosition, varLength, var.getRefAllele());
            byte[] avgQualities = new byte[var.getRefAllele().length];
            byte avgQ = getAverageOfBytes(qualitiesArray);
            Arrays.fill(avgQualities, avgQ);
            newQualitiesArray = removeDeletion(qualitiesArray, inArrayPosition, varLength, avgQualities);
            addedOffset = varLength;
        }
        else {
            // Placeholder for other types of variants
            newSequenceArray = sequenceArray;
            newQualitiesArray = qualitiesArray;
        }
        assert (newSequenceArray.length == newQualitiesArray.length): "Length of the modified qualities does not match the length of the modified sequence" +
                " seqLength=" + newSequenceArray.length + " qualLength=" + qualitiesArray.length;
        sequenceArray = newSequenceArray;
        qualitiesArray = newQualitiesArray;
        return addedOffset;
    }

    public byte getAverageOfBytes(byte[] array){
        int sum = 0;
        int denom = array.length;
        int answer;
        for (byte b : array) {
            sum += (int) b;
        }
        answer = sum / denom;
        return (byte) answer;
    }

    private byte[] removeInsertion(byte[] original, int inArrayPosition, int varLength){
        byte[] answer = new byte[original.length - varLength];
        //System.arraycopy(original, 0, answer, 0, inArrayPosition);
        System.arraycopy(original, 0, answer, 0, inArrayPosition + 1);
        System.arraycopy(original, inArrayPosition + varLength, answer, inArrayPosition + 1, answer.length - inArrayPosition - 1);
        //System.arraycopy(original, inArrayPosition + varLength + 1, answer, inArrayPosition + 1, answer.length - inArrayPosition - 1);
        return answer;
    }

    private byte[] removeDeletion(byte[] original, int inArrayPosition, int varLength, byte[] newContent){
        byte[] answer = new byte[original.length + varLength];
        System.arraycopy(original, 0, answer, 0, inArrayPosition+1);
        //System.arraycopy(original, 0, answer, 0, inArrayPosition + 1);
        assert (varLength == newContent.length): "Length of reference is not equal to varLength";
        //System.arraycopy(newContent, 1, answer, inArrayPosition + 1, varLength-1);
        System.arraycopy(newContent, 0, answer, inArrayPosition, varLength);
        //System.arraycopy(original, inArrayPosition + 1, answer, inArrayPosition + varLength-1, original.length - inArrayPosition - 1);
        System.arraycopy(original, inArrayPosition + 1, answer, inArrayPosition + varLength, original.length - inArrayPosition - 1);
        return answer;
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

    public boolean addAllVariantsToAnonymize(List<CalledVariation> variations){
        boolean added = false;
        for (CalledVariation var : variations){
            added = addVariantToAnonymize(var);
            if (!added) return false;
        }
        return added;
    }

    public boolean addVariantToAnonymize(CalledVariation variation){
        // Can change when dealing with SVs
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

    public boolean variantsToAnonymizeIsEmpty(){
        return variantsToAnonymize.isEmpty();
    }

    /**
     * Should only be called when the anonymized read is going to be written
     * @return FastqRecord representing the anonymized read
     */
    @Override
    public FastqRecord getFastqRecord() {
        if (isReverse){
            SequenceUtil.reverseComplement(sequenceArray);
            SequenceUtil.reverseQualities(qualitiesArray);
        }
        String name = this.readName;
        //String readSequence = new String(sequenceArray, StandardCharsets.UTF_8);
        //String qualitySequence = new String(qualitiesArray, StandardCharsets.UTF_8);
        String comment = "";
        return new FastqRecord(name, sequenceArray, comment, qualitiesArray);
    }

    public void updateIfPossible(SAMRecord samRecord){
        if(isSupplementaryOrSecondary && !samRecord.isSecondaryOrSupplementary()){
            setSequenceArray(samRecord);
        }
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

    public void setVariantsToAnonymize(Map<String, List<CalledVariation>> variantsToAnonymize) {
        this.variantsToAnonymize = variantsToAnonymize;
    }

    @Override
    public String toString(){
        return "name=" + this.readName + " seq=" + this.contig + " pos=" + this.start + " end=" +
                this.end + " length=" + this.length + " sequence=" +  new String(this.sequenceArray, StandardCharsets.UTF_8)
                + " qualities=" + Arrays.toString(this.qualitiesArray) + " pairIdx=" + this.pair +
                " isSupplementaryOrSecondary=" + this.isSupplementaryOrSecondary + " isAnonymized=" + this.isAnonymized + " isReverse=" + this.isReverse;
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
