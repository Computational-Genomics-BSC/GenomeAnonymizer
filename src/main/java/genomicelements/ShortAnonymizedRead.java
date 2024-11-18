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
                //DEBUG
                if(getReadName().equals("6a129c2336afd2745c10a3b7ca407903") || getReadName().equals("f1f575b69219b650a8e9900e12c9acf8")){
                    System.out.println("Read before SNV masking: " + this);
                    //if(samRecord.isSecondaryOrSupplementary() && calledVar.getSomaticVariationType().equals(SomaticVariationType.TUMORAL_NORMAL_VARIANT)){
                    System.out.println("$Found var: " + var + " in suppl.=" + " in_read_pos=" + inReadPos);
                    //}
                }
                //DEBUG
                modifyBaseInRead(inReadPos, var.getRefAllele()[0]);
                //DEBUG
                if(getReadName().equals("6a129c2336afd2745c10a3b7ca407903") || getReadName().equals("f1f575b69219b650a8e9900e12c9acf8")){
                    System.out.println("Read after SNV masking: " + this);
                    //if(samRecord.isSecondaryOrSupplementary() && calledVar.getSomaticVariationType().equals(SomaticVariationType.TUMORAL_NORMAL_VARIANT)){
                    //System.out.println("$Found var: " + var + " in suppl.=" + " in_read_pos=" + inReadPos);
                    //}
                }
                //DEBUG
            }
        }
        if(!indelsToAnonymize.isEmpty()){
            indelsToAnonymize.sort(Comparator.comparing(var -> var.getInReadPosition(this)));
            int offset = 0;
            for (int i = 0; i < indelsToAnonymize.size(); i++) {
                CalledVariation var = indelsToAnonymize.get(i);
                int inReadPos = var.getInReadPosition(this);
                //offset += modifyIndel(inReadPos + offset + 1, var);
                if(i==0) offset += modifyIndel(inReadPos, var);
                else offset += modifyIndel((inReadPos + offset)+1, var);
                //else offset += modifyIndel((inReadPos + offset)-1, var);
                //if(i==0) offset--;
                //DEBUG
                if(getReadName().equals("6a129c2336afd2745c10a3b7ca407903") || getReadName().equals("f1f575b69219b650a8e9900e12c9acf8")){
                    System.out.println("Read after INDEL masking: " + this);
                    System.out.println("$Found var: " + var + " in suppl.=" + " org_in_read_pos=" + inReadPos + " new_offset=" + offset);
                }
                //DEBUG
            }
        }
        isAnonymized = true;
    }

    public void modifyBaseInRead(int inReadPosition, byte asciiBase){
        int inArrayPosition = inReadPosition - 1;
        assert inArrayPosition < sequenceArray.length: "In read position is bigger than the length of the read sequence: readpos=" + inReadPosition + " seq_length=" + sequenceArray.length;
        sequenceArray[inArrayPosition] = asciiBase;
    }

    public void modifyBaseAndQualityInRead(int inReadPosition, byte asciiBase, byte asciiBaseQuality){
        modifyBaseInRead(inReadPosition, asciiBase);
        qualitiesArray[inReadPosition] = asciiBaseQuality;
    }

    private int modifyIndel(int inReadPosition, CalledVariation var) {
        int addedOffset = 0;
        int inArrayPosition = inReadPosition - 1;
        //DEBUG
//        if(readName.equals("e064a13f1bd6b357f76262ae9b7bf68f") && pair == PAIR_2_IDX){
//            System.out.println("original_pos=" + inArrayPosition);
//        }
        //DEBUG
        int varLength = var.getLength();
        byte [] newSequenceArray;
        byte[] newQualitiesArray;
        if (CalledVariation.VariantType.INS.equals(var.getVariantType())){
            // Deletes the insertion array
            //DEBUG
            //if(readName.equals("7990d52d9859724cb58d1188107eea6a") && pair == PAIR_2_IDX){
            //    System.out.println("before remove call: inArrayPosition=" + inArrayPosition + " -(varLength-1)=" + -(varLength-1)
            //    );//+" (varLength-1)=" + (varLength-1));
            //}
            //DEBUG
            newSequenceArray = removeInsertion(sequenceArray, inArrayPosition, varLength, var);
            newQualitiesArray = removeInsertion(qualitiesArray, inArrayPosition, varLength, var);
            //addedOffset = -(varLength-1);
            addedOffset = -(varLength);
            //DEBUG
            // Deletes the insertion array
            //DEBUG
//            if(readName.equals("7990d52d9859724cb58d1188107eea6a") && pair == PAIR_2_IDX){
//                System.out.println("after remove call: inArrayPosition=" + inArrayPosition + " -(varLength)=" + -(varLength-1) +
//                        " (varLength-1)+1=" + (varLength-1)+1);
//            }
            //DEBUG
        }
        else if (CalledVariation.VariantType.DEL.equals(var.getVariantType())){
            newSequenceArray = removeDeletion(sequenceArray, inArrayPosition, varLength, var.getRefAllele());
            byte[] avgQualities = new byte[var.getRefAllele().length];
            byte avgQ = getAverageOfBytes(qualitiesArray);
            Arrays.fill(avgQualities, avgQ);
            newQualitiesArray = removeDeletion(qualitiesArray, inArrayPosition, varLength, avgQualities);
            addedOffset = varLength;
            //addedOffset = (varLength+1);
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

    private byte[] removeInsertion(byte[] original, int inArrayPosition, int varLength, CalledVariation debugParam){
        byte[] answer = new byte[original.length - varLength];
        //System.arraycopy(original, 0, answer, 0, inArrayPosition);
        //DEBUG
        //if(inArrayPosition + 1 >= answer.length){
        //    System.out.println(debugParam + " inArrayPosition=" + inArrayPosition + " read=" + readName + " pair=" + pair);
        //    return original;
        //}
        //DEBUG
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
