package genomicelements;

import java.util.*;

public class ShortAnonymizedReadAlignment implements AnonymizedRead{

    private ShortReadAlignment readAlignment;
    private boolean isAnonymized;
    private byte[] anonymizedSequenceArray;
    private byte[] anonymizedQualitiesArray;
    private List<CalledVariation> SNVsToAnonymize;
    private List<CalledVariation> indelsToAnonymize;
    private List<CalledVariation> SVsToAnonymize;

    public ShortAnonymizedReadAlignment(ShortReadAlignment readAlignment) {
        this.isAnonymized = false;
        this.anonymizedSequenceArray = new byte[0];
        this.anonymizedQualitiesArray = new byte[0];
        SNVsToAnonymize = new ArrayList<>();
        indelsToAnonymize = new ArrayList<>();
        SVsToAnonymize = new ArrayList<>();
    }


    public void anonymizeVariants() {
        if (SNVsToAnonymize.isEmpty()){
            isAnonymized = true;
            return;
        }
        if(!SNVsToAnonymize.isEmpty()) {
            SNVsToAnonymize.sort(Comparator.comparing(var -> var.getInReadPosition(this)));
            for (CalledVariation var : SNVsToAnonymize) {
                int inReadPos = var.getInReadPosition(this);
                modifyBaseInRead(inReadPos, var.getRefAllele()[0]);
            }
        }
        if(!indelsToAnonymize.isEmpty()){
            indelsToAnonymize.sort(Comparator.comparing(var -> var.getInReadPosition(this)));
            int offset = 0;
            for (int i = 0; i < indelsToAnonymize.size(); i++) {
                CalledVariation var = indelsToAnonymize.get(i);
                int inReadPos = var.getInReadPosition(this);
                if(i==0) offset += modifyIndel(inReadPos, var);
                else offset += modifyIndel((inReadPos + offset)+1, var);
            }
        }
        isAnonymized = true;
    }

    public void modifyBaseInRead(int inReadPosition, byte asciiBase){
        int inArrayPosition = inReadPosition - 1;
        anonymizedSequenceArray[inArrayPosition] = asciiBase;
    }

    public void modifyBaseAndQualityInRead(int inReadPosition, byte asciiBase, byte asciiBaseQuality){
        modifyBaseInRead(inReadPosition, asciiBase);
        anonymizedQualitiesArray[inReadPosition] = asciiBaseQuality;
    }

//    private int modifyIndel(int inReadPosition, CalledVariation var) {
//        int addedOffset = 0;
//        int inArrayPosition = inReadPosition - 1;
//        int varLength = var.getLength();
//        byte [] newSequenceArray;
//        byte[] newQualitiesArray;
//        if (CalledVariation.VariantType.INS.equals(var.getVariantType())){
//            newSequenceArray = removeInsertion(sequenceArray, inArrayPosition, varLength, var);
//            newQualitiesArray = removeInsertion(qualitiesArray, inArrayPosition, varLength, var);
//            addedOffset = -(varLength);
//        }
//        else if (CalledVariation.VariantType.DEL.equals(var.getVariantType())){
//            newSequenceArray = removeDeletion(sequenceArray, inArrayPosition, varLength, var.getRefAllele());
//            byte[] avgQualities = new byte[var.getRefAllele().length];
//            byte avgQ = getAverageOfBytes(qualitiesArray);
//            Arrays.fill(avgQualities, avgQ);
//            newQualitiesArray = removeDeletion(qualitiesArray, inArrayPosition, varLength, avgQualities);
//            addedOffset = varLength;
//        }
//        else {
//            // Placeholder for other types of variants
//            newSequenceArray = sequenceArray;
//            newQualitiesArray = qualitiesArray;
//        }
//        return addedOffset;
//    }

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

    public String getReadName() {
        return readAlignment.getReadName();
    }

    public int getLength() {
        return readAlignment.getLength();
    }

    public int getPairIdx() {
        return readAlignment.getPairIdx();
    }

    public void setVariantsToAnonymize(List<CalledVariation> variants) {
        addAllVariantsToAnonymize(variants);
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
        boolean added = false;
        String varType = CalledVariation.VariantType.SNV.equals(variation.getVariantType()) ?
                CalledVariation.GENERIC_TYPE_SNV : CalledVariation.GENERIC_TYPE_INDEL;
        if (CalledVariation.GENERIC_TYPE_SNV.equals(varType)) {
            added = SNVsToAnonymize.add(variation);
        }
        if(CalledVariation.GENERIC_TYPE_INDEL.equals(varType)){
            added = indelsToAnonymize.add(variation);
        }
        //Add SVs also
        return added;
    }

    public String getReadId(){
        return readAlignment.getReadAlignmentId();
    }

    public String getSequenceName() {
        return readAlignment.getSequenceName();
    }

    public int getStart() {
        return readAlignment.getStart();
    }

    public int getEnd() {
        return readAlignment.getEnd();
    }

    public boolean isSupplementary() {
        return readAlignment.isSupplementary;
    }

    public boolean isAnonymized() {
        return isAnonymized;
    }

    /**
     * TO ERASE: DEPRECATED
     */
//    @Override
//    public FastqRecord getFastqRecord() {
//        if (isReverse){
//            SequenceUtil.reverseComplement(sequenceArray);
//            SequenceUtil.reverseQualities(qualitiesArray);
//        }
//        String name = this.readName;
//        //String readSequence = new String(sequenceArray, StandardCharsets.UTF_8);
//        //String qualitySequence = new String(qualitiesArray, StandardCharsets.UTF_8);
//        String comment = "";
//        return new FastqRecord(name, sequenceArray, comment, qualitiesArray);
//    }
}
