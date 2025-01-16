package genomicelements;

import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;

import java.nio.charset.StandardCharsets;
import java.util.*;

/**
 * Class that generates anonymized versions of short read alignments, ready for writing
 * @author Nicolas Gaitan
 * @author Rodrigo Martin
 */
public class ShortAnonymizedReadAlignment implements AnonymizedRead{

    private final ShortReadAlignment readAlignment;
    private byte[] referenceContigSequence;
    private boolean isAnonymized;
    private byte[] anonymizedSequenceArray;
    private byte[] anonymizedQualitiesArray;
    List<CigarElement> anonymizedCigarElements;
    private Cigar anonymizedCigar;
    private List<CalledVariation> SNVsToAnonymize;
    private List<CalledVariation> indelsToAnonymize;
    private List<CalledVariation> SVsToAnonymize;

    public ShortAnonymizedReadAlignment(ShortReadAlignment readAlignment) {
        this.readAlignment = readAlignment;
        this.isAnonymized = false;
        this.anonymizedSequenceArray = new byte[0];
        this.anonymizedQualitiesArray = new byte[0];
        this.anonymizedCigarElements = new ArrayList<>();
        this.SNVsToAnonymize = new ArrayList<>();
        this.indelsToAnonymize = new ArrayList<>();
        this.SVsToAnonymize = new ArrayList<>();
    }

    public void setReferenceContigSequence(byte[] referenceContigSequence) {
        // 1-based memoized reference sequence, corresponding to the contig to which this read is mapped
        this.referenceContigSequence = referenceContigSequence;
    }

    public void anonymizeVariants() throws IllegalStateException{
        if(referenceContigSequence ==null){
            throw new IllegalStateException("The ShortAnonymizedReadAlignment.anonymizeVariants was called" +
                    "without setting the contigReferenceSequence first. setContigReferenceSequence, should always" +
                    " be alled after the constructor.");
        }
        int originalSeqLength = getOriginalSequenceArray().length;
        int expectedSize = estimateNewReadSize(originalSeqLength);
        anonymizedSequenceArray = new byte[expectedSize];
        anonymizedQualitiesArray = new byte[expectedSize];
        byte avgQual = getAverageOfBytes(getOriginalQualitiesArray());
        List<CigarElement> originalCigarElements = readAlignment.getCigarElements();
        //Invariant as the anonymized read will map to the exact start coordinate as the original
        final int alnStart = getStart();
        int alnOrgEnd = getEnd();
        //1 - based coordinate
        int currentRefAlnPos = alnStart-1;
        //Should be the same if we fill or remove bases from the end of the read
        //int alnNewEnd = 0;
        //Holds the position over the original read
        int i = 0;
        //Holds the position over the anonymized read
        int j = 0;
        //Holds the current index of the CIGAR alignment elements
        int c = 0;
        //SNV Operations (value: byte as new base) to perform using the original read index
        byte[] SNVops = processSNVoperations(originalSeqLength);
        //Indel Operations (value: number of bases to remove or add) to perform using the CIGAR index
        int[] indelOps = processIndelOps(originalCigarElements.size());
        while(c < originalCigarElements.size()){
            CigarElement cigarElem = originalCigarElements.get(c);
            CigarOperator currentCigarOp = cigarElem.getOperator();
            int opLength = cigarElem.getLength();
            int indelOp = indelOps[c];
            if(indelOp > 0){
                makeAdditiveChange(j, currentRefAlnPos, avgQual, indelOp);
                //i++;
                j += indelOp;
                currentRefAlnPos += indelOp;
            }
            else if(indelOp < 0){
                //makeSubstractiveChange();
                i += Math.abs(indelOp);
                //j++;
                //currentRefAlnPos++;
            }
            else{
                if(CigarOperator.M.equals(currentCigarOp)){
                    for(int x = 0; x < opLength; x++){
                        //DEBUG
//                        System.out.println("x="+ x);
//                        System.out.println("i="+ i);
//                        if(i==SNVops.length) {
//                            System.out.println(originalCigarElements + " "+ currentCigarOp);
//                        }
                        //DEBUG
                        byte snvOp = SNVops[i];
                        if(snvOp > 0){
                            anonymizedSequenceArray[j] = snvOp;
                        }
                        else{
                            anonymizedSequenceArray[j] = getOriginalSequenceArray()[i];
                        }
                        anonymizedQualitiesArray[j] = getOriginalQualitiesArray()[i];
                        i++;
                        j++;
                        currentRefAlnPos++;
                    }
                    //anonymizedCigarElements.add(new CigarElement(opLength, currentCigarOp));
                }
                else if(currentCigarOp.consumesReadBases()){
                    for(int x = 0; x < opLength; x++){
                        anonymizedSequenceArray[j] = getOriginalSequenceArray()[i];
                        anonymizedQualitiesArray[j] = getOriginalQualitiesArray()[i];
                        i++;
                        j++;
                    }
                    //anonymizedCigarElements.add(new CigarElement(opLength, currentCigarOp));
                    //i++;
                    //j++;
                    //currentRefAlnPos += opLength;
                }
                else{
                    currentRefAlnPos += opLength;
                    //anonymizedCigarElements.add(new CigarElement(opLength, currentCigarOp));
                }
                anonymizedCigarElements.add(new CigarElement(opLength, currentCigarOp));
            }
            c++;
        }
        generateDefinitiveCigar();
        isAnonymized = true;
    }

    private void makeAdditiveChange(int initPos, int refInitPos, byte avgQual, int length) {
        int j = initPos;
        int r = refInitPos;
        for(int l = 0; l < length; l++){
            anonymizedSequenceArray[j] = referenceContigSequence[r];
            anonymizedQualitiesArray[j] = avgQual;
            j++;
            r++;
        }
        anonymizedCigarElements.add(new CigarElement(length, CigarOperator.M));
    }

    private int estimateNewReadSize(int originalSeqLength) {
        int newSize = originalSeqLength;
        for(CalledVariation indel : indelsToAnonymize) {
            CalledVariation.VariantType variantType = indel.getVariantType();
            if (CalledVariation.VariantType.DEL.equals(variantType)) {
                newSize += indel.getLength();
            }
            if (CalledVariation.VariantType.INS.equals(variantType)) {
                newSize -= indel.getLength();
            }
            //TODO: Account for SVs (SoftClips at first)
        }
        return newSize;
    }

    private byte[] processSNVoperations(int originalSeqLength) {
        byte[] snvOps = new byte[originalSeqLength];
        for (CalledVariation snv : SNVsToAnonymize){
            int snvOpPosition = snv.getInReadPosition(this)-1;
            //Change to retrieve from memoized ref genome
            byte op = snv.getRefAllele()[0];
            snvOps[snvOpPosition] = op;
        }
        return snvOps;
    }

    private int[] processIndelOps(int originalCigarLength) {
        int[] indelOps = new int[originalCigarLength];
        for (CalledVariation indel : indelsToAnonymize){
            //TODO: Account for SVs (SoftClips at first)
            int indelOpPos = indel.getInReadPosition(this);
            int op = CalledVariation.VariantType.INS.equals(indel.getVariantType()) ?
                    -(indel.getLength()) : indel.getLength();
            indelOps[indelOpPos] = op;
        }
        return indelOps;
    }

    private void generateDefinitiveCigar() {
        List<CigarElement> fixedCigarElements = new ArrayList<>();
        boolean previousMerged = false;
        CigarElement currentElement = null;
        CigarElement nextElement = anonymizedCigarElements.get(0);
        CigarOperator currentOp;
        //If there is only 1 element, nextElement will hold it and the cycle does not happen
        for(int i = 0; i < anonymizedCigarElements.size()-1; i++){
            if(!previousMerged) currentElement = anonymizedCigarElements.get(i);
            else previousMerged = false;
            currentOp = currentElement.getOperator();
            nextElement = anonymizedCigarElements.get(i+1);
            CigarOperator nextOp = nextElement.getOperator();
            if(nextOp.equals(currentOp)){
                currentElement = new CigarElement(currentElement.getLength() + nextElement.getLength(),
                        currentOp);
                previousMerged = true;
            }
            else{
                fixedCigarElements.add(currentElement);
            }
        }
        if (!previousMerged) fixedCigarElements.add(nextElement);
        else fixedCigarElements.add(currentElement);
        anonymizedCigarElements = fixedCigarElements;
        anonymizedCigar = new Cigar(anonymizedCigarElements);
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

    public SAMRecord getAnonymizedSamRecord(){
        SAMRecord answer = readAlignment.cloneRecord();
        answer.setReadBases(anonymizedSequenceArray);
        answer.setBaseQualities(anonymizedQualitiesArray);
        answer.setCigar(anonymizedCigar);
        // getStart and getEnd are 1-based, adjust accordingly, currently referenceSequence is 0-based
        byte[] refSequenceAln = Arrays.copyOfRange(referenceContigSequence, answer.getStart()-1, answer.getEnd());
        SequenceUtil.calculateMdAndNmTags(answer, refSequenceAln, true, true);
        return answer;
    }

    public String getReadAlignmentId() {
        return readAlignment.getReadAlignmentId();
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

    public byte[] getOriginalSequenceArray(){
        return readAlignment.getSequenceArray();
    }

    public byte[] getOriginalQualitiesArray(){
        return readAlignment.getQualitiesArray();
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

    public int getMappingQuality(){
        return readAlignment.getMappingQuality();
    }

    public boolean isSupplementary() {
        return readAlignment.isSupplementary();
    }

    public boolean isAnonymized() {
        return isAnonymized;
    }

    @Override
    public String toString(){
        StringBuilder builder = new StringBuilder();
        builder.append("ReadName=").append(getReadName());
        builder.append(" ReadID=").append(getReadId());
        builder.append(" Start=").append(getStart());
        builder.append(" End=").append(getEnd());
        builder.append(" OrgSeq=").append(new String(getOriginalSequenceArray(), StandardCharsets.UTF_8));
        builder.append(" AnonSeq=").append(new String(anonymizedSequenceArray, StandardCharsets.UTF_8));
        builder.append(" OrgQual=").append(Arrays.toString(getOriginalQualitiesArray()));
        builder.append(" AnonQual=").append(Arrays.toString(anonymizedQualitiesArray));
        builder.append(" CIGAR=").append(anonymizedCigar.toString());
        builder.append(" TAGS=").append(readAlignment.getTags());
        return builder.toString();
    }
}
