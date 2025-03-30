package genomicelements;

import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;

import java.nio.charset.StandardCharsets;
import java.util.*;

import static utils.Operations.compare;

/**
 * Class that generates anonymized versions of short read alignments, ready for writing
 * @author Nicolas Gaitan
 * @author Rodrigo Martin
 */
public class ShortAnonymizedReadAlignment implements AnonymizedRead, GenomicRegion{

    public static final int PAIR_1_IDX = 0;
    public static final int PAIR_2_IDX = 1;
    public static final String DEFAULT_ID_NAME_SEPARATOR = ";";

    private final SAMRecord readAlignment;
    private String readAlnId;
    private byte[] referenceContigSequence;
    private int alnStart;
    private boolean isAnonymized;
    private byte[] anonymizedSequenceArray;
    private byte[] anonymizedQualitiesArray;
    List<CigarElement> anonymizedCigarElements;
    private Cigar anonymizedCigar;
    private List<PairCalledVariation> SNVsimpleSignals;
    private List<PairCalledVariation> indelSimpleSignals;
    private List<Signal> complexSignals;
    private boolean isNormalDataset = true;

    public ShortAnonymizedReadAlignment(SAMRecord readAlignment) {
        this.readAlignment = readAlignment;
        this.readAlnId = generateReadAlnId(readAlignment);
        this.alnStart = getStart();
        this.isAnonymized = false;
        this.anonymizedSequenceArray = readAlignment.getReadBases();
        this.anonymizedQualitiesArray = readAlignment.getBaseQualities();
        this.anonymizedCigar = readAlignment.getCigar();
        this.anonymizedCigarElements = new ArrayList<>();
        this.SNVsimpleSignals = new ArrayList<>();
        this.indelSimpleSignals = new ArrayList<>();
        this.complexSignals = new ArrayList<>();
    }

    public ShortAnonymizedReadAlignment(SAMRecord readAlignment, boolean isNormalDataset) {
        this(readAlignment);
        this.isNormalDataset = isNormalDataset;
    }

    public SAMRecord getAnonymizedSamRecord() {
        if (!isAnonymized()) {
            anonymizeVariants();
        }
        SAMRecord answer = this.cloneRecord();
        answer.setAlignmentStart(alnStart);
        answer.setReadBases(anonymizedSequenceArray);
        answer.setBaseQualities(anonymizedQualitiesArray);
        answer.setCigar(anonymizedCigar);
        SequenceUtil.calculateMdAndNmTags(answer, referenceContigSequence, true, true);
        return answer;
    }

    public void setReferenceContigSequence(byte[] referenceContigSequence) {
        // 0-based memoized reference sequence, corresponding to the contig to which this read is mapped
        this.referenceContigSequence = referenceContigSequence;
    }

    private void anonymizeVariants() throws IllegalStateException{
        if(referenceContigSequence == null){
            throw new IllegalStateException("The ShortAnonymizedReadAlignment.anonymizeVariants was called" +
                    "without setting the contigReferenceSequence first. setContigReferenceSequence, should always" +
                    " be alled after the constructor.");
        }
        int originalSeqLength = getOriginalSequenceArray().length;
        int expectedSize = estimateNewReadSize(originalSeqLength);
        anonymizedSequenceArray = new byte[expectedSize];
        anonymizedQualitiesArray = new byte[expectedSize];
        byte avgQual = getAverageOfBytes(getOriginalQualitiesArray());
        List<CigarElement> originalCigarElements = readAlignment.getCigar().getCigarElements();
        //Except a starting PG SoftClip, alignment start is invariant, as the anonymized read will map to the exact start coordinate as the original
        //final int alnStart = getStart();
        int alnOrgEnd = getEnd();
        //0 based coordinate
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
        int[] indelOps = processCigarOperations(originalCigarElements.size());
        while(c < originalCigarElements.size()){
            CigarElement cigarElem = originalCigarElements.get(c);
            CigarOperator currentCigarOp = cigarElem.getOperator();
            int opLength = cigarElem.getLength();
            int indelOp = indelOps[c];
            if(indelOp > 0){
                //For softclips, if it begins the read, substract exactly the length from currentRefAlnPos, if it ends add also exactly the length
                if(CigarOperator.S == currentCigarOp){
                    currentRefAlnPos -= indelOp;
                    //Update alignment start in 1-based coordinates
                    alnStart = currentRefAlnPos+1;
                    i += indelOp;
                }
                makeAdditiveChange(j, currentRefAlnPos, avgQual, indelOp);
                j += indelOp;
                currentRefAlnPos += indelOp;
            }
            else if(indelOp < 0){
                i += Math.abs(indelOp);
            }
            else{
                if(CigarOperator.M.equals(currentCigarOp)){
                    for(int x = 0; x < opLength; x++){
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
                }
                else if(currentCigarOp.consumesReadBases()){
                    for(int x = 0; x < opLength; x++){
                        anonymizedSequenceArray[j] = getOriginalSequenceArray()[i];
                        anonymizedQualitiesArray[j] = getOriginalQualitiesArray()[i];
                        i++;
                        j++;
                    }
                }
                else{
                    currentRefAlnPos += opLength;
                }
                anonymizedCigarElements.add(new CigarElement(opLength, currentCigarOp));
            }
            c++;
        }
        // Add reference bases to fill read to its original length (or more) for base-removing operations, exclude reads that would fall out of reference bounds
        if(j < expectedSize){
            makeAdditiveChange(j, currentRefAlnPos, avgQual, expectedSize-j);
        }
        // Cut the read if the new size is larger than the original
        int cutLength;
        if(anonymizedSequenceArray.length > originalSeqLength){
            cutLength = anonymizedSequenceArray.length - originalSeqLength;
            anonymizedSequenceArray = Arrays.copyOf(anonymizedSequenceArray, originalSeqLength);
            anonymizedQualitiesArray = Arrays.copyOf(anonymizedQualitiesArray, originalSeqLength);
            // Adjust CIGAR operators according to the cut length
            for(int k = anonymizedCigarElements.size() - 1; k >= 0; k--){
                CigarElement cigarElement = anonymizedCigarElements.get(k);
                int cigarLength = cigarElement.getLength();
                int diff = cigarLength - cutLength;
                if(!cigarElement.getOperator().consumesReadBases()){
                    anonymizedCigarElements.remove(k);
                    continue;
                }
                if(diff > 0){
                    CigarOperator op = cigarElement.getOperator();
                    if(op == CigarOperator.I){
                        op = CigarOperator.S;
                    }
                    anonymizedCigarElements.set(k, new CigarElement(diff, op));
                    break;
                }
                else if(diff == 0){
                    anonymizedCigarElements.remove(k);
                    break;
                }
                else{
                    anonymizedCigarElements.remove(k);
                    cutLength = Math.abs(diff);
                }
            }
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
        for(PairCalledVariation indel : indelSimpleSignals) {
            PairCalledVariation.VariantType variantType = indel.getVariantType();
            if (PairCalledVariation.VariantType.DEL.equals(variantType)) {
                newSize += indel.getLength();
            }
            if (PairCalledVariation.VariantType.INS.equals(variantType)) {
                newSize -= indel.getLength();
            }
        }
        return newSize;
    }

    private byte[] processSNVoperations(int originalSeqLength) {
        byte[] snvOps = new byte[originalSeqLength];
        for (PairCalledVariation snv : SNVsimpleSignals){
            int snvOpPosition = snv.getInReadPosition(this);
            //Change to retrieve from memoized ref genome
            byte op = snv.getRefAllele()[0];
            snvOps[snvOpPosition] = op;
        }
        return snvOps;
    }

    private int[] processCigarOperations(int originalCigarLength) {
        int[] indelOps = new int[originalCigarLength];
        // A negative operation (op) value, causes an elimination of the signal, whereas a positive value generates
        // an additive change with base pair filling from the reference genome
        for (PairCalledVariation indel : indelSimpleSignals){
            int indelOpPos = indel.getInReadPosition(this);
            int op = PairCalledVariation.VariantType.INS == indel.getVariantType() ?
                    -(indel.getLength()) : indel.getLength();
            indelOps[indelOpPos] = op;
        }
        for (Signal signal : complexSignals){
            if(Signal.Source.SOFT_CLIP == signal.getSource()){
                int softClipSignalPos = signal.getInReadPosition();
                int op = signal.getLength();
                //SoftClips are always filled, but the filling is done either at the beginning or end of the read, such that
                //starting softclips are treated as deletions, and ending softclips as insertions
                indelOps[softClipSignalPos] = softClipSignalPos == 0 ? op : -op;
            }
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

    public void setSignalsToAnonymize(List<Signal> signals) {
        addAllSignalsToAnonymize(signals);
    }

    public boolean addAllSignalsToAnonymize(List<Signal> signals){
        boolean added = false;
        for (Signal signal : signals){
            added = addSignalToAnonymize(signal);
            if (!added) return false;
        }
        return added;
    }

    public boolean addSignalToAnonymize(Signal signal){
        // Can change when dealing with SVs
        boolean added = false;
        if (Signal.Source.SIMPLE_VARIATION == signal.getSource()){
            PairCalledVariation variation = signal.getCalledVariation();
            String varType = PairCalledVariation.VariantType.SNV.equals(variation.getVariantType()) ?
                    PairCalledVariation.GENERIC_TYPE_SNV : PairCalledVariation.GENERIC_TYPE_INDEL;
            if (PairCalledVariation.GENERIC_TYPE_SNV.equals(varType)) {
                added = SNVsimpleSignals.add(variation);
            }
            if(PairCalledVariation.GENERIC_TYPE_INDEL.equals(varType)){
                added = indelSimpleSignals.add(variation);
            }
        } else if (Signal.Source.SOFT_CLIP == signal.getSource()) {
            added = complexSignals.add(signal);
        }
        //TODO: Add other signal types
        return added;
    }

    public String getReadAlignmentId(){
        return readAlnId;
    }

    public String getSequenceName() {
        return readAlignment.getContig();
    }

    @Override
    public int getSequenceIdx() {
        return 0;
    }

    public int getStart() {
        return readAlignment.getStart();
    }

    public int getEnd() {
        return readAlignment.getEnd();
    }

    @Override
    public void setSequenceIdx(int sequenceIdx) {
    }

    public int getMappingQuality(){
        return readAlignment.getMappingQuality();
    }

    public boolean isSupplementary() {
        return readAlignment.getSupplementaryAlignmentFlag();
    }

    public boolean isReverse() {
        return readAlignment.getReadNegativeStrandFlag();
    }

    public void setReverse(boolean reverse) {
        readAlignment.setReadNegativeStrandFlag(reverse);
    }

    public boolean isAnonymized() {
        if(SNVsimpleSignals.isEmpty() && indelSimpleSignals.isEmpty() && complexSignals.isEmpty()){
            //This is a read that does not have to be anonymized
            isAnonymized = true;
        }
        return isAnonymized;
    }

    public boolean isFromNormalDataset(){
        return isNormalDataset;
    }

    public boolean isFromTumoralDataset(){
        return !isNormalDataset;
    }

    public String getReadName() {
        return readAlignment.getReadName();
    }

    public int getLength() {
        return readAlignment.getReadLength();
    }

    public boolean isPair1(){
        return readAlignment.getFirstOfPairFlag();
    }

    public boolean isPair2(){
        return readAlignment.getSecondOfPairFlag();
    }

    public int getPairIdx() {
        return isPair1() ? PAIR_1_IDX : PAIR_2_IDX;
    }

    public void setPairIdx(int pairIdx) {
        if(pairIdx == PAIR_1_IDX){
            readAlignment.setFirstOfPairFlag(true);
            readAlignment.setSecondOfPairFlag(false);
        }else if(pairIdx == PAIR_2_IDX){
            readAlignment.setFirstOfPairFlag(false);
            readAlignment.setSecondOfPairFlag(true);
        }
    }

    public byte[] getOriginalSequenceArray(){
        return readAlignment.getReadBases();
    }

    public byte[] getOriginalQualitiesArray(){
        return readAlignment.getBaseQualities();
    }

    public SAMRecord cloneRecord(){
        try {
            return (SAMRecord) readAlignment.clone();
        } catch (CloneNotSupportedException e) {
            e.printStackTrace();
            throw new RuntimeException("Error cloning record " + readAlignment.getReadName(), e);
        }
    }

    @Override
    public String toString(){
        StringBuilder builder = new StringBuilder();
        builder.append("ReadName=").append(getReadName());
        builder.append(" ReadID=").append(getReadAlignmentId());
        builder.append(" Start=").append(getStart());
        builder.append(" End=").append(getEnd());
        builder.append(" OrgSeq=").append(new String(getOriginalSequenceArray(), StandardCharsets.UTF_8));
        builder.append(" AnonSeq=").append(new String(anonymizedSequenceArray, StandardCharsets.UTF_8));
        builder.append(" OrgQual=").append(Arrays.toString(getOriginalQualitiesArray()));
        builder.append(" AnonQual=").append(Arrays.toString(anonymizedQualitiesArray));
        builder.append(" CIGAR=").append(anonymizedCigar.toString());
        builder.append(" TAGS=").append(readAlignment.getAttributes());
        return builder.toString();
    }

    //Provides a unique ID for each read alignment, directly from a SAMRecord
    public static String generateReadAlnId(SAMRecord alignment){
        StringBuilder builder = new StringBuilder();
        int pairIdx = alignment.getFirstOfPairFlag() ? PAIR_1_IDX : PAIR_2_IDX;
        builder.append(alignment.getReadName());
        builder.append(DEFAULT_ID_NAME_SEPARATOR);
        builder.append(pairIdx);
        builder.append(DEFAULT_ID_NAME_SEPARATOR);
        builder.append(alignment.getAlignmentStart());
        builder.append(alignment.getCigar().toString());
        builder.append(alignment.getBaseQualityString());
        builder.append(alignment.getReadString());
        return builder.toString();
    }

    @Override
    public int compareTo(GenomicRegion genomicRegion) {
        return compare(this, genomicRegion);
    }
}
