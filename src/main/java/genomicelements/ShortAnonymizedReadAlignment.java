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
    public static final String ORIGIN_PAIR_TAG = "op";

    private final SAMRecord readAlignment;
    private String readAlnId;
    private byte[] referenceContigSequence;
    private int alnStart;
    private byte averageBaseQuality;
    private boolean isAnonymized;
    private byte[] anonymizedSequenceArray;
    private byte[] anonymizedQualitiesArray;
    private boolean isNormalDataset;
    private boolean mateOriginalPosIsEqual;

    List<CigarElement> anonymizedCigarElements;
    private Cigar anonymizedCigar;
    private List<PairCalledVariation> SNVsimpleSignals;
    private List<PairCalledVariation> indelSimpleSignals;
    private List<Signal> complexSignals;

    //Anonymization behaviour modifiers
    private boolean fixOrientation = false;
    private boolean hasDestructiveSignal = false;
    private boolean hasChromChangeSignal = false;
    private boolean updateInfoForMate = false;
    private boolean extendLeft = false;

    public ShortAnonymizedReadAlignment(SAMRecord readAlignment) {
        this.readAlignment = readAlignment;
        this.readAlnId = generateReadAlnId(readAlignment);
        this.alnStart = getStart();
        this.isAnonymized = false;
        this.mateOriginalPosIsEqual = readAlignment.getMateAlignmentStart() == getStart();
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

    public SAMRecord getAnonymizedSamRecord() throws IllegalStateException {
        if (!isAnonymized()) {
            throw new IllegalStateException("Read is not anonymized, so the anonymized SAMRecord cannot be generated");
        }
        return readAlignment;
    }

    private void correctOrientation(SAMRecord answer) {
        boolean comesFirst = answer.getAlignmentStart() <= answer.getMateAlignmentStart();
        if (comesFirst){
            answer.setReadNegativeStrandFlag(false);
            answer.setMateNegativeStrandFlag(true);
            answer.setInferredInsertSize(Math.abs(answer.getInferredInsertSize()));
        }
        else {
            answer.setReadNegativeStrandFlag(true);
            answer.setMateNegativeStrandFlag(false);
            answer.setInferredInsertSize(-Math.abs(answer.getInferredInsertSize()));
        }
    }

    public SAMRecord getNewPair(int insertSize) {
        SAMRecord newPair = new SAMRecord(readAlignment.getHeader());
        //Determine pair number and gather related information
        int newPairIdx = 1 - getPairIdx();
        boolean newPairIsFirstOfPair = newPairIdx == PAIR_1_IDX;
        int newPairStart;
        int newPairLength = getLength();
        int distanceFromMate = insertSize - (getLength()+getLength());
        boolean newMapsFirst = hasChromChangeSignal ? readAlignment.getReadNegativeStrandFlag() : readAlignment.getMateAlignmentStart() <= readAlignment.getAlignmentStart();
        if(newMapsFirst) {
            newPairStart = Math.max(alnStart - distanceFromMate - newPairLength + 1, 1);
            if(newPairStart + newPairLength >= referenceContigSequence.length){
                newPairStart = referenceContigSequence.length - newPairLength;
            }
        }
        else {
            newPairStart = Math.min(getEnd() + distanceFromMate + 1, referenceContigSequence.length - newPairLength);
        }
        //Generate read sequence from the reference, according complying with the input insert size
        byte[] newPairSequenceArray = new byte[newPairLength];
        int refPos = newPairStart - 1;
        for(int i = 0; i < newPairLength; i++){
            newPairSequenceArray[i] = referenceContigSequence[refPos];
            refPos++;
        }
        //Generate read qualities from the average base quality of the mate
        byte[] newPairQualitiesArray = new byte[newPairLength];
        Arrays.fill(newPairQualitiesArray, averageBaseQuality);
        //Generate CIGAR
        List<CigarElement> newPairCigarElements = new ArrayList<>();
        newPairCigarElements.add(new CigarElement(newPairLength, CigarOperator.M));
        Cigar newPairCigar = new Cigar(newPairCigarElements);
        //Set mate information (this)
        newPair.setMateReferenceName(getSequenceName());
        newPair.setMateAlignmentStart(alnStart);
        newPair.setMateReferenceIndex(readAlignment.getReferenceIndex());
        //Set the new pair information
        newPair.setReadName(getReadName());
        newPair.setReferenceName(getSequenceName());
        newPair.setAlignmentStart(newPairStart);
        newPair.setReadBases(newPairSequenceArray);
        newPair.setBaseQualities(newPairQualitiesArray);
        newPair.setCigar(newPairCigar);
        newPair.setMappingQuality(getMappingQuality());
        newPair.setReferenceIndex(readAlignment.getReferenceIndex());
        //Set all the flags
        newPair.setFirstOfPairFlag(newPairIsFirstOfPair);
        newPair.setSecondOfPairFlag(!newPairIsFirstOfPair);
        newPair.setReadPairedFlag(true);
        newPair.setMateUnmappedFlag(false);
        newPair.setProperPairFlag(true);
        newPair.setDuplicateReadFlag(readAlignment.getDuplicateReadFlag());
        newPair.setReadFailsVendorQualityCheckFlag(readAlignment.getReadFailsVendorQualityCheckFlag());
        newPair.setSupplementaryAlignmentFlag(false);
        //Correct insert size
        int correctedInsertSize = newMapsFirst ?
                this.getEnd() - newPair.getStart() + 1:
                newPair.getEnd() - this.getStart() + 1;
        if(newMapsFirst){
            newPair.setInferredInsertSize(correctedInsertSize);
            readAlignment.setInferredInsertSize(-correctedInsertSize);
        }
        else{
            newPair.setInferredInsertSize(-correctedInsertSize);
            readAlignment.setInferredInsertSize(correctedInsertSize);
        }
        //Update this pair information based on new pair
        readAlignment.setProperPairFlag(true);
        readAlignment.setMateReferenceName(getSequenceName());
        readAlignment.setMateAlignmentStart(newPairStart);
        readAlignment.setMateReferenceIndex(readAlignment.getReferenceIndex());
        //Correct orientations
        correctOrientation(newPair);
        correctOrientation(readAlignment);
        //Set TAGs for new pair: NM, MD, AS, XN, RG
        SequenceUtil.calculateMdAndNmTags(newPair, referenceContigSequence, true, true);
        newPair.setAttribute("AS", newPairLength);
        newPair.setAttribute("RG", readAlignment.getAttribute("RG"));
        newPair.setAttribute("MQ", getMappingQuality());
        //Custom tag to identify the origin pair (will be removed before writing)
        newPair.setAttribute(ORIGIN_PAIR_TAG, getPairIdx());
        readAlignment.setAttribute(ORIGIN_PAIR_TAG, getPairIdx());
        //Set TAGs for this pair: NM, MD, AS, XN, RG
        return newPair;
    }

    public boolean fixOrientation() {
        return fixOrientation;
    }

    public void setReferenceContigSequence(byte[] referenceContigSequence) {
        // 0-based memoized reference sequence, corresponding to the contig to which this read is mapped
        this.referenceContigSequence = referenceContigSequence;
    }

    public void setFixOrientation(boolean fixOrientation) {
        this.fixOrientation = fixOrientation;
    }

    public void anonymizeRead() throws IllegalStateException{
        if(referenceContigSequence == null){
            throw new IllegalStateException("The ShortAnonymizedReadAlignment.anonymizeVariants was called" +
                    "without setting the contigReferenceSequence first. setContigReferenceSequence, should always" +
                    " be alled after the constructor.");
        }
        int originalSeqLength = getOriginalSequenceArray().length;
        int expectedSize = Math.max(originalSeqLength, estimateNewReadSize(originalSeqLength));
        anonymizedSequenceArray = new byte[expectedSize];
        anonymizedQualitiesArray = new byte[expectedSize];
        averageBaseQuality = getAverageOfBytes(getOriginalQualitiesArray());
        List<CigarElement> initialCigarElements = new ArrayList<>(readAlignment.getCigar().getCigarElements());
        //Except a starting PG SoftClip, alignment start is invariant, as the anonymized read will map to the exact start coordinate as the original
        //0 based coordinate
        int currentRefAlnPos = alnStart-1;
        //Should be the same if we fill or remove bases from the end of the read
        //Holds the position over the original read
        int i = 0;
        //Holds the position over the anonymized read
        int j = 0;
        //Holds the current index of the CIGAR alignment elements
        int c = 0;
        //SNV Operations (value: byte as new base) to perform using the original read index
        byte[] SNVops = processSNVoperations(originalSeqLength);
        //Indel Operations (value: number of bases to remove or add) to perform using the CIGAR index
        int[] operations = processComplexOperations(initialCigarElements);
        while(c < initialCigarElements.size()){
            CigarElement cigarElem = initialCigarElements.get(c);
            CigarOperator currentCigarOp = cigarElem.getOperator();
            int opLength = cigarElem.getLength();
            int op = operations[c];
            if(op > 0){
                makeAdditiveChange(j, currentRefAlnPos, op);
                j += op;
                currentRefAlnPos += op;
            }
            else if(op < 0){
                i += Math.abs(op);
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
            if(extendLeft){
                byte[] tempSequenceArray = new byte[anonymizedSequenceArray.length];
                byte[] tempQualitiesArray = new byte[anonymizedQualitiesArray.length];
                int m = tempSequenceArray.length - 1;
                for(int k = j-1; k >= 0; k--){
                    tempSequenceArray[m] = anonymizedSequenceArray[k];
                    tempQualitiesArray[m] = anonymizedQualitiesArray[k];
                    m--;
                }
                anonymizedSequenceArray = tempSequenceArray;
                anonymizedQualitiesArray = tempQualitiesArray;
                alnStart -= (m+1);
                makeAdditiveChange(0, alnStart - 1, m + 1, true);
            }
            else {
                makeAdditiveChange(j, currentRefAlnPos, expectedSize-j);
            }
        }
        // Cut the read if the new size is larger than the original
        if(anonymizedSequenceArray.length > originalSeqLength) cutRead(originalSeqLength);
        // Process the CIGAR elements to remove redundant ones
        generateDefinitiveCigar();
        //Set the read alignment start position if modified
        if(alnStart != getStart()) updateInfoForMate = true;
        //Anonymize the SAMRecord read alignment
        readAlignment.setAlignmentStart(alnStart);
        readAlignment.setReadBases(anonymizedSequenceArray);
        readAlignment.setBaseQualities(anonymizedQualitiesArray);
        readAlignment.setCigar(anonymizedCigar);
        SequenceUtil.calculateMdAndNmTags(readAlignment, referenceContigSequence, true, true);
        // Fix orientation if needed
        if(fixOrientation) correctOrientation(readAlignment);
        isAnonymized = true;
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

    private void makeAdditiveChange(int initPos, int refInitPos, int length) {
        makeAdditiveChange(initPos, refInitPos, length, false);
    }

    /**
     * Fills the read with bases from the reference genome, starting at the given position
     * @param initPos Position in the read to start filling
     * @param refInitPos Position in the reference genome to start filling
     * @param length Number of bases to fill
     */
    private void makeAdditiveChange(int initPos, int refInitPos, int length, boolean extendLeft) {
        int j = initPos;
        int r = refInitPos;
        for(int l = 0; l < length; l++){
            anonymizedSequenceArray[j] = referenceContigSequence[r];
            anonymizedQualitiesArray[j] = averageBaseQuality;
            j++;
            r++;
        }
        if(extendLeft){
            List<CigarElement> tempCigarElements = new ArrayList<>();
            tempCigarElements.add(new CigarElement(length, CigarOperator.M));
            tempCigarElements.addAll(anonymizedCigarElements);
            anonymizedCigarElements = tempCigarElements;
        }
        else{
            anonymizedCigarElements.add(new CigarElement(length, CigarOperator.M));
        }
    }

    private void cutRead(int originalSeqLength) {
        int cutLength = anonymizedSequenceArray.length - originalSeqLength;
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

    private int[] processComplexOperations(List<CigarElement> originalCigarElems) {
        int n = originalCigarElems.size();
        int[] operations = new int[n];
        // A negative operation (op) value, causes an elimination of the signal, whereas a positive value generates
        // an additive change with base pair filling from the reference genome
        for (PairCalledVariation indel : indelSimpleSignals){
            int indelOpPos = indel.getInReadPosition(this);
            int op = PairCalledVariation.VariantType.INS == indel.getVariantType() ?
                    -(indel.getLength()) : indel.getLength();
            operations[indelOpPos] = op;
        }
        for (Signal signal : complexSignals){
            if(Signal.Source.SOFT_CLIP == signal.getSource()){
                //SoftClips are always filled, but the filling is done either at the beginning or end of the read, such that
                //starting softclips are treated as deletions, and ending softclips as insertions
                int softClipSignalPos = signal.getInReadPosition();
                int op = signal.getLength();
                boolean mapsFirstInPair = readAlignment.getStart() <= readAlignment.getMateAlignmentStart();
                boolean isSoftClipAtStart = softClipSignalPos == 0;
                operations[softClipSignalPos] = -op;
                if (mapsFirstInPair){
                    extendLeft = true;
                    CigarElement firstCigarElement = readAlignment.getCigar().getFirstCigarElement();
                    if (!isSoftClipAtStart && CigarOperator.S == firstCigarElement.getOperator()) {
                        operations[0] = -firstCigarElement.getLength();
                    }
                }
                else{
                    extendLeft = false;
                    CigarElement lastCigarElement = readAlignment.getCigar().getLastCigarElement();
                    if (isSoftClipAtStart && CigarOperator.S == lastCigarElement.getOperator()){
                        operations[operations.length - 1] = -lastCigarElement.getLength();
                    }
                }
            }
            if(Signal.Source.STRAND_ORIENTATION == signal.getSource()){
                setFixOrientation(true);
            }
            if(signal.isDestructiveSignal()){
                if(Signal.Source.CHROM_CHANGE == signal.getSource()) hasChromChangeSignal = true;
                hasDestructiveSignal = true;
            }
        }
        return operations;
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
        }
        else {
            // Add complex signals: SOFT_CLIP, STRAND_ORIENTATION, INSERT_SIZE, CHROM_CHANGE
            added = complexSignals.add(signal);
        }
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

    public Cigar getCigar() {
        return this.anonymizedCigar;
    }

    @Override
    public void setSequenceIdx(int sequenceIdx) {
    }

    public int getMappingQuality(){
        return readAlignment.getMappingQuality();
    }

    public int getPairUpdatedPos() {
        if (!isAnonymized()) {
            throw new IllegalStateException("Read is not anonymized, so this position may not be updated");
        }
        return alnStart;
    }

    public boolean mateOriginalPosIsEqual() {
        return mateOriginalPosIsEqual;
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

    public boolean hasDestructiveSignal(){
        return hasDestructiveSignal;
    }

    public boolean updateInfoForMate(){
        return updateInfoForMate;
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
