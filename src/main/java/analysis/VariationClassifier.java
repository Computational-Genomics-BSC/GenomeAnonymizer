package analysis;
import genomicelements.CalledVariation;
import genomicelements.CalledVariation.SomaticVariationType;
import genomicelements.CalledVariation.VariantType;
import genomicelements.PairedPileup;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;
import genomicelements.ShortAnonymizedRead;

import java.util.*;

import static genomicelements.ShortAnonymizedReadPair.*;

/**
 * Class that offers methods to classify all variation from a sample over pileup positions
 * @author Nicolás Gaitán
 */

public class VariationClassifier {

    public static final char NULL_BASE = 'N';
    public static final Set<Character> ALPHABET = new HashSet<>(
            Arrays.asList(
                    'A', 'T', 'C', 'G'
            )
    );
    /*public static final Set<CigarOperator> REF_CONSUMING_OPS = new HashSet<>(
            Arrays.asList(
                    CigarOperator.M, CigarOperator.D, CigarOperator.N, CigarOperator.EQ, CigarOperator.X
            )
    );*/

    public VariationClassifier(){
    }

    //TODO: Flush all variants from the variationPerPos map in the calling method at certain intervals, also diffuse
    /**
     * @param variationPerPos Map where keys are coordinates of variants, and the variants are values of different types
     * @param pairedPileup
     * @param seenReads
     * @param referenceWalker
     *
     */
    public void classifyVariationInPairedPileup(Map<Integer, List<CalledVariation>> variationPerPos, PairedPileup pairedPileup, Set<String> seenReads,
                                                                         IndexedFastaSequenceFile referenceWalker){
        // Map<Integer, List<CalledVariation>> variationPerPos = new HashMap<>();
        List<RecordAndOffset> normalPileup = pairedPileup.getNormalPileup();
        List<RecordAndOffset> tumorPileup = pairedPileup.getTumorPileup();
        String contig = pairedPileup.getRefenceSequenceName();
        int refPos = pairedPileup.getReferencePos();
        byte refBase = referenceWalker.getSubsequenceAt(contig, refPos, refPos).getBases()[0];
        classifyPileupVariation(contig, refPos, normalPileup, seenReads, refBase, variationPerPos, referenceWalker,true);
        classifyPileupVariation(contig, refPos, tumorPileup, seenReads, refBase, variationPerPos, referenceWalker,false);
        // return variationPerPos;
    }

    private void classifyPileupVariation(String sequenceName, int refPosition, List<RecordAndOffset> pileup, Set<String> seenReads, byte referenceBase,
                                         Map<Integer, List<CalledVariation>> variationPerPos, IndexedFastaSequenceFile referenceWalker, boolean isNormalDataset) {
        // May be removing the read;pair name of seenReads after it reaches he last position in pileup
        // or adding the CIGAR to seen reads string
        for (RecordAndOffset pileupRecord : pileup){
            String readName = pileupRecord.getReadName();
            SAMRecord samRecord = pileupRecord.getRecord();
            variationPerPos.computeIfAbsent(refPosition, v -> new ArrayList<>());
            if (samRecord.getReadUnmappedFlag()) continue;
            //This may be extended to support other types of reads (e.g. long reads)
            int pairIdx = samRecord.getFirstOfPairFlag() ? PAIR_1_IDX : PAIR_2_IDX;
            // pairReadName represents the name of the read and which pair it is only
            String pairReadName = getShortReadPairName(readName, pairIdx);
            // specificReadName represents the name of the read , pair, and which partial alignment (if any) it comes from
            String specificReadName = getSpecificShortReadPairName(samRecord, pairIdx);
            //
            if (!seenReads.contains(specificReadName)){
                discoverIndels(samRecord, pairReadName, variationPerPos, referenceWalker, isNormalDataset);
                seenReads.add(specificReadName);
            }
            int inReadPosition = samRecord.getReadPositionAtReferencePosition(refPosition);
            if (inReadPosition==0) continue;
            char referenceBaseUpper = Character.toUpperCase((char) referenceBase);
            char readBaseUpper = Character.toUpperCase((char) pileupRecord.getReadBase());
            discoverSNVs(pairReadName, variationPerPos, readBaseUpper, referenceBaseUpper, sequenceName, refPosition,
                    inReadPosition, isNormalDataset);
        }
    }

    public void discoverIndels(SAMRecord samRecord, String pairReadName, Map<Integer, List<CalledVariation>> variationPerPos,
                               IndexedFastaSequenceFile referenceWalker, boolean isNormalDataset){
        List<CigarElement> cigarElems = samRecord.getCigar().getCigarElements();
        int initRefPos = samRecord.getAlignmentStart();
        int currentCigarLength = 0;
        int readConsumedBaseNumber = 0;
        String sequenceName = samRecord.getContig();
        byte[] sequenceBases = samRecord.getReadBases();
        for (CigarElement cigarElement : cigarElems){
            CigarOperator op = cigarElement.getOperator();
            if (op.isIndel()){
                int currentRefPos = initRefPos + currentCigarLength;
                int inReadPos = readConsumedBaseNumber;// + currentCigarLength;
                int length = cigarElement.getLength();
                VariantType indelType;
                int vcfStdEnd;
                int inRefend;
                int inReadEnd;
                if (CigarOperator.I.equals(op)){
                    indelType = VariantType.INS;
                    inRefend = currentRefPos;// + 1;
                    vcfStdEnd = inRefend + 1;
                    inReadEnd = inReadPos + length;
                }
                else{
                    indelType = VariantType.DEL;
                    inRefend = currentRefPos + length - 1;
                    vcfStdEnd = inRefend;
                    inReadEnd = inReadPos + 1;
                }
                // Ends vary based on the functions to recover the alleles, whether they are inclusive or exclusive on interval ends
                byte[] altAllele = Arrays.copyOfRange(sequenceBases, inReadPos, inReadEnd);
                byte[] refAllele = referenceWalker.getSubsequenceAt(sequenceName, currentRefPos, inRefend).getBases();
                CalledVariation calledVar = new CalledVariation(sequenceName, currentRefPos, vcfStdEnd, indelType, length,
                        altAllele, refAllele);
                List<CalledVariation> variationInPos = variationPerPos.computeIfAbsent(currentRefPos, v -> new ArrayList<>());
                int indexSearch = variationInPos.indexOf(calledVar);
                boolean variationExists = indexSearch != -1;
                if (variationExists) calledVar = variationInPos.get(indexSearch);
                calledVar.addSupportingRead(pairReadName, inReadPos);
                processSomaticType(variationInPos, calledVar, variationExists, isNormalDataset);
            }
            // Fix error with estimating ends of slices for read sequence
            if(op.consumesReferenceBases()){
                currentCigarLength += cigarElement.getLength();
            }
            if(op.consumesReadBases()){
                readConsumedBaseNumber += cigarElement.getLength();
            }
            // if(CigarOperator.D.equals(op)){
            //    readConsumedBaseNumber -= cigarElement.getLength();
            //}
        }
    }

    private void discoverSNVs(String pairReadName, Map<Integer, List<CalledVariation>> variationPerPos,
                              char readBase, char referenceBase, String sequenceName, int refPosition, int inReadPosition,
                              boolean isNormalDataset) {
        if (readBase == NULL_BASE || readBase == referenceBase || !ALPHABET.contains(referenceBase)) return;
        byte[] altAllele = new byte[1];
        altAllele[0] = (byte) readBase;
        byte[] refAllele = new byte[1];
        refAllele[0] = (byte) referenceBase;
        CalledVariation calledVar = new CalledVariation(sequenceName, refPosition, refPosition, CalledVariation.VariantType.SNV, 1,
                altAllele, refAllele);
        List<CalledVariation> variationInPos = variationPerPos.computeIfAbsent(refPosition, v -> new ArrayList<>());
        int indexSearch = variationInPos.indexOf(calledVar);
        boolean variationExists = indexSearch != -1;
        if (variationExists) calledVar = variationInPos.get(indexSearch);
        calledVar.addSupportingRead(pairReadName, inReadPosition);
        processSomaticType(variationInPos, calledVar, variationExists, isNormalDataset);
    }

    private void processSomaticType(List<CalledVariation> variationInPos, CalledVariation calledVar, boolean variationExists, boolean isNormalDataset) {
        if (!variationExists){
            if (isNormalDataset){
                calledVar.setSomaticVariationType(SomaticVariationType.NORMAL_SINGLE_READ_VARIANT);
            }
            else{
                calledVar.setSomaticVariationType(SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT);
            }
            variationInPos.add(calledVar);
        }
        else{
            SomaticVariationType calledVarType = calledVar.getSomaticVariationType();
            if(isNormalDataset){
                if (SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT.equals(calledVarType) || SomaticVariationType.TUMORAL_ONLY_VARIANT.equals(calledVarType)){
                    calledVar.setSomaticVariationType(SomaticVariationType.TUMORAL_NORMAL_VARIANT);
                }
                if (SomaticVariationType.NORMAL_SINGLE_READ_VARIANT.equals(calledVarType)){
                    calledVar.setSomaticVariationType(SomaticVariationType.NORMAL_ONLY_VARIANT);
                }
            }
            else{
                if (SomaticVariationType.NORMAL_SINGLE_READ_VARIANT.equals(calledVarType) || SomaticVariationType.NORMAL_ONLY_VARIANT.equals(calledVarType)){
                    calledVar.setSomaticVariationType(SomaticVariationType.TUMORAL_NORMAL_VARIANT);
                }
                if (SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT.equals(calledVarType)){
                    calledVar.setSomaticVariationType(SomaticVariationType.TUMORAL_ONLY_VARIANT);
                }
            }
        }
    }

    public static String getSpecificShortReadPairName(SAMRecord samRec, int pairIdx) {
        if (samRec.getAttribute("SA") != null){
            return samRec.getReadName() + READ_PAIR_NAME_SEPARATOR + pairIdx + READ_PAIR_NAME_SEPARATOR + generateAlignmentHash(samRec);
        }
        else{
            return samRec.getReadName() + READ_PAIR_NAME_SEPARATOR + pairIdx;
        }
    }

    public static String getShortReadPairName(String readName, int pairIdx) {
        return readName + READ_PAIR_NAME_SEPARATOR + pairIdx;

    }

    private static String generateAlignmentHash(SAMRecord samRec) {
        long answer = 17;
        answer = 37*answer + samRec.getCigar().toString().hashCode(); //samRec.getCigarString().hashCode();
        answer = 37*answer + samRec.getBaseQualityString().hashCode();
        answer = 37*answer + samRec.getReadString().hashCode();
        return String.valueOf(answer);
    }
}
