package analysis;
import genomicelements.CalledVariation;
import genomicelements.CalledVariation.SomaticVariationType;
import genomicelements.CalledVariation.VariantType;
import genomicelements.GenomicRegion;
import genomicelements.PairedPileup;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;
import io.SamplePairReadAlignmentReader;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.logging.Logger;
import static genomicelements.ShortReadAlignment.*;


/**
 * Class used to classify all variation from a paired normal-tumor sample from each pileup position
 * @author Nicolas Gaitan
 */

public class VariationClassifier {

    private static final Logger LOGGER = Logger.getLogger(VariationClassifier.class.getName());

    private static final double SOFTCLIP_DISTANCE_THRESHOLD = 10;

    public static final char NULL_BASE = 'N';
    public static final Set<Character> ALPHABET = new HashSet<>(
            Arrays.asList(
                    'A', 'T', 'C', 'G'
            )
    );

    public static final int SLIDING_WINDOW_LIMIT = 200;

    Set<String> readsToExclude;
    Map<String, List<CalledVariation>> potentialGermlinesPerRead;
    Map<String, Map<Integer,CalledVariation>> somaticVariantsToKeep;
    boolean diffuseIndelCalls;
    boolean anonymizePotentialLeaksInVCFSomatics;

    public VariationClassifier(){
        potentialGermlinesPerRead = new HashMap<>();
        diffuseIndelCalls = false;
        setAnonymizePotentialLeaks(false);
        readsToExclude = new HashSet<>();
        somaticVariantsToKeep = new HashMap<>();
    }

    public Map<String, List<CalledVariation>> getPotentialGermlinesPerRead() {
        return potentialGermlinesPerRead;
    }

    public void setReadsToExclude(Set<String> readsToExclude){
        this.readsToExclude = readsToExclude;
    }

    public void setVCFVariantsToKeep(Map<String, Map<Integer,CalledVariation>> somaticVariantsToKeep){
        this.somaticVariantsToKeep = somaticVariantsToKeep;
    }

    public void setAnonymizePotentialLeaks(boolean anonymize){
        this.anonymizePotentialLeaksInVCFSomatics = anonymize;
    }
    /**
     * Call for discovering variation over a specific genomic region
     * @param normalPath
     * @param tumorPath
     * @param refGenome
     * @param region
     * @throws IOException
     */
    public void callVariation(String normalPath, String tumorPath, String refGenome, GenomicRegion region) throws IOException {
        try(SamplePairReadAlignmentReader pairPileupReader = new SamplePairReadAlignmentReader(normalPath, tumorPath, refGenome, region);
            IndexedFastaSequenceFile referenceWalker = new IndexedFastaSequenceFile(new File(refGenome))){
            //Retrieve signals from their normal sample even if there is no coverage in the tumor sample
            pairPileupReader.setReturnNormal(true);
            callVariation(pairPileupReader, referenceWalker);
        }
    }

    public void callVariation(SamplePairReadAlignmentReader pairPileupReader, IndexedFastaSequenceFile referenceWalker){
        Map<Integer, List<CalledVariation>> variationPerPos = new HashMap<>();
        Set<String> seenReads = new HashSet<>();
        int p = 1;
        for (PairedPileup pileup : pairPileupReader){
            int pos = pileup.getReferencePos();
            classifyVariationInPairedPileup(variationPerPos, pileup, seenReads, referenceWalker);
            processPotentialGermlines(variationPerPos.get(pos));
            if (p==SLIDING_WINDOW_LIMIT) {
                p = 0;
//                if (diffuseIndelCalls){
//
//                }
            }
            variationPerPos.remove(pos-SLIDING_WINDOW_LIMIT);
            p++;
        }
    }

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
    }

    private void classifyPileupVariation(String sequenceName, int refPosition, List<RecordAndOffset> pileup, Set<String> seenReads, byte referenceBase,
                                         Map<Integer, List<CalledVariation>> variationPerPos, IndexedFastaSequenceFile referenceWalker, boolean isNormalDataset) {
        // In case single normal pileup is queried alone, to guarantee that null tumor pileups are not accessed
        if(pileup==null) return;
        // May be removing the read;pair name of seenReads after it reaches he last position in pileup
        // or adding the CIGAR to seen reads string
        for (RecordAndOffset pileupRecord : pileup){
            String readName = pileupRecord.getReadName();
            variationPerPos.computeIfAbsent(refPosition, v -> new ArrayList<>());
            if(readsToExclude.contains(readName)) continue;
            SAMRecord samRecord = pileupRecord.getRecord();
            //This may be extended to support other types of reads (e.g. long reads)
            int pairIdx = samRecord.getFirstOfPairFlag() ? PAIR_1_IDX : PAIR_2_IDX;
            // pairReadName represents the name of the read, the pair, and the reference position of the alignment
            //String pairReadName = getShortReadPairName(readName, pairIdx);
            //String pairReadName = getShortReadAlignmentId(readName, pairIdx, refPosition);
            String pairReadName = generateReadId(samRecord);
            // specificReadName represents the name of the read , pair, and which partial alignment (if any) it comes from
            // TODO: Delete the now unnecesary additional Id, merge with pairReadName
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
        for (int i = 0; i < cigarElems.size(); i++){
            CigarElement cigarElement = cigarElems.get(i);
            CigarOperator op = cigarElement.getOperator();
            if (op.isIndel()){
                int currentRefPos = initRefPos + currentCigarLength-1;
                //int inReadPos = readConsumedBaseNumber == 0 ? 1 : readConsumedBaseNumber;
                int inReadPos = samRecord.getReadPositionAtReferencePosition(currentRefPos);
                int length = cigarElement.getLength();
                VariantType indelType;
                int vcfStdEnd;
                int inRefend;
                int inReadEnd;
                byte[] altAllele;
                if (CigarOperator.I.equals(op)){
                    indelType = VariantType.INS;
                    inRefend = currentRefPos;// + 1;
                    vcfStdEnd = inRefend + 1;
                    inReadEnd = inReadPos + length + 1;
                    altAllele = new byte[1+length];
                }
                else{
                    indelType = VariantType.DEL;
                    inRefend = currentRefPos + length;
                    vcfStdEnd = inRefend;
                    inReadEnd = inReadPos + 1;
                    altAllele = new byte[1];
                }
                // Ends vary based on the functions to recover the alleles, whether they are inclusive or exclusive on interval ends
                //byte[] altAllele = Arrays.copyOfRange(sequenceBases, inReadPos-1, inReadEnd-1);
                byte[] refAllele = referenceWalker.getSubsequenceAt(sequenceName, currentRefPos, inRefend).getBases();
                altAllele[0] = refAllele[0];
                if (CigarOperator.I.equals(op)){
                    System.arraycopy(sequenceBases, inReadPos, altAllele, 1, altAllele.length - 1);
                }
                //altAllele = Arrays.copyOfRange(sequenceBases, inReadPos, inReadEnd);
                CalledVariation calledVar = new CalledVariation(sequenceName, currentRefPos, vcfStdEnd, indelType, length,
                        altAllele, refAllele);
                List<CalledVariation> variationInPos = variationPerPos.computeIfAbsent(currentRefPos, v -> new ArrayList<>());
                int indexSearch = variationInPos.indexOf(calledVar);
                boolean variationExists = indexSearch != -1;
                if (variationExists) calledVar = variationInPos.get(indexSearch);
                // Saves the CIGAR index of the INDEL signal
                calledVar.addSupportingRead(pairReadName, i);
                processSomaticType(variationInPos, calledVar, variationExists, isNormalDataset);
            }
            if(op.consumesReferenceBases()){
                currentCigarLength += cigarElement.getLength();
            }
            if(op.consumesReadBases()){
                readConsumedBaseNumber += cigarElement.getLength();
            }
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
        CalledVariation calledVar = new CalledVariation(sequenceName, refPosition, refPosition, VariantType.SNV, 1,
                altAllele, refAllele);
        List<CalledVariation> variationInPos = variationPerPos.computeIfAbsent(refPosition, v -> new ArrayList<>());
        int indexSearch = variationInPos.indexOf(calledVar);
        boolean variationExists = indexSearch != -1;
        if (variationExists) calledVar = variationInPos.get(indexSearch);
        //TODO: Check and fix inReadPosition if wrongly estimated in supplementaries
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

    /**
     * Retrieves the calls to be anonymized for each read alignment, uniquely by read name, pair and position
     * @param variationInPos
     */
    private void processPotentialGermlines(List<CalledVariation> variationInPos) {
        for (CalledVariation var : variationInPos){
            if(!anonymizeThisVariant(var)) continue;
            Map<String, Integer> supportingReads = var.getSupportingReads();
            for (Map.Entry<String, Integer> entry : supportingReads.entrySet()){
                String readAlnId = entry.getKey();
                List<CalledVariation> potentialGermlinesInReadAlignment = potentialGermlinesPerRead
                        .computeIfAbsent(readAlnId, v -> new ArrayList<>());
                potentialGermlinesInReadAlignment.add(var);
            }
        }
    }

    private boolean anonymizeThisVariant(CalledVariation variant) {
        boolean isValidatedSomatic = false;
        // Anonymize only potential germlines if seen in both datasets, at least once in each, or more than once if only found in the normal tissue mappings
        boolean isPotentialGermline = SomaticVariationType.TUMORAL_NORMAL_VARIANT.equals(variant.getSomaticVariationType()) ||
                SomaticVariationType.NORMAL_ONLY_VARIANT.equals(variant.getSomaticVariationType());
        //Avoid anonymizing somatic variants recorded in the VCF file
        if (!somaticVariantsToKeep.isEmpty()){
            if(somaticVariantsToKeep.containsKey(variant.getSeqName())){
                Map<Integer, CalledVariation> validatedSomaticsAtSeq = somaticVariantsToKeep.get(variant.getSeqName());
                if(validatedSomaticsAtSeq.containsKey(variant.getPos())){
                    CalledVariation validatedSomaticAtPos = validatedSomaticsAtSeq.get(variant.getPos());
                    if(validatedSomaticAtPos.equals(variant)) isValidatedSomatic = true;
                }
            }
        }
        if(isValidatedSomatic){
            return isPotentialGermline && anonymizePotentialLeaksInVCFSomatics;
        }
        return isPotentialGermline;
    }

    public void setDiffuseIndelCalls(boolean diffuseIndelCalls) {
        this.diffuseIndelCalls = diffuseIndelCalls;
    }


    public static String getSpecificShortReadPairName(SAMRecord samRec, int pairIdx) {
        if (samRec.getAttribute("SA") != null){
            return samRec.getReadName() + DEFAULT_ID_NAME_SEPARATOR + pairIdx + DEFAULT_ID_NAME_SEPARATOR + generateAlignmentHash(samRec);
        }
        else{
            return samRec.getReadName() + DEFAULT_ID_NAME_SEPARATOR + pairIdx;
        }
    }
}
