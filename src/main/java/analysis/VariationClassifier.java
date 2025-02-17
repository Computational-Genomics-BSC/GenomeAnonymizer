package analysis;
import genomicelements.CalledVariation;
import genomicelements.CalledVariation.SomaticVariationType;
import genomicelements.CalledVariation.VariantType;
import genomicelements.GenomicRegion;
import genomicelements.PairedPileup;
import genomicelements.Signal;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;
import io.SamplePairReadAlignmentReader;
import utils.Operations;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.logging.Logger;
import java.lang.Math;
import static genomicelements.ShortReadAlignment.*;


/**
 * Class used to classify all variation from a paired normal-tumor sample from each pileup position
 * @author Nicolas Gaitan
 */

public class VariationClassifier {

    private static final Logger LOGGER = Logger.getLogger(VariationClassifier.class.getName());

    public static final char NULL_BASE = 'N';
    public static final Set<Character> ALPHABET = new HashSet<>(
            Arrays.asList(
                    'A', 'T', 'C', 'G'
            )
    );

    public static final int SLIDING_WINDOW_LIMIT = 200;
    public static final int MAX_LOCATION_DISTANCE_THRESHOLD = 400;
    public static final int SIGNAL_REGION_LIMIT = 10000;
    // Assuming a maximum position distance of 10, and 15 of length
    public static final int COMPLEX_SIGNAL_THRESHOLD = 18;


    private Set<String> readsToExclude;
    private int insertSizeMinThreshold;
    private int insertSizeMaxThreshold;
    private Map<String, List<Signal>> potentialGermlinesPerRead;
    private Map<String, Map<Integer,CalledVariation>> somaticVariantsToKeep;
    private boolean diffuseIndelCalls;
    private boolean anonymizePotentialLeaksInVCFSomatics;

    //DEBUG
    public Map<String, Long> METHOD_TIME_MAP = new HashMap<>();
    //DEBUG

    public VariationClassifier(){
        potentialGermlinesPerRead = new HashMap<>();
        diffuseIndelCalls = false;
        setAnonymizePotentialLeaks(false);
        readsToExclude = new HashSet<>();
        somaticVariantsToKeep = new HashMap<>();
    }

    public Map<String, List<Signal>> getPotentialGermlinesPerRead() {
        return potentialGermlinesPerRead;
    }

    public void setInsertSizeMinThreshold(int insertSizeMinThreshold) {
        this.insertSizeMinThreshold = insertSizeMinThreshold;
    }

    public void setInsertSizeMaxThreshold(int insertSizeMaxThreshold) {
        this.insertSizeMaxThreshold = insertSizeMaxThreshold;
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
        List<Signal> signalsInRegion = new ArrayList<>();
        Iterator<PairedPileup> pileupIterator = pairPileupReader.iterator();
        while (pileupIterator.hasNext()) {
            PairedPileup pileup = pileupIterator.next();
            int pos = pileup.getReferencePos();
            long startclassifyVariationInPairedPileup = System.currentTimeMillis();
            classifyVariationInPairedPileup(variationPerPos, signalsInRegion, pileup, seenReads, referenceWalker);
            long endclassifyVariationInPairedPileup = System.currentTimeMillis();
            METHOD_TIME_MAP.compute("classifyVariationInPairedPileup", (k,v) -> v == null ?
                    endclassifyVariationInPairedPileup-startclassifyVariationInPairedPileup :
                    v + endclassifyVariationInPairedPileup-startclassifyVariationInPairedPileup);
            long startprocessSimpleSignals = System.currentTimeMillis();
            processSimpleSignals(variationPerPos.get(pos));
            long endprocessSimpleSignals = System.currentTimeMillis();
            METHOD_TIME_MAP.compute("processSimpleSignals",  (k,v) -> v == null ?
                    endprocessSimpleSignals-startprocessSimpleSignals :
                    v + endprocessSimpleSignals-startprocessSimpleSignals);
            if (!pileupIterator.hasNext() || signalsInRegion.size() >= SIGNAL_REGION_LIMIT) {
                long startprocessComplexSignals = System.currentTimeMillis();
                processComplexSignals(signalsInRegion);
                long endprocessComplexSignals = System.currentTimeMillis();
                METHOD_TIME_MAP.compute("processComplexSignals",  (k,v) -> v == null ?
                        endprocessComplexSignals-startprocessComplexSignals :
                        v + endprocessComplexSignals-startprocessComplexSignals);
//                diffuseIndelCalls(); -> here?
                signalsInRegion = new ArrayList<>();
                p = 0;
            }
            variationPerPos.remove(pos - SLIDING_WINDOW_LIMIT);
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
    public void classifyVariationInPairedPileup(Map<Integer, List<CalledVariation>> variationPerPos, List<Signal> signalsInRegion,
                                                PairedPileup pairedPileup, Set<String> seenReads, IndexedFastaSequenceFile referenceWalker){
        // Map<Integer, List<CalledVariation>> variationPerPos = new HashMap<>();
        List<RecordAndOffset> normalPileup = pairedPileup.getNormalPileup();
        List<RecordAndOffset> tumorPileup = pairedPileup.getTumorPileup();
        String contig = pairedPileup.getRefenceSequenceName();
        int refPos = pairedPileup.getReferencePos();
        byte refBase = referenceWalker.getSubsequenceAt(contig, refPos, refPos).getBases()[0];
        classifyPileupVariation(contig, refPos, normalPileup, seenReads, refBase, variationPerPos, signalsInRegion, referenceWalker,true);
        classifyPileupVariation(contig, refPos, tumorPileup, seenReads, refBase, variationPerPos, signalsInRegion, referenceWalker,false);
    }

    private void classifyPileupVariation(String sequenceName, int refPosition, List<RecordAndOffset> pileup, Set<String> seenReads, byte referenceBase,
                                         Map<Integer, List<CalledVariation>> variationPerPos, List<Signal> signalsInRegion,
                                         IndexedFastaSequenceFile referenceWalker, boolean isNormalDataset) {
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
                long startdiscoverIndelsAndSignalsFromCIGAR = System.currentTimeMillis();
                discoverIndelsAndSignalsFromCIGAR(samRecord, pairReadName, variationPerPos, signalsInRegion, referenceWalker, isNormalDataset);
                long enddiscoverIndelsAndSignalsFromCIGAR = System.currentTimeMillis();
                METHOD_TIME_MAP.compute("discoverIndelsAndSignalsFromCIGAR",  (k,v) -> v == null ?
                        enddiscoverIndelsAndSignalsFromCIGAR-startdiscoverIndelsAndSignalsFromCIGAR :
                        v + enddiscoverIndelsAndSignalsFromCIGAR-startdiscoverIndelsAndSignalsFromCIGAR);
                // One signal per mate
                complexSignalsFromMates(samRecord, signalsInRegion);
                seenReads.add(specificReadName);
            }
            int inReadPosition = samRecord.getReadPositionAtReferencePosition(refPosition);
            if (inReadPosition==0) continue;
            char referenceBaseUpper = Character.toUpperCase((char) referenceBase);
            char readBaseUpper = Character.toUpperCase((char) pileupRecord.getReadBase());
            long startdiscoverSNVs = System.currentTimeMillis();
            discoverSNVs(pairReadName, variationPerPos, readBaseUpper, referenceBaseUpper, sequenceName, refPosition,
                    inReadPosition, isNormalDataset);
            long enddiscoverSNVs = System.currentTimeMillis();
            METHOD_TIME_MAP.compute("discoverSNVs",  (k,v) -> v == null ?
                    enddiscoverSNVs-startdiscoverSNVs :
                    v + enddiscoverSNVs-startdiscoverSNVs);
        }
    }

    public void discoverIndelsAndSignalsFromCIGAR(SAMRecord samRecord, String pairReadName, Map<Integer, List<CalledVariation>> variationPerPos,
                                                  List<Signal> signalsInRegion, IndexedFastaSequenceFile referenceWalker, boolean isNormalDataset){
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
                int inReadPos = samRecord.getReadPositionAtReferencePosition(currentRefPos);
                int length = cigarElement.getLength();
                VariantType indelType;
                int vcfStdEnd;
                int inRefend;
                int inReadEnd;
                byte[] altAllele;
                if (CigarOperator.I == op){
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
                byte[] refAllele = referenceWalker.getSubsequenceAt(sequenceName, currentRefPos, inRefend).getBases();
                altAllele[0] = refAllele[0];
                if (CigarOperator.I.equals(op)){
                    System.arraycopy(sequenceBases, inReadPos, altAllele, 1, altAllele.length - 1);
                }
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
            if(op.isClipping()){
                //int currentRefPos = initRefPos + currentCigarLength-1;
                int currentRefPos = currentCigarLength == 0 ? initRefPos : initRefPos + currentCigarLength-1;
                int inReadPos = samRecord.getReadPositionAtReferencePosition(currentRefPos);
                int length = cigarElement.getLength();
                if(CigarOperator.S == op){
                    Signal calledSignal = new Signal(sequenceName, currentRefPos, generateReadId(samRecord), i, length,
                            Signal.Source.SOFT_CLIP);
                    calledSignal.setIsFromNormalDataset(isNormalDataset);
                    signalsInRegion.add(calledSignal);
                }
            }
            if(op.consumesReferenceBases()){
                currentCigarLength += cigarElement.getLength();
            }
            if(op.consumesReadBases()){
                readConsumedBaseNumber += cigarElement.getLength();
            }
        }
    }

    public void complexSignalsFromMates(SAMRecord samRecord, List<Signal> signalsInRegion) {
        // Check if reads are in different chromosomes
        if (!samRecord.getReferenceIndex().equals(samRecord.getMateReferenceIndex())) {
            Signal calledSignal = new Signal(samRecord.getContig(), samRecord.getAlignmentStart(), generateReadId(samRecord), 0, 0, Signal.Source.CHROM_CHANGE);
            signalsInRegion.add(calledSignal);
            return;
        }
        // Check signal strands: FF, RF and RR
        // Check if this is the first or second pair (assume both are mapped)
        boolean firstRead = samRecord.getAlignmentStart() <= samRecord.getMateAlignmentStart();
        boolean firstForward, secondForward;
        if (firstRead) {
            firstForward = !samRecord.getReadNegativeStrandFlag();
            secondForward = !samRecord.getMateNegativeStrandFlag();
        } else {
            firstForward = !samRecord.getMateNegativeStrandFlag();
            secondForward = !samRecord.getReadNegativeStrandFlag();
        }
        int insertSize = Math.abs(samRecord.getInferredInsertSize());
        if ((firstForward && secondForward) || (!firstForward && !secondForward) || (!firstForward && secondForward)) {
            Signal calledSignal = new Signal(samRecord.getContig(), samRecord.getAlignmentStart(), generateReadId(samRecord), 0, insertSize, Signal.Source.STRAND_ORIENTATION);
            signalsInRegion.add(calledSignal);
            return;
        }
        // Check insert size
        if (insertSize < insertSizeMinThreshold || insertSize > insertSizeMaxThreshold) {
            Signal calledSignal = new Signal(samRecord.getContig(), samRecord.getAlignmentStart(), generateReadId(samRecord), 0, insertSize, Signal.Source.INSERT_SIZE);
            signalsInRegion.add(calledSignal);
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
    private void processSimpleSignals(List<CalledVariation> variationInPos) {
        for (CalledVariation var : variationInPos){
            if(!isPotentialGermline(var)) continue;
            Map<String, Integer> supportingReads = var.getSupportingReads();
            for (Map.Entry<String, Integer> entry : supportingReads.entrySet()){
                String readAlnId = entry.getKey();
                List<Signal> potentialGermlinesInReadAlignment = potentialGermlinesPerRead
                        .computeIfAbsent(readAlnId, v -> new ArrayList<>());
                potentialGermlinesInReadAlignment.add(new Signal(readAlnId, var));
            }
        }
    }

    /**
     * Classify complex signals into potential germline variations by virtually inferring either a complete graph of
     * signals from the normal dataset, or a bipartite graph from the normal-tumor pair
     *
     * @param signalsInRegion
     */
    private void processComplexSignals(List<Signal> signalsInRegion) {
//        List<List<Integer>> adjacencyGraph = new ArrayList<>();
//        Collections.fill(adjacencyGraph, new ArrayList<>());
        int n = signalsInRegion.size();
        List<Signal> currentPartition = new ArrayList<>();
        for(int i = 0; i < n-1; i++){
            Signal currentSignal = signalsInRegion.get(i);
            currentPartition.add(currentSignal);
            Signal nextSignal = signalsInRegion.get(i+1);
            int locationDistance = Math.abs(nextSignal.getLocation() - currentSignal.getLocation());
            boolean lastSignalUnreachable = locationDistance > MAX_LOCATION_DISTANCE_THRESHOLD;
            if(lastSignalUnreachable || currentPartition.size() >= SIGNAL_REGION_LIMIT  || (i==n-2)){
                if(!lastSignalUnreachable){
                    currentPartition.add(nextSignal);
                    i++;
                }
                processPartition(currentPartition);
                currentPartition = new ArrayList<>();
            }
        }
    }

    private void processPartition(List<Signal> signals) {
        int n = signals.size();
        boolean[] isClassifiedPG = new boolean[n];
        for (int i = 0; i < n; i++){
            Signal firstSignal = signals.get(i);
            for (int j = i + 1; j < n; j++){
                //if (i == j) continue;
                Signal secondSignal = signals.get(j);
                double signalDistance = Operations.computeTwoDimEuclideanDistance(firstSignal.getLocation(), secondSignal.getLocation(),
                        firstSignal.getLength(), secondSignal.getLength());
                if (signalDistance > COMPLEX_SIGNAL_THRESHOLD ||
                        firstSignal.isFromTumoralDataset() && secondSignal.isFromTumoralDataset()) continue;
                long startclassifyPGcomplexSignal = System.currentTimeMillis();
                classifyPGcomplexSignal(firstSignal, i, isClassifiedPG);
                classifyPGcomplexSignal(secondSignal, j, isClassifiedPG);
                long endclassifyPGcomplexSignal = System.currentTimeMillis();
                METHOD_TIME_MAP.compute("classifyPGcomplexSignal",  (k,v) -> v == null ?
                        endclassifyPGcomplexSignal-startclassifyPGcomplexSignal :
                        v + endclassifyPGcomplexSignal-startclassifyPGcomplexSignal);
//                adjacencyGraph.get(i).add(j);
//                adjacencyGraph.get(j).add(i);
            }
        }
    }

    private void classifyPGcomplexSignal(Signal signal, int pos, boolean[] isClassifiedPG) {
        if(!isClassifiedPG[pos]) {
            String readAlnId = signal.getReadAlnName();
            List<Signal> potentialGermlinesInReadAlignment = potentialGermlinesPerRead
                    .computeIfAbsent(readAlnId, v -> new ArrayList<>());
            potentialGermlinesInReadAlignment.add(signal);
            isClassifiedPG[pos] = true;
        }
    }

    private boolean isPotentialGermline(CalledVariation variation) {
        boolean isValidatedSomatic = false;
        // Anonymize only potential germlines if seen in both datasets, at least once in each, or more than once if only found in the normal tissue mappings
        boolean isPotentialGermline = SomaticVariationType.TUMORAL_NORMAL_VARIANT.equals(variation.getSomaticVariationType()) ||
                SomaticVariationType.NORMAL_ONLY_VARIANT.equals(variation.getSomaticVariationType());
        //Avoid anonymizing somatic variants recorded in the VCF file
        if (!somaticVariantsToKeep.isEmpty()){
            if(somaticVariantsToKeep.containsKey(variation.getSeqName())){
                Map<Integer, CalledVariation> validatedSomaticsAtSeq = somaticVariantsToKeep.get(variation.getSeqName());
                if(validatedSomaticsAtSeq.containsKey(variation.getPos())){
                    CalledVariation validatedSomaticAtPos = validatedSomaticsAtSeq.get(variation.getPos());
                    if(validatedSomaticAtPos.equals(variation)) isValidatedSomatic = true;
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
            return samRec.getReadName() + DEFAULT_ID_NAME_SEPARATOR + pairIdx + DEFAULT_ID_NAME_SEPARATOR + generateComplement(samRec);
        }
        else{
            return samRec.getReadName() + DEFAULT_ID_NAME_SEPARATOR + pairIdx;
        }
    }
}
