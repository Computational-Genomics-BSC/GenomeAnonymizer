package analysis;
import genomicelements.*;
import genomicelements.PairCalledVariation.SomaticVariationType;
import genomicelements.PairCalledVariation.VariantType;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import io.SamplePairReadAlignmentReader;
import utils.MapCacheFIFO;
import utils.Operations;

import java.io.Closeable;
import java.io.IOException;
import java.util.*;
import java.util.function.Consumer;
import java.util.logging.Logger;

import static analysis.GenomeAnonymizer.DEFAULT_MIN_MAPPING_QUALITY;
import static utils.Operations.overlap;


/**
 * Class used to classify all variation from a paired normal-tumor sample from each pileup position, and provide
 * variated and not read alignments
 * @author Nicolas Gaitan
 */

public class AnonymizedReadAlignmentProvider implements Iterable<AnonymizedRead>, Closeable {

    private static final Logger LOGGER = Logger.getLogger(AnonymizedReadAlignmentProvider.class.getName());

    public static final List<Byte> ALPHABET_AS_BYTES = Arrays.asList(
            (byte) 'A', (byte) 'T', (byte) 'C', (byte) 'G'
    );

    public static final int SHORT_READ_LEFT_PILEUP_REGION_EXTENSION = 500;
    public static final int SLIDING_WINDOW_LIMIT = 200;
    public static final int MAX_LOCATION_DISTANCE_THRESHOLD = 400;
    public static final int SIGNAL_PER_REGION_LIMIT = 1000;

    // Assuming a maximum position distance of 10, and 15 of length
    public static final int COMPLEX_SIGNAL_THRESHOLD = 18;

    SamplePairReadAlignmentReader pairPileupReader;
    private int currentPileupPosition = 0;

    private Queue<AnonymizedRead> anonymizedReadQueue;
    private MapCacheFIFO<String, AnonymizedRead> anonymizedReadCache;
    private Set<String> onHoldReads;

    private Map<Integer, List<PairCalledVariation>> snvsPerPos;
    private Map<Integer, List<PairCalledVariation>> indelsPerPos;
    private List<Signal> complexSignals;

    private Set<String> readsToExclude;
    private Map<String, Map<Integer, PairCalledVariation>> somaticVariantsToKeep;
    private boolean diffuseIndelCalls;
    private boolean anonymizePotentialLeaksInVCFSomatics;

    private byte[] refSequence;
    private GenomicRegion genomicRegion;

    private int minMappingQuality = DEFAULT_MIN_MAPPING_QUALITY;

    private int numProcessedPileups = 0;

    //DEBUG
    public Map<String, Long> METHOD_TIME_MAP = new HashMap<>();
    //DEBUG

    public AnonymizedReadAlignmentProvider(){
        anonymizedReadQueue = new LinkedList<>();
        anonymizedReadCache = new MapCacheFIFO<>(10_000_000);
        onHoldReads = new HashSet<>();
        diffuseIndelCalls = false;
        setAnonymizePotentialLeaks(false);
        readsToExclude = new HashSet<>();
        somaticVariantsToKeep = new HashMap<>();
        snvsPerPos = new HashMap<>();
        indelsPerPos = new HashMap<>();
        complexSignals = new ArrayList<>();
        refSequence = new byte[0];
    }

    public void setReadsToExclude(Set<String> readsToExclude){
        this.readsToExclude = readsToExclude;
    }

    public void setVCFVariantsToKeep(Map<String, Map<Integer, PairCalledVariation>> somaticVariantsToKeep){
        this.somaticVariantsToKeep = somaticVariantsToKeep;
    }

    public void setMinMappingQuality(int minMappingQuality) {
        this.minMappingQuality = minMappingQuality;
    }

    public void setAnonymizePotentialLeaks(boolean anonymize){
        this.anonymizePotentialLeaksInVCFSomatics = anonymize;
    }

    public void setRefSequence(byte[] refSequence) {
        this.refSequence = refSequence;
    }

    /**
     * Call for discovering variation over a specific genomic region
     * @param normalPath
     * @param tumorPath
     * @param refGenome
     * @param region
     * @throws IOException
     */
    public void init(String normalPath, String tumorPath, String refGenome, GenomicRegion region) throws IOException {
        this.genomicRegion = region;
        int leftLimit = Math.max(genomicRegion.getStart() - SHORT_READ_LEFT_PILEUP_REGION_EXTENSION, 1);
        GenomicRegion leftExtendedRegion = new GenomicRegionBaseImpl(genomicRegion.getSequenceName(),
                leftLimit, genomicRegion.getEnd());
        leftExtendedRegion.setSequenceIdx(region.getSequenceIdx());
        pairPileupReader = new SamplePairReadAlignmentReader(normalPath, tumorPath, refGenome, refSequence, leftExtendedRegion);
        //Remove uncovered reads
        pairPileupReader.setRemoveUncovered(true);
        pairPileupReader.setIncludeDuplicates(true);
        pairPileupReader.setReadsToExclude(readsToExclude);
    }

    public void processNextPairedPileup(PairedPileup pileup){
        int refPosition = pileup.getReferencePos();
        currentPileupPosition = refPosition;
        initializeSignalsInPos(refPosition);
        long startclassifyVariationInPairedPileup = System.currentTimeMillis();
        classifyVariationInPairedPileup(pileup);
        long endclassifyVariationInPairedPileup = System.currentTimeMillis();
        METHOD_TIME_MAP.compute("classifyVariationInPairedPileup", (k,v) -> v == null ?
                endclassifyVariationInPairedPileup-startclassifyVariationInPairedPileup :
                v + endclassifyVariationInPairedPileup-startclassifyVariationInPairedPileup);
        long startprocessSimpleSignals = System.currentTimeMillis();
        processSimpleSignals(refPosition, snvsPerPos);
        processSimpleSignals(refPosition, indelsPerPos);
        long endprocessSimpleSignals = System.currentTimeMillis();
        METHOD_TIME_MAP.compute("processSimpleSNVSignals",  (k,v) -> v == null ?
                endprocessSimpleSignals-startprocessSimpleSignals :
                v + endprocessSimpleSignals-startprocessSimpleSignals);
        if (genomicRegion.getEnd() == currentPileupPosition || complexSignals.size() >= SIGNAL_PER_REGION_LIMIT) {
            long startprocessComplexSignals = System.currentTimeMillis();
            processComplexSignals();
            long endprocessComplexSignals = System.currentTimeMillis();
            METHOD_TIME_MAP.compute("processComplexSignals",  (k,v) -> v == null ?
                    endprocessComplexSignals-startprocessComplexSignals :
                    v + endprocessComplexSignals-startprocessComplexSignals);
//                diffuseIndelCalls(); -> here?
            complexSignals = new ArrayList<>();
        }
        clearUnusedSimpleSignals(refPosition - SLIDING_WINDOW_LIMIT);
        numProcessedPileups++;
        logProcessedPileups();
    }

    private void initializeSignalsInPos(int pos) {
        snvsPerPos.computeIfAbsent(pos, k -> new ArrayList<>());
        indelsPerPos.computeIfAbsent(pos, k -> new ArrayList<>());
    }

    private void clearUnusedSimpleSignals(int posToFlush) {
        snvsPerPos.remove(posToFlush);
        indelsPerPos.remove(posToFlush);
    }

    /**
     * @param pairedPileup
     *
     */
    public void classifyVariationInPairedPileup(PairedPileup pairedPileup){
        LocusPileUp normalPileup = pairedPileup.getNormalPileup();
        LocusPileUp tumorPileup = pairedPileup.getTumorPileup();
        classifyPileupVariation(normalPileup, true);
        classifyPileupVariation(tumorPileup, false);
    }

    private void classifyPileupVariation(LocusPileUp pileup, boolean isNormalDataset) {
        // In case single normal pileup is queried alone, to guarantee that null tumor pileups are not accessed
        if(pileup==null){
            return;
        }
        //Get only the new reads that appear on this pileup, or those that vary from the reference at this pileup
        List<PileupRead> pileupReads = pileup.claimReadsOnPileup();
        //Check for potential indels or other complex signals
        for (PileupRead pileupRead : pileupReads) {
            if(pileupRead.isNew()){
                //Avoid queuing reads that fall in the region offset, but are not part of the pileup region
                if(overlap(pileupRead, genomicRegion)){
                    AnonymizedRead anonymizedRead = new ShortAnonymizedReadAlignment(pileupRead.getRead(), isNormalDataset);
                    anonymizedRead.setReferenceContigSequence(refSequence);
                    anonymizedReadQueue.offer(anonymizedRead);
                    anonymizedReadCache.put(anonymizedRead.getReadAlignmentId(), anonymizedRead);
                }
                long startdiscoverIndelsAndSignalsFromCIGAR = System.currentTimeMillis();
                discoverIndelsAndComplexSignalsFromCIGAR(pileupRead, isNormalDataset);
                long enddiscoverIndelsAndSignalsFromCIGAR = System.currentTimeMillis();
                METHOD_TIME_MAP.compute("discoverIndelsAndComplexSignalsFromCIGAR", (k, v) -> v == null ?
                        enddiscoverIndelsAndSignalsFromCIGAR - startdiscoverIndelsAndSignalsFromCIGAR :
                        v + enddiscoverIndelsAndSignalsFromCIGAR - startdiscoverIndelsAndSignalsFromCIGAR);
            }
            //Check for potential SNVs
            if(pileupRead.differsFromReferenceAtPileup()){
                long startdiscoverSNVsFromRead = System.currentTimeMillis();
                discoverSNVs(pileupRead, isNormalDataset);
                long enddiscoverSNVsFromRead = System.currentTimeMillis();
                METHOD_TIME_MAP.compute("discoverSNVsFromRead", (k, v) -> v == null ?
                        enddiscoverSNVsFromRead - startdiscoverSNVsFromRead :
                        v + enddiscoverSNVsFromRead - startdiscoverSNVsFromRead);
            }
        }
    }

    public void discoverIndelsAndComplexSignalsFromCIGAR(PileupRead pileupRead, boolean isNormalDataset) {
        List<CigarElement> cigarElems = pileupRead.getCigar().getCigarElements();
        String readAlnId = pileupRead.getReadAlignmentId();
        int initRefPos = pileupRead.getStart();
        int cigarPos = 0;
        int readPos = 0;
        String sequenceName = pileupRead.getSequenceName();
        byte[] sequenceBases = pileupRead.getReadBases();
        for (int i = 0; i < cigarElems.size(); i++){
            CigarElement cigarElement = cigarElems.get(i);
            CigarOperator op = cigarElement.getOperator();
            if (op.isIndel()){
                int currentRefPos = initRefPos + cigarPos-1;
                int inReadPos = pileupRead.getRead().getReadPositionAtReferencePosition(currentRefPos);
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
                byte[] refAllele = Arrays.copyOfRange(refSequence, currentRefPos-1, inRefend);
                altAllele[0] = refAllele[0];
                if (CigarOperator.I.equals(op)){
                    System.arraycopy(sequenceBases, inReadPos, altAllele, 1, altAllele.length - 1);
                }
                PairCalledVariation calledVar = new PairCalledVariation(sequenceName, currentRefPos, vcfStdEnd, indelType, length,
                        altAllele, refAllele);
                List<PairCalledVariation> variationInPos = indelsPerPos.computeIfAbsent(currentRefPos, v -> new ArrayList<>());
                int indexSearch = variationInPos.indexOf(calledVar);
                boolean variationExists = indexSearch != -1;
                if (variationExists) calledVar = variationInPos.get(indexSearch);
                // Saves the CIGAR index of the INDEL signal
                calledVar.addSupportingRead(readAlnId, i);
                processSomaticType(variationInPos, calledVar, variationExists, isNormalDataset);
            }
            if(op.isClipping()){
                int currentRefPos = cigarPos == 0 ? initRefPos : initRefPos + cigarPos-1;
                int length = cigarElement.getLength();
                if(CigarOperator.S == op){
                    //Temp. solution: Avoid soft-clipping signals that fall outside the beginning of the reference sequence,
                    //  or from reads that overlap with the last position
                    if(initRefPos-length-1 >= 0 && !overlap(pileupRead.getStart(), pileupRead.getEnd()+length+1, refSequence.length-1)){
                        Signal calledSignal = new Signal(sequenceName, currentRefPos, readAlnId, i, length,
                                Signal.Source.SOFT_CLIP);
                        calledSignal.setIsFromNormalDataset(isNormalDataset);
                        complexSignals.add(calledSignal);
                        if (anonymizedReadCache.containsKey(readAlnId)) onHoldReads.add(readAlnId);
                    }
                }
            }
            if(op.consumesReferenceBases()){
                cigarPos += cigarElement.getLength();
            }
            if(op.consumesReadBases()){
                readPos += cigarElement.getLength();
            }
        }
    }

    private void discoverSNVs(PileupRead pileupRead, boolean isNormalDataset) {
        int readBase = pileupRead.getBaseAtPileup();
        int referenceBase = pileupRead.getReferenceBase();
        String sequenceName = pileupRead.getSequenceName();
        int refPosition = pileupRead.getLocation();
        String readAlnId = pileupRead.getReadAlignmentId();
        int inReadPosition = pileupRead.getReadPosition();
        byte[] altAllele = new byte[1];
        altAllele[0] = (byte) readBase;
        byte[] refAllele = new byte[1];
        refAllele[0] = (byte) referenceBase;
        PairCalledVariation calledVar = new PairCalledVariation(sequenceName, refPosition, refPosition, VariantType.SNV, 1,
                altAllele, refAllele);
        List<PairCalledVariation> variationInPos = snvsPerPos.computeIfAbsent(refPosition, v -> new ArrayList<>());
        int indexSearch = variationInPos.indexOf(calledVar);
        boolean variationExists = indexSearch != -1;
        if (variationExists) calledVar = variationInPos.get(indexSearch);
        calledVar.addSupportingRead(readAlnId, inReadPosition);
        processSomaticType(variationInPos, calledVar, variationExists, isNormalDataset);
    }

    private void processSomaticType(List<PairCalledVariation> variationInPos, PairCalledVariation calledVar,
                                    boolean variationExists, boolean isNormalDataset) {
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

    private void processSimpleSignals(int refPos, Map<Integer, List<PairCalledVariation>> variationTypePerPos) {
        List<PairCalledVariation> variationInPos = variationTypePerPos.get(refPos);
        for (PairCalledVariation var : variationInPos){
            if(isPotentialGermline(var)) {
                Map<String, Integer> supportingReads = var.getSupportingReadPositions();
                for (Map.Entry<String, Integer> entry : supportingReads.entrySet()){
                    String readAlnId = entry.getKey();
                    AnonymizedRead anonymizedRead = anonymizedReadCache.get(readAlnId);
                    //The anonymizedRead is null if it comes from the offset before the first pileup position,
                    //it is used for classification but is left to be returned by other thread
                    if(anonymizedRead != null) anonymizedRead.addSignalToAnonymize(new Signal(readAlnId, var));
                }
            }
        }
    }

    private boolean isPotentialGermline(PairCalledVariation variation) {
        boolean isValidatedSomatic = false;
        // Anonymize only potential germlines if seen in both datasets, at least once in each, or more than once if only found in the normal tissue mappings
        boolean isPotentialGermline = SomaticVariationType.TUMORAL_NORMAL_VARIANT.equals(variation.getSomaticVariationType()) ||
                SomaticVariationType.NORMAL_ONLY_VARIANT.equals(variation.getSomaticVariationType());
        //Avoid anonymizing somatic variants recorded in the VCF file
        if (!somaticVariantsToKeep.isEmpty()){
            if(somaticVariantsToKeep.containsKey(variation.getSeqName())){
                Map<Integer, PairCalledVariation> validatedSomaticsAtSeq = somaticVariantsToKeep.get(variation.getSeqName());
                if(validatedSomaticsAtSeq.containsKey(variation.getPos())){
                    PairCalledVariation validatedSomaticAtPos = validatedSomaticsAtSeq.get(variation.getPos());
                    if(validatedSomaticAtPos.equals(variation)) isValidatedSomatic = true;
                }
            }
        }
        if(isValidatedSomatic){
            return isPotentialGermline && anonymizePotentialLeaksInVCFSomatics;
        }
        return isPotentialGermline;
    }

    /**
     * Classify complex signals into potential germline variations by virtually inferring either a complete graph of
     * signals from the normal dataset, or a bipartite graph from the normal-tumor pair
     *
     */
    private void processComplexSignals() {
//        List<List<Integer>> adjacencyGraph = new ArrayList<>();
//        Collections.fill(adjacencyGraph, new ArrayList<>());
        int n = complexSignals.size();
        List<Signal> currentPartition = new ArrayList<>();
        for(int i = 0; i < n-1; i++){
            Signal currentSignal = complexSignals.get(i);
            currentPartition.add(currentSignal);
            Signal nextSignal = complexSignals.get(i+1);
            int locationDistance = Math.abs(nextSignal.getLocation() - currentSignal.getLocation());
            boolean lastSignalUnreachable = locationDistance > MAX_LOCATION_DISTANCE_THRESHOLD;
            if(lastSignalUnreachable || currentPartition.size() >= SIGNAL_PER_REGION_LIMIT || (i==n-2)){
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
            onHoldReads.remove(firstSignal.getReadAlnName());
        }
    }

    private void classifyPGcomplexSignal(Signal signal, int pos, boolean[] isClassifiedPG) {
        if(!isClassifiedPG[pos]) {
            String readAlnId = signal.getReadAlnName();
            AnonymizedRead anonymizedRead = anonymizedReadCache.get(readAlnId);
            //The anonymizedRead is null if it comes from the offset before the first pileup position,
            //it is used for classification but is left to be returned by other thread
            if(anonymizedRead != null) anonymizedRead.addSignalToAnonymize(signal);
            isClassifiedPG[pos] = true;
            //DEBUG
            //System.out.println("& " + signal.toString());
            //DEBUG
        }
    }

    private void logProcessedPileups() {
        if (numProcessedPileups % 1_000_000 == 0) LOGGER.info("GenomicRegion["+ genomicRegion.toString()+"]: "+ "Processed " + numProcessedPileups + " pileup positions");
    }

    public void setDiffuseIndelCalls(boolean diffuseIndelCalls) {
        this.diffuseIndelCalls = diffuseIndelCalls;
    }

    @Override
    public Iterator<AnonymizedRead> iterator() {
        return new Iterator<AnonymizedRead>() {

            final Iterator<PairedPileup> pairedPileupIterator = pairPileupReader.iterator();
            AnonymizedRead next = getNextAnonymizedRead();

            @Override
            public boolean hasNext() {
                return next != null;
            }

            @Override
            public AnonymizedRead next() {
                if(next == null) throw new NoSuchElementException();
                AnonymizedRead current = next;
                next = getNextAnonymizedRead();
                return current;
            }

            private AnonymizedRead getNextAnonymizedRead() {
                //Advance the pileup or get the next fully processed read
                AnonymizedRead answer;
                while(pairedPileupIterator.hasNext()){
                    answer = anonymizedReadQueue.peek();
                    if(answer != null && answer.getEnd() < currentPileupPosition){
                        answer = anonymizedReadQueue.remove();
                        if (!onHoldReads.contains(answer.getReadAlignmentId())){
                            anonymizedReadCache.remove(answer.getReadAlignmentId());
                            return answer;
                        }
                        else{
                            anonymizedReadQueue.offer(answer);
                        }
                    }
                    processNextPairedPileup(pairedPileupIterator.next());
                }
                //Get the next remaining reads from the queue
                if(!anonymizedReadQueue.isEmpty()){
                    answer = anonymizedReadQueue.poll();
                    anonymizedReadCache.remove(answer.getReadAlignmentId());
                    return answer;
                }
                return null;
            }
        };
    }

    @Override
    public void forEach(Consumer<? super AnonymizedRead> action) {
        Iterable.super.forEach(action);
    }

    @Override
    public Spliterator<AnonymizedRead> spliterator() {
        return Iterable.super.spliterator();
    }

    @Override
    public void close() throws IOException {
        pairPileupReader.close();
    }
}
