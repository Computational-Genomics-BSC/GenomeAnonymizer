package analysis;
import genomicelements.*;
import genomicelements.PairCalledVariation.SomaticVariationType;
import genomicelements.PairCalledVariation.VariantType;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import io.SamplePairReadAlignmentReader;
import utils.MapCacheFIFO;
import utils.Operations;

import java.io.Closeable;
import java.io.IOException;
import java.util.*;
import java.util.function.Consumer;
import java.util.logging.Logger;

import static analysis.GenomeAnonymizer.DEFAULT_MIN_MAPPING_QUALITY;
import static utils.Operations.isContained;
import static utils.Operations.overlap;


/**
 * Class used to classify all variation from a paired normal-tumor sample from each pileup position, and provide
 * variated and not read alignments
 * @author Nicolas Gaitan
 */

public class AnonymizedReadAlignmentProvider implements Iterable<AnonymizedRead>, Closeable {

    private static final Logger LOGGER = Logger.getLogger(AnonymizedReadAlignmentProvider.class.getName());

    public static final int SHORT_READ_LEFT_PILEUP_REGION_EXTENSION = 500;
    public static final int SLIDING_WINDOW_LIMIT = 200;
    public static final int MAX_LOCATION_DISTANCE_THRESHOLD = 400;

//    private static final int INDEL_SIGNAL_PER_REGION_LIMIT = 100;
    public static final int MAX_SIGNAL_PER_REGION_LIMIT = 5000;

    public static final int INDEL_SIGNAL_PER_REGION_LIMIT = 500;
    public static final int COMPLEX_SIGNAL_PER_REGION_LIMIT = 1000;

    // Assuming a maximum position distance of 10, and 10 of length difference
    public static final double INDEL_SIGNAL_THRESHOLD = 14.14;
    // Assuming a maximum position distance of 15, and 25 of length difference
    public static final double COMPLEX_SIGNAL_THRESHOLD = 29.15;

    SamplePairReadAlignmentReader pairPileupReader;
    private int currentPileupPosition = 0;

    private Queue<AnonymizedRead> anonymizedReadQueue;
    private MapCacheFIFO<String, AnonymizedRead> anonymizedReadCache;
    private Map<String, Integer> onHoldReads;

    private Map<Integer, List<PairCalledVariation>> snvsPerPos;
    private Map<Integer, List<PairCalledVariation>> indelsPerPos;
    private List<Signal> signals;
    private Set<Integer> partiallyUncoveredPositions;

    private Set<String> readsToExclude;
    private int insertSizeMinThreshold;
    private int insertSizeMaxThreshold;

    private byte[] refSequence;
    private GenomicRegion genomicRegion;

    private int minMappingQuality = DEFAULT_MIN_MAPPING_QUALITY;

    private int numProcessedPileups = 0;

    //DEBUG
    private GenomicRegion extendedRegion;
    public Map<String, Long> METHOD_TIME_MAP = new HashMap<>();
    //DEBUG

    public AnonymizedReadAlignmentProvider(){
        anonymizedReadQueue = new LinkedList<>();
        anonymizedReadCache = new MapCacheFIFO<>(10_000_000);
        onHoldReads = new HashMap<>();
        readsToExclude = new HashSet<>();
        snvsPerPos = new HashMap<>();
        indelsPerPos = new HashMap<>();
        signals = new ArrayList<>();
        partiallyUncoveredPositions = new HashSet<>();
        refSequence = new byte[0];
    }

    public void setReadsToExclude(Set<String> readsToExclude){
        this.readsToExclude = readsToExclude;
    }

    public void setMinMappingQuality(int minMappingQuality) {
        this.minMappingQuality = minMappingQuality;
    }

    public void setRefSequence(byte[] refSequence) {
        this.refSequence = refSequence;
    }

    public void setInsertSizeMinThreshold(int insertSizeMinThreshold) {
        this.insertSizeMinThreshold = insertSizeMinThreshold;
    }

    public void setInsertSizeMaxThreshold(int insertSizeMaxThreshold) {
        this.insertSizeMaxThreshold = insertSizeMaxThreshold;
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
        pairPileupReader.setIncludeDuplicates(true);
        pairPileupReader.setReadsToExclude(readsToExclude);
        //DEBUG
        extendedRegion = leftExtendedRegion;
        //DEBUG
    }

    public void processNextPairedPileup(PairedPileup pileup, boolean hasNext){
        int refPosition = pileup.getLocation();
        currentPileupPosition = refPosition;
        initializeSignalsInPos(refPosition);
        long startclassifyVariationInPairedPileup = System.currentTimeMillis();
        classifyVariationInPairedPileup(pileup);
        long endclassifyVariationInPairedPileup = System.currentTimeMillis();
        METHOD_TIME_MAP.compute("classifyVariationInPairedPileup", (k,v) -> v == null ?
                endclassifyVariationInPairedPileup-startclassifyVariationInPairedPileup :
                v + endclassifyVariationInPairedPileup-startclassifyVariationInPairedPileup);
        long startprocessSimpleSignals = System.currentTimeMillis();
        processSNVSignals(refPosition, snvsPerPos);
        long endprocessSimpleSignals = System.currentTimeMillis();
        METHOD_TIME_MAP.compute("processSimpleSNVSignals",  (k,v) -> v == null ?
                endprocessSimpleSignals-startprocessSimpleSignals :
                v + endprocessSimpleSignals-startprocessSimpleSignals);
        if ( !hasNext || signals.size() >= MAX_SIGNAL_PER_REGION_LIMIT ) {
            long startprocessComplexSignals = System.currentTimeMillis();
            processSignals(signals);
            long endprocessComplexSignals = System.currentTimeMillis();
            METHOD_TIME_MAP.compute("processComplexSignals",  (k,v) -> v == null ?
                    endprocessComplexSignals-startprocessComplexSignals :
                    v + endprocessComplexSignals-startprocessComplexSignals);
            signals = new ArrayList<>();
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
        if(normalPileup == null || tumorPileup == null){
            partiallyUncoveredPositions.add(pairedPileup.getLocation());
        }
        classifyPileupVariation(normalPileup, true);
        classifyPileupVariation(tumorPileup, false);
    }

    private void classifyPileupVariation(LocusPileUp pileup, boolean isNormalDataset) {
        // Guarantee that null pileups are not accessed
        if(pileup==null){
            return;
        }
        //Get only the new reads that appear on this pileup, or those that vary from the reference at this pileup
        List<PileupRead> pileupReads = pileup.claimReadsOnPileup();
        //Check for potential indels or other complex signals
        for (PileupRead pileupRead : pileupReads) {
            if(pileupRead.isNew()){
                //Avoid queuing reads that fall in the region offset, but are not part of the pileup region
                if(belongsToRegion(pileupRead, genomicRegion)){
                    AnonymizedRead anonymizedRead = new ShortAnonymizedReadAlignment(pileupRead.getRead(), isNormalDataset);
                    anonymizedRead.setReferenceContigSequence(refSequence);
                    anonymizedReadQueue.offer(anonymizedRead);
                    anonymizedReadCache.put(anonymizedRead.getReadAlignmentId(), anonymizedRead);
                }
                long startdiscoverIndelsAndSignalsFromCIGAR = System.currentTimeMillis();
                discoverSignalsFromCIGAR(pileupRead, isNormalDataset);
                long enddiscoverIndelsAndSignalsFromCIGAR = System.currentTimeMillis();
                METHOD_TIME_MAP.compute("discoverIndelsAndComplexSignalsFromCIGAR", (k, v) -> v == null ?
                        enddiscoverIndelsAndSignalsFromCIGAR - startdiscoverIndelsAndSignalsFromCIGAR :
                        v + enddiscoverIndelsAndSignalsFromCIGAR - startdiscoverIndelsAndSignalsFromCIGAR);
                // One signal per mate
                discoverSignalsFromMates(pileupRead);
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

    private boolean belongsToRegion(PileupRead pileupRead, GenomicRegion genomicRegion) {
        if (isContained(pileupRead, genomicRegion)) return true;
        if (!overlap(pileupRead, genomicRegion)) return false;
        int distanceToRegionStart = Math.abs(pileupRead.getStart() - genomicRegion.getStart());
        int distanceToRegionEnd = Math.abs(pileupRead.getStart() - genomicRegion.getEnd());
        int closestDistance = Math.min(distanceToRegionStart, distanceToRegionEnd);
        return closestDistance == distanceToRegionStart;
    }

    public void discoverSignalsFromCIGAR(PileupRead pileupRead, boolean isNormalDataset) {
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
                Signal calledSignal = new Signal(readAlnId, calledVar);
                calledSignal.setIsFromNormalDataset(isNormalDataset);
                signals.add(calledSignal);
                if (anonymizedReadCache.containsKey(readAlnId)) {
                    onHoldReads.compute(readAlnId, (k, v) -> v == null ? 1 : v + 1);
                }
            }
            if(op.isClipping()){
                int currentRefPos = cigarPos == 0 ? initRefPos : initRefPos + cigarPos-1;
                int length = cigarElement.getLength();
                if(CigarOperator.S == op){
                    //Avoid soft-clipping signals that fall outside the beginning of the reference sequence,
                    // or from reads that overlap with the last position
                    if(initRefPos-length-1 >= 0 && !overlap(pileupRead.getStart(), pileupRead.getEnd()+length+1, refSequence.length-1)){
                        Signal calledSignal = new Signal(sequenceName, currentRefPos, readAlnId, i, length,
                                Signal.Source.SOFT_CLIP);
                        calledSignal.setIsFromNormalDataset(isNormalDataset);
                        signals.add(calledSignal);
                        if (anonymizedReadCache.containsKey(readAlnId)) {
                            onHoldReads.compute(readAlnId, (k, v) -> v == null ? 1 : v + 1);
                        }
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

    public void discoverSignalsFromMates(PileupRead pileupRead) {
        SAMRecord samRecord = pileupRead.getRead();
        // Check if reads are in different chromosomes
        if (!samRecord.getReferenceIndex().equals(samRecord.getMateReferenceIndex())) {
            Signal calledSignal = new Signal(samRecord.getContig(), samRecord.getAlignmentStart(), pileupRead.getReadAlignmentId(), 0, 0, Signal.Source.CHROM_CHANGE);
            signals.add(calledSignal);
            return;
        }
        // Check signal strands: FF, RF and RR
        // Check if this is the first or second pair (assume both are mapped)
        boolean firstRead = samRecord.getAlignmentStart() <= samRecord.getMateAlignmentStart();
        boolean firstForward;
        boolean secondForward;
        if (firstRead) {
            firstForward = !samRecord.getReadNegativeStrandFlag();
            secondForward = !samRecord.getMateNegativeStrandFlag();
        }
        else {
            firstForward = !samRecord.getMateNegativeStrandFlag();
            secondForward = !samRecord.getReadNegativeStrandFlag();
        }
        int insertSize = Math.abs(samRecord.getInferredInsertSize());
        if ((firstForward && secondForward) || (!firstForward && !secondForward) || (!firstForward && secondForward)) {
            Signal calledSignal = new Signal(samRecord.getContig(), samRecord.getAlignmentStart(), pileupRead.getReadAlignmentId(), 0, insertSize, Signal.Source.STRAND_ORIENTATION);
            signals.add(calledSignal);
            return;
        }
        // Check insert size
        if (insertSize < insertSizeMinThreshold || insertSize > insertSizeMaxThreshold) {
            Signal calledSignal = new Signal(samRecord.getContig(), samRecord.getAlignmentStart(), pileupRead.getReadAlignmentId(), 0, insertSize, Signal.Source.INSERT_SIZE);
            signals.add(calledSignal);
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

    private void processSNVSignals(int refPos, Map<Integer, List<PairCalledVariation>> variationTypePerPos) {
        List<PairCalledVariation> variationInPos = variationTypePerPos.get(refPos);
        for (PairCalledVariation var : variationInPos){
            if( var.isGermline() || overlapsUncoveredPosition(var) ) {
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

    /**
     * Classify complex signals into potential germline variations by virtually inferring either a complete graph of
     * signals from the normal dataset, or a bipartite graph from the normal-tumor pair
     *
     */
    private void processSignals(List<Signal> signals) {
        Map<String, List<Signal>> signalTypeMap = new HashMap<>();
        // For each source of signal, create an empty list
        for(Signal.Source source : Signal.Source.values()){
            signalTypeMap.put(source.name(), new ArrayList<>());
        }
        for(Signal signal : signals){
            // Get the signal list
            List<Signal> signalList = signalTypeMap.get(signal.getSource().name());
            signalList.add(signal);
        }
        processTypeSignals(signalTypeMap.get(Signal.Source.SIMPLE_VARIATION.name()), INDEL_SIGNAL_PER_REGION_LIMIT, INDEL_SIGNAL_THRESHOLD);
        processTypeSignals(signalTypeMap.get(Signal.Source.SOFT_CLIP.name()), COMPLEX_SIGNAL_PER_REGION_LIMIT, COMPLEX_SIGNAL_THRESHOLD);
        processTypeSignals(signalTypeMap.get(Signal.Source.INSERT_SIZE.name()), COMPLEX_SIGNAL_PER_REGION_LIMIT, COMPLEX_SIGNAL_THRESHOLD);
        processTypeSignals(signalTypeMap.get(Signal.Source.STRAND_ORIENTATION.name()), COMPLEX_SIGNAL_PER_REGION_LIMIT, COMPLEX_SIGNAL_THRESHOLD);
        processTypeSignals(signalTypeMap.get(Signal.Source.CHROM_CHANGE.name()), COMPLEX_SIGNAL_PER_REGION_LIMIT, COMPLEX_SIGNAL_THRESHOLD);
        //Compare remaining cigar signals from different concordant types
        List<Signal> remainingSignals = new ArrayList<>();
        List<String> remainingSignalTypes = new ArrayList<>();
        remainingSignalTypes.add(Signal.Source.SIMPLE_VARIATION.name());
        remainingSignalTypes.add(Signal.Source.SOFT_CLIP.name());
        for (String type : remainingSignalTypes){
            for (Signal signal : signalTypeMap.get(type)){
                if( !signal.isGermline() ){
                    remainingSignals.add(signal);
                }
            }
        }
        processTypeSignals(remainingSignals, COMPLEX_SIGNAL_PER_REGION_LIMIT, COMPLEX_SIGNAL_THRESHOLD);
    }

    private void processTypeSignals(List<Signal> typeSignals, int signalPerRegionLimit, double signalTypeThreshold) {
        int n = typeSignals.size();
        List<Signal> currentPartition = new ArrayList<>();
        for(int i = 0; i < n-1; i++){
            Signal currentSignal = typeSignals.get(i);
            currentPartition.add(currentSignal);
            Signal nextSignal = typeSignals.get(i+1);
            int locationDistance = Math.abs(nextSignal.getLocation() - currentSignal.getLocation());
            boolean lastSignalUnreachable = locationDistance > MAX_LOCATION_DISTANCE_THRESHOLD;
            if(lastSignalUnreachable || currentPartition.size() >= signalPerRegionLimit || (i==n-2)){
                if(!lastSignalUnreachable){
                    currentPartition.add(nextSignal);
                    i++;
                }
                processPartition(currentPartition, signalTypeThreshold);
                currentPartition = new ArrayList<>();
            }
        }
    }

    private void processPartition(List<Signal> signals, double signalTypeThreshold) {
        int n = signals.size();
        for (int i = 0; i < n; i++){
            Signal firstSignal = signals.get(i);
            for (int j = i + 1; j < n; j++){
                Signal secondSignal = signals.get(j);
                double signalDistance = Operations.computeTwoDimEuclideanDistance(firstSignal.getLocation(), secondSignal.getLocation(),
                        firstSignal.getLength(), secondSignal.getLength());
                //Classify signals that are close enough to be considered as a germline signal
                if ( signalDistance < signalTypeThreshold &&
                        (firstSignal.isFromNormalDataset() || secondSignal.isFromNormalDataset()) ) {
                    long startclassifyPGcomplexSignal = System.currentTimeMillis();
                    classifyGermlineComplexSignal(firstSignal);
                    classifyGermlineComplexSignal(secondSignal);
                    long endclassifyPGcomplexSignal = System.currentTimeMillis();
                    METHOD_TIME_MAP.compute("classifyPGcomplexSignal", (k, v) -> v == null ?
                            endclassifyPGcomplexSignal - startclassifyPGcomplexSignal :
                            v + endclassifyPGcomplexSignal - startclassifyPGcomplexSignal);
                }
                else{
                    //Classify signals that are not close enough to be considered as a germline signal,
                    // but are undecidable as they are not overlapped by reads on the other dataset,
                    // or SIMPLE_VARIATION that is already classified as germline
                    if(overlapsUncoveredPosition(firstSignal) ||
                            ( firstSignal.getSource() == Signal.Source.SIMPLE_VARIATION && firstSignal.getCalledVariation().isGermline() )){
                        classifyGermlineComplexSignal(firstSignal);
                    }
                    if(overlapsUncoveredPosition(secondSignal) ||
                            ( secondSignal.getSource() == Signal.Source.SIMPLE_VARIATION && secondSignal.getCalledVariation().isGermline() )){
                        classifyGermlineComplexSignal(secondSignal);
                    }
                }
            }
            if(onHoldReads.containsKey(firstSignal.getReadAlnName())){
                int remainingReadSignals = onHoldReads.get(firstSignal.getReadAlnName());
                if(remainingReadSignals == 1){
                    onHoldReads.remove(firstSignal.getReadAlnName());
                }
                else{
                    onHoldReads.put(firstSignal.getReadAlnName(), remainingReadSignals-1);
                }
            }
        }
    }

    private void classifyGermlineComplexSignal(Signal signal) {
        if( !signal.isGermline() ) {
            String readAlnId = signal.getReadAlnName();
            AnonymizedRead anonymizedRead = anonymizedReadCache.get(readAlnId);
            //The anonymizedRead is null if it comes from the offset before the first pileup position,
            //it is used for classification but is left to be returned by other thread
            if(anonymizedRead != null) anonymizedRead.addSignalToAnonymize(signal);
            signal.setIsGermline(true);
        }
    }

    private boolean overlapsUncoveredPosition(GenomicRegion varSignal) {
        for(int i = varSignal.getStart(); i <= varSignal.getEnd(); i++){
            if(partiallyUncoveredPositions.contains(i)) return true;
        }
        return false;
    }

    private void logProcessedPileups() {
        if (numProcessedPileups % 1_000_000 == 0) LOGGER.info("GenomicRegion["+ genomicRegion.toString()+"]: "+ "Processed " + numProcessedPileups + " pileup positions");
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
                        if (!onHoldReads.containsKey(answer.getReadAlignmentId())){
                            anonymizedReadCache.remove(answer.getReadAlignmentId());
                            return answer;
                        }
                        else{
                            anonymizedReadQueue.offer(answer);
                        }
                    }
                    processNextPairedPileup(pairedPileupIterator.next(), pairedPileupIterator.hasNext());
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
