package analysis;
import genomicelements.*;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import io.SamplePairReadAlignmentReader;
import smile.clustering.*;
import smile.clustering.linkage.*;
import smile.graph.AdjacencyList;
import smile.math.distance.Distance;
import static smile.math.MathEx.*;

import utils.MapCacheFIFO;
import utils.Operations;

import java.io.Closeable;
import java.io.IOException;
import java.util.*;
import java.util.function.Consumer;
import java.util.logging.Logger;

import static analysis.GenomeAnonymizer.DEFAULT_MIN_MAPPING_QUALITY;
import static utils.Operations.*;


/**
 * Class used to classify all variation from a paired normal-tumor sample from each pileup position, and provide
 * variated and not read alignments
 * @author Nicolas Gaitan
 */

public class AnonymizedReadAlignmentProvider implements Iterable<AnonymizedRead>, Closeable {

    private static final Logger LOGGER = Logger.getLogger(AnonymizedReadAlignmentProvider.class.getName());

    public static final int SHORT_READ_LEFT_PILEUP_REGION_EXTENSION = 500;
    public static final int MAX_LOCATION_DISTANCE_THRESHOLD = 400;

    public static final int MAX_SIGNAL_PER_REGION_LIMIT = 5000;

    public static final int SIGNAL_PER_PARTITION_LIMIT = 1500;

    // Assuming a maximum position distance of 5, and 5 of length difference
    public static final double INDEL_SIGNAL_THRESHOLD = 7.07;
    // Assuming a maximum position distance of 3, and 3 of length difference
    public static final double INDEL_SOFTCLIP_SIGNAL_THRESHOLD = 4.24;
    // Assuming a maximum position distance of 15, and 25 of length difference
    public static final double SOFTCLIP_SIGNAL_THRESHOLD = 29.15;
    // Assuming a maximum position distance of 50, and 250 of length difference
    public static final double INSERT_SIZE_SIGNAL_THRESHOLD = 250;
    // Assuming a maximum position distance of 150
    public static final double STRAND_ORIENTATION_SIGNAL_THRESHOLD = 150;


    SamplePairReadAlignmentReader pairPileupReader;
    private int currentPileupPosition = 0;

    private Queue<AnonymizedRead> anonymizedReadQueue;
    private MapCacheFIFO<String, AnonymizedRead> anonymizedReadCache;
    private Map<String, Integer> onHoldReads;

    private List<Signal> snvSignals;
    private SignalCollection signals;
    private Set<Integer> partiallyUncoveredPositions;

    private Set<String> readsToExclude;
    private int insertSizeMinThreshold;
    private int insertSizeMaxThreshold;

    private byte[] refSequence;
    private GenomicRegion genomicRegion;

    private int minMappingQuality = DEFAULT_MIN_MAPPING_QUALITY;
    private boolean includeDuplicates = false;

    private int numProcessedPileups = 0;

    private String sampleType = GenomeAnonymizer.SAMPLE_TYPE_WGS;
    private int minDepthForVAFCorrection = GenomeAnonymizer.DEFAULT_MIN_DEPTH_FOR_VAF_CORRECTION;

    // Map of chromosome -> (sorted position -> original VAF data) for potential somatic sites
    private Map<String, TreeMap<Integer, SomaticSiteVAF>> somaticSiteVAFs = new HashMap<>();

    //DEBUG
    public Map<String, Long> METHOD_TIME_MAP = new HashMap<>();
    //DEBUG

    public AnonymizedReadAlignmentProvider(){
        anonymizedReadQueue = new LinkedList<>();
        anonymizedReadCache = new MapCacheFIFO<>(10_000_000);
        onHoldReads = new HashMap<>();
        readsToExclude = new HashSet<>();
        snvSignals = new ArrayList<>();
        signals = new SignalCollection();
        partiallyUncoveredPositions = new HashSet<>();
        refSequence = new byte[0];
    }

    public void setReadsToExclude(Set<String> readsToExclude){
        this.readsToExclude = readsToExclude;
    }

    public void setMinMappingQuality(int minMappingQuality) {
        this.minMappingQuality = minMappingQuality;
    }

    public void setIncludeDuplicates(boolean includeDuplicates) {
        this.includeDuplicates = includeDuplicates;
    }

    public void setSampleType(String sampleType) {
        this.sampleType = sampleType;
    }

    public void setMinDepthForVAFCorrection(int minDepthForVAFCorrection) {
        this.minDepthForVAFCorrection = minDepthForVAFCorrection;
    }

    public Map<String, TreeMap<Integer, SomaticSiteVAF>> getSomaticSiteVAFs() {
        return somaticSiteVAFs;
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
        pairPileupReader.setIncludeDuplicates(includeDuplicates);
        pairPileupReader.setReadsToExclude(readsToExclude);
        pairPileupReader.setMinimumMappingQuality(minMappingQuality);
    }

    public void processNextPairedPileup(PairedPileup pileup, boolean hasNext){
        currentPileupPosition = pileup.getLocation();
        long startclassifyVariationInPairedPileup = System.currentTimeMillis();
        classifyVariationInPairedPileup(pileup);
        long endclassifyVariationInPairedPileup = System.currentTimeMillis();
        METHOD_TIME_MAP.compute("classifyVariationInPairedPileup", (k,v) -> v == null ?
                endclassifyVariationInPairedPileup-startclassifyVariationInPairedPileup :
                v + endclassifyVariationInPairedPileup-startclassifyVariationInPairedPileup);
        long startprocessSimpleSignals = System.currentTimeMillis();
        if (GenomeAnonymizer.SAMPLE_TYPE_GENE_PANEL.equals(sampleType)) {
            int normalDepth = pileup.getNormalPileup() != null ? pileup.getNormalPileup().getPropperSize() : 0;
            int tumorDepth = pileup.getTumorPileup() != null ? pileup.getTumorPileup().getPropperSize() : 0;
            processGenePanelSNVSignals(normalDepth, tumorDepth);
        }
        else {
            // WGS mode - use original stringent logic
            processWGSSNVSignals();
        }
        snvSignals.clear();
        long endprocessSimpleSignals = System.currentTimeMillis();
        METHOD_TIME_MAP.compute("processSimpleSNVSignals",  (k,v) -> v == null ?
                endprocessSimpleSignals-startprocessSimpleSignals :
                v + endprocessSimpleSignals-startprocessSimpleSignals);
        if ( !hasNext || signals.size() >= MAX_SIGNAL_PER_REGION_LIMIT ) {
            long startprocessComplexSignals = System.currentTimeMillis();
            processSignals();
            long endprocessComplexSignals = System.currentTimeMillis();
            METHOD_TIME_MAP.compute("processComplexSignals",  (k,v) -> v == null ?
                    endprocessComplexSignals-startprocessComplexSignals :
                    v + endprocessComplexSignals-startprocessComplexSignals);
            signals.clear();
        }
        numProcessedPileups++;
        logProcessedPileups();
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
                discoverSignalsFromMates(pileupRead, isNormalDataset);
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
                Signal.IndelSignalType indelType;
                int vcfStdEnd;
                int inRefend;
                int inReadEnd;
                byte[] altAllele;
                if (CigarOperator.I == op){
                    indelType = Signal.IndelSignalType.INSERTION;
                    inRefend = currentRefPos;// + 1;
                    vcfStdEnd = inRefend + 1;
                    inReadEnd = inReadPos + length + 1;
                    altAllele = new byte[1+length];
                }
                else{
                    indelType = Signal.IndelSignalType.DELETION;
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
                Signal signal = new Signal(sequenceName, currentRefPos, readAlnId, i, length, Signal.Source.INDEL);
                signal.setIndelSignalType(indelType);
                signal.setAltAllele(altAllele);
                signal.setRefAllele(refAllele);
                signal.setIsFromNormalDataset(isNormalDataset);
                signals.add(signal);
                putReadOnHold(readAlnId);
            }
            if(op.isClipping()){
                int currentRefPos = cigarPos == 0 ? initRefPos : initRefPos + cigarPos-1;
                int length = cigarElement.getLength();
                if(CigarOperator.S == op){
                    //Avoid soft-clip signals that fall outside the beginning of the reference sequence,
                    // or from reads that overlap with the last position
                    if(initRefPos - pileupRead.getReadLength() - 1 >= 0 && !overlap(pileupRead.getStart(), pileupRead.getEnd()+pileupRead.getReadLength()+1, refSequence.length-1)){
                        Signal calledSignal = new Signal(sequenceName, currentRefPos, readAlnId, i, length, Signal.Source.SOFT_CLIP);
                        calledSignal.setIsFromNormalDataset(isNormalDataset);
                        signals.add(calledSignal);
                        putReadOnHold(readAlnId);
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

    public void discoverSignalsFromMates(PileupRead pileupRead, boolean isNormalDataset) {
        String readAlnId = pileupRead.getReadAlignmentId();
        SAMRecord samRecord = pileupRead.getRead();
        // Check if the read maps before the mate
        int cmp = compare(samRecord.getReferenceIndex(), samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(),
                samRecord.getMateReferenceIndex(), samRecord.getMateAlignmentStart(),
                samRecord.getMateAlignmentStart() + samRecord.getReadLength() - 1);
        boolean thisPairMapsFirst = cmp <= 0;
        // Check if the read is the pair mapped first
        if (samRecord.isSecondaryOrSupplementary() || !thisPairMapsFirst) {
            return;
        }
        // Check insert size
        int insertSize = Math.abs(samRecord.getInferredInsertSize());
        // Check if reads are in different chromosomes
        if (!samRecord.getReferenceIndex().equals(samRecord.getMateReferenceIndex())) {
            Signal calledSignal = new Signal(samRecord.getContig(), samRecord.getAlignmentStart(), pileupRead.getReadAlignmentId(), 0, 0, Signal.Source.CHROM_CHANGE);
            calledSignal.setIsFromNormalDataset(isNormalDataset);
            signals.add(calledSignal);
            putReadOnHold(readAlnId);
            return;
        }
        if (insertSize < insertSizeMinThreshold || insertSize > insertSizeMaxThreshold) {
            Signal calledSignal = new Signal(samRecord.getContig(), samRecord.getAlignmentStart(), pileupRead.getReadAlignmentId(), 0, insertSize, Signal.Source.INSERT_SIZE);
            calledSignal.setIsFromNormalDataset(isNormalDataset);
            signals.add(calledSignal);
            putReadOnHold(readAlnId);
            return;
        }
        // Check if this is the first or second pair (assume both are mapped)
        boolean firstRead = samRecord.getAlignmentStart() <= samRecord.getMateAlignmentStart();
        // Check signal strands: FF, RF and RR
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
        if ((firstForward && secondForward) || (!firstForward && !secondForward) || (!firstForward && secondForward)) {
            Signal calledSignal = new Signal(samRecord.getContig(), samRecord.getAlignmentStart(), pileupRead.getReadAlignmentId(), 0, insertSize, Signal.Source.STRAND_ORIENTATION);
            calledSignal.setIsFromNormalDataset(isNormalDataset);
            signals.add(calledSignal);
            putReadOnHold(readAlnId);
        }
    }

    private void putReadOnHold(String readAlnId) {
        if (anonymizedReadCache.containsKey(readAlnId)) {
            onHoldReads.compute(readAlnId, (k, v) -> v == null ? 1 : v + 1);
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
        Signal snvSignal = new Signal(sequenceName, refPosition, readAlnId, inReadPosition, 1, Signal.Source.SNV);
        snvSignal.setIsFromNormalDataset(isNormalDataset);
        snvSignal.setAltAllele(altAllele);
        snvSignal.setRefAllele(refAllele);
        if(includeDuplicates && pileupRead.getRead().getDuplicateReadFlag()) snvSignal.setComesFromDuplicate(true);
        snvSignal.setBaseQuality(pileupRead.getRead().getBaseQualities()[inReadPosition]);
        snvSignals.add(snvSignal);
    }

    private void processWGSSNVSignals() {
        boolean[] seen = new boolean[256];
        boolean[] isGermline = new boolean[256];
        for (Signal var : snvSignals) {
            seen[var.getAltAllele()[0]] = true;
            if ( (seen[var.getAltAllele()[0]] && var.isFromNormalDataset())
                    || overlapsUncoveredPosition(var) ) {
                isGermline[var.getAltAllele()[0]] = true;
            }
        }
        for (Signal var : snvSignals){
            if(isGermline[var.getAltAllele()[0]]) {
                String readAlnId = var.getReadAlnName();
                AnonymizedRead anonymizedRead = anonymizedReadCache.get(readAlnId);
                //The anonymizedRead is null if it comes from the offset before the first pileup position,
                //it is used for classification but is left to be returned by other thread
                if(anonymizedRead != null) anonymizedRead.addSignalToAnonymize(var);
            }
        }
    }

    /**
     * Process SNV signals for gene panel mode with less stringent germline determination.
     * Only marks signals as germline if they appear in more normal reads than the calculated threshold.
     */
    private void processGenePanelSNVSignals(int normalDepth, int tumorDepth) {
        int readThreshold = calculateGenePanelThreshold(normalDepth);

        // Count how many normal/tumor reads have each alt allele and collect their base qualities
        // Using int arrays indexed by byte values for more efficient access than a HashMap
        int[] altAlleleNormalReadCount = new int[256];
        int[] altAlleleTumoralReadCount = new int[256];
        List<Byte>[] normalAltQualities = new List[256];
        List<Byte>[] tumorAltQualities = new List[256];
        for (Signal var : snvSignals) {
            byte altAllele = var.getAltAllele()[0];
            if(var.comesFromDuplicateRead()) continue;
            int alleleIdx = altAllele & 0xFF;
            if (var.isFromNormalDataset()) {
                altAlleleNormalReadCount[altAllele]++;
                if (normalAltQualities[alleleIdx] == null) normalAltQualities[alleleIdx] = new ArrayList<>();
                normalAltQualities[alleleIdx].add(var.getBaseQuality());
            }
            else {
                altAlleleTumoralReadCount[altAllele]++;
                if (tumorAltQualities[alleleIdx] == null) tumorAltQualities[alleleIdx] = new ArrayList<>();
                tumorAltQualities[alleleIdx].add(var.getBaseQuality());
            }
        }

        // Mark as germline only if the alt allele appears in enough normal reads to exceed the threshold
        // This prevents marking error-generated variation as germline (thus allowing somatics in output)
        // In gene panel mode, we're less stringent: variants present in few reads are allowed through
        // Track which alleles have already been captured at this position to avoid duplicates
        boolean[] alleleCaptured = new boolean[256];
        double[] normalALTReadFractions = new double[256];
        double[] tumorALTReadFractions = new double[256];
        List<Signal> candidateSomaticSignals = new ArrayList<>();
        for (Signal var : snvSignals) {
            byte altAllele = var.getAltAllele()[0];
            int normalALTReadCount = altAlleleNormalReadCount[altAllele];
            int tumorALTReadCount = altAlleleTumoralReadCount[altAllele];
            double normalALTReadFraction = normalALTReadCount > 0 ? (double) normalALTReadCount / normalDepth : 0.0;
            double tumorALTReadFraction = tumorALTReadCount > 0 ? (double) tumorALTReadCount / tumorDepth : 0.0;
            boolean hasGermlineLikeNormalTumorRatio = tumorALTReadCount == 0 || normalALTReadFraction / tumorALTReadFraction > 0.33;
            // Only count as germline if it exceeds the threshold in normal reads or overlaps uncovered position
            // This allows somatic variants (only present in reads affected by sequencing error) to passthrough intact
            boolean isGermline = (normalALTReadCount > readThreshold || hasGermlineLikeNormalTumorRatio) || overlapsUncoveredPosition(var);
            if (isGermline) {
                String readAlnId = var.getReadAlnName();
                AnonymizedRead anonymizedRead = anonymizedReadCache.get(readAlnId);
                if(anonymizedRead != null) anonymizedRead.addSignalToAnonymize(var);
            } else if (!alleleCaptured[altAllele & 0xFF] && tumorALTReadCount > 0) {
                // Capture original VAF data for this potential somatic site (used by fixVAF correction)
                alleleCaptured[altAllele & 0xFF] = true;
                normalALTReadFractions[altAllele & 0xFF] = normalALTReadFraction;
                tumorALTReadFractions[altAllele & 0xFF] = tumorALTReadFraction;
                candidateSomaticSignals.add(var);
            }
        }
        // After processing all signals, add the best potential somatic ALT to the somaticSiteVAFs map
        // Skip VAF correction entirely if original tumor depth is below the minimum threshold to avoid
        // false variant calls at low-coverage sites in the anonymized BAM
        if (tumorDepth < minDepthForVAFCorrection) return;
        double highestTumorVAF = 0.0;
        for  (Signal var : candidateSomaticSignals) {
            byte altAllele = var.getAltAllele()[0];
            int alleleIdx = altAllele & 0xFF;
            double normalVAF = normalALTReadFractions[alleleIdx];
            double tumorVAF = tumorALTReadFractions[alleleIdx];
            if (tumorVAF > highestTumorVAF) {
                byte medTumorQuality = computeMedianQuality(tumorAltQualities[alleleIdx]);
                // Fall back to tumor quality if no normal ALT reads exist at this (somatic) site
                byte medNormalQuality = normalAltQualities[alleleIdx] != null
                        ? computeMedianQuality(normalAltQualities[alleleIdx])
                        : medTumorQuality;
                somaticSiteVAFs.computeIfAbsent(var.getSequenceName(), k -> new TreeMap<>())
                    .put(var.getLocation(), new SomaticSiteVAF(altAllele, (float) normalVAF, (float) tumorVAF, medNormalQuality, medTumorQuality));
                highestTumorVAF = tumorVAF; // Only add one potential somatic site per position (the one with the highest tumor VAF)
            }
        }
    }

    private byte computeMedianQuality(List<Byte> qualities) {
        if (qualities == null || qualities.isEmpty()) return 0;
        List<Byte> sorted = new ArrayList<>(qualities);
        Collections.sort(sorted);
        int size = sorted.size();
        if (size % 2 == 1) return sorted.get(size / 2);
        return (byte) (((sorted.get((size - 1) / 2) & 0xFF) + (sorted.get(size / 2) & 0xFF)) / 2);
    }

    /**
     * Classify complex signals into potential germline variations by virtually inferring either a complete graph of
     * signals from the normal dataset, or a bipartite graph from the normal-tumor pair
     *
     */
    private void processSignals() {
        List<ReadSignal> cigarSignals = new ArrayList<>();
        List<ReadSignal> insertSizeSignals = new ArrayList<>();
        List<ReadSignal> strandOrientationSignals = new ArrayList<>();
        List<ReadSignal> chromChangeSignals = new ArrayList<>();
        // For each source of signal, create a partition, mixed for cigarSignals
        List<ReadSignal> signalsList = signals.getSignalsList();
        for(ReadSignal signal : signalsList){
            if(signal.isIndel()){
                cigarSignals.add(signal);
            }
            else if(signal.isSoftClip()){
                cigarSignals.add(signal);
            }
            else if(signal.getSource() == Signal.Source.INSERT_SIZE){
                insertSizeSignals.add(signal);
            }
            else if(signal.getSource() == Signal.Source.STRAND_ORIENTATION){
                strandOrientationSignals.add(signal);
            }
            else if(signal.getSource() == Signal.Source.CHROM_CHANGE){
                chromChangeSignals.add(signal);
            }
        }
        processTypeSignals(insertSizeSignals);
        processTypeSignals(strandOrientationSignals);
        processTypeSignals(chromChangeSignals);
        List<List<ReadSignal>> indelGraph = processTypeSignals(cigarSignals, SIGNAL_PER_PARTITION_LIMIT*2, true);
        rescueMissedSomaticIndelSignals(indelGraph);
    }

    private void rescueMissedSomaticIndelSignals(List<List<ReadSignal>> indelGraphs) {
        for(List<ReadSignal> indelPartitionGraph : indelGraphs) {
            if(indelPartitionGraph.size() < 2) {
                // If there are not enough signals to form a graph, skip processing
                continue;
            }
            // Handle large graphs by spatial partitioning with overlap
            if(indelPartitionGraph.size() > 10000) {

                int maxPartitionSize = 8000; // Leave room for overlap
                int overlapSize = 1000; // Overlap between partitions to catch boundary signals

                List<List<ReadSignal>> partitions = new ArrayList<>();
                int start = 0;

                while (start < indelPartitionGraph.size()) {
                    int end = Math.min(start + maxPartitionSize, indelPartitionGraph.size());

                    // If this isn't the last partition, extend to include overlap
                    if (end < indelPartitionGraph.size()) {
                        int overlapEnd = Math.min(end + overlapSize, indelPartitionGraph.size());
                        List<ReadSignal> partition = new ArrayList<>(indelPartitionGraph.subList(start, overlapEnd));
                        partitions.add(partition);
                        start += maxPartitionSize; // Move start by partition size, not including overlap
                    } else {
                        // Last partition - take everything remaining
                        List<ReadSignal> partition = new ArrayList<>(indelPartitionGraph.subList(start, end));
                        partitions.add(partition);
                        break;
                    }
                }

                // Process each partition with deduplication
                Set<String> processedSignals = new HashSet<>(); // Track processed signals to avoid duplicates from overlaps
                for (List<ReadSignal> partition : partitions) {
                    processIndelGraphPartition(partition, processedSignals);
                }
                continue;
            }
            // Process normally for smaller graphs
            processIndelGraphPartition(indelPartitionGraph, new HashSet<>());
        }
    }

    private void processIndelGraphPartition(List<ReadSignal> indelPartitionGraph, Set<String> processedSignals) {
        if (indelPartitionGraph.size() < 2) return;

        // Process the indel graph to declassify leaked somatic signals using hierarchical clustering
        ReadSignal[] indelGraphArray = indelPartitionGraph.toArray(new ReadSignal[0]);

        // Generate the linkage for hierarchical clustering, using the Euclidean distance between signals (location and length)
        Linkage linkage = WardLinkage.of(indelGraphArray, new SignalDistance());
        HierarchicalClustering hierarchicClusters = HierarchicalClustering.fit(linkage);

        // Check if there is a clear first cluster split, and rescue the signals if they are only somatic
        double[] heights = hierarchicClusters.height().clone();

        // Save signals from the two biggest clusters, if they are too far apart and one is only somatic
        double heightCutoff = mean(heights);
        if(heightCutoff == 0){
            // If the height cutoff is zero, it means that all signals are at identical, so skip processing
            return;
        }

        boolean isNonMonotonicTree = false;
        for (int i = 0; i < heights.length - 1; i++) {
            if (heights[i] > heights[i + 1]) {
                isNonMonotonicTree = true;
                break;
            }
        }
        if(isNonMonotonicTree) return;

        int[] clusters = hierarchicClusters.partition(heightCutoff);
        boolean[] isSomaticCluster = new boolean[clusters.length];
        Arrays.fill(isSomaticCluster, true);
        List<List<ReadSignal>> clusteredSignals = new ArrayList<>(clusters.length);
        for(int i = 0; i < clusters.length; i++) {
            clusteredSignals.add(new ArrayList<>());
        }
        Set<Integer> clusterSet = new HashSet<>();
        for(int i = 0; i < clusters.length; i++) {
            ReadSignal signal = indelGraphArray[i];
            int clusterId = clusters[i];
            if(signal.isFromNormalDataset()){
                isSomaticCluster[clusterId] = false;
            }
            clusteredSignals.get(clusterId).add(signal);
            clusterSet.add(clusterId);
        }

        for(int clusterId : clusterSet) {
            if(isSomaticCluster[clusterId]){
                List<ReadSignal> rescuedSomaticSignals = clusteredSignals.get(clusterId);
                for(ReadSignal signal : rescuedSomaticSignals) {
                    if (!signal.isClassifiedByDistance()) continue;
                    // If the signal is somatic, rescue it from anonymization
                    String readAlnName = signal.getReadAlnName();
                    // Use readAlnName to avoid processing the same signal multiple times from overlaps
                    if (processedSignals.contains(readAlnName)) continue;
                    processedSignals.add(readAlnName);
                    AnonymizedRead anonymizedRead = anonymizedReadCache.get(readAlnName);
                    if (anonymizedRead != null) anonymizedRead.rescueIndelSignalToAnonymize((Signal) signal);
                    signal.setIsGermline(false);
                    // If the signal is somatic, remove it from the onHoldReads map
                    manageAnalyzedSignalsInHoldedReadAln(signal);
                }
            }
        }
    }

    static class SignalDistance implements Distance<ReadSignal> {
        @Override
        public double d(ReadSignal a, ReadSignal b) {
            return Operations.computeTwoDimEuclideanDistance(a.getLocation(), b.getLocation(),
                    a.getLength(), b.getLength());
        }
    }

    private void processTypeSignals(List<ReadSignal> typeSignals) {
        processTypeSignals(typeSignals, AnonymizedReadAlignmentProvider.SIGNAL_PER_PARTITION_LIMIT, false);
    }

    private List<List<ReadSignal>> processTypeSignals(List<ReadSignal> typeSignals, int signalPerPartitionLimit, boolean returnConnectedSignals) {
        int n = typeSignals.size();
        List<ReadSignal> currentPartition = new ArrayList<>();
        List<List<ReadSignal>> connectedPartitionSignals = new ArrayList<>();
        for(int i = 0; i < n-1; i++){
            ReadSignal currentSignal = typeSignals.get(i);
            currentPartition.add(currentSignal);
            ReadSignal nextSignal = typeSignals.get(i+1);
            int locationDistance = Math.abs(nextSignal.getLocation() - currentSignal.getLocation());
            boolean lastSignalUnreachable = locationDistance > MAX_LOCATION_DISTANCE_THRESHOLD;
            if(lastSignalUnreachable || currentPartition.size() >= signalPerPartitionLimit || (i==n-2)){
                if(!lastSignalUnreachable){
                    currentPartition.add(nextSignal);
                    i++;
                }
                if(returnConnectedSignals) {
                    processPartition(currentPartition, true, connectedPartitionSignals);
                }
                else{
                    processPartition(currentPartition);
                }
                currentPartition = new ArrayList<>();
            }
        }
        return connectedPartitionSignals;
    }

    private void processPartition(List<ReadSignal> signals){
        processPartition(signals, false, new ArrayList<>());
    }

    /**
     * Process a partition of signals, classifying them as germline complex signals
     * @param signals
     */
    private void processPartition(List<ReadSignal> signals, boolean returnConnectedComponents, List<List<ReadSignal>> connectedSignals) {
        int n = signals.size();
        AdjacencyList signalsGraph = new AdjacencyList(n, false);
        for (int i = 0; i < n; i++){
            ReadSignal firstSignal = signals.get(i);
            for (int j = i + 1; j < n; j++){
                ReadSignal secondSignal = signals.get(j);
                //Avoid comparing signals from small deletions against insertions
                if(isDifferentIndelVariation(firstSignal, secondSignal)) continue;
                //Avoid comparing signals from the same read
                if(firstSignal.getReadAlnName().contains(secondSignal.getReadAlnName())) continue;
                double signalTypeThreshold = getSignalThreshold(firstSignal, secondSignal);
                double signalDistance = Operations.computeTwoDimEuclideanDistance(firstSignal.getLocation(), secondSignal.getLocation(),
                        firstSignal.getLength(), secondSignal.getLength());
                //Classify signals that are close enough to be considered as a germline signal
                if(signalDistance < signalTypeThreshold &&
                        (firstSignal.isFromNormalDataset() || secondSignal.isFromNormalDataset()) ) {
                    if(returnConnectedComponents){
                        signalsGraph.addEdge(i, j);
                        signalsGraph.addEdge(j, i);
                    }
                    long startclassifyPGcomplexSignal = System.currentTimeMillis();
                    classifyGermlineComplexSignal(firstSignal);
                    classifyGermlineComplexSignal(secondSignal);
                    firstSignal.setClassifiedByDistance(true);
                    secondSignal.setClassifiedByDistance(true);
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
                            ( firstSignal.isIntraAlignmentSignal() && firstSignal.isGermline() )){
                        classifyGermlineComplexSignal(firstSignal);
                    }
                    if(overlapsUncoveredPosition(secondSignal) ||
                            ( secondSignal.isIntraAlignmentSignal() && secondSignal.isGermline() )){
                        classifyGermlineComplexSignal(secondSignal);
                    }
                }
            }
            // Manage the read that has been analyzed, removing it from the onHoldReads map
            manageAnalyzedSignalsInHoldedReadAln(firstSignal);
        }
        if(returnConnectedComponents){
            // Add the connected signals to the list of connected components
            // Get the connected components from the graph
            int[][] connectedComponentsArray = signalsGraph.bfcc();
            for (int[] component : connectedComponentsArray) {
                List<ReadSignal> connectedComponentSignals = new ArrayList<>();
                for (int idx : component) {
                    ReadSignal signal = signals.get(idx);
                    if (signal instanceof  MultiSourceCIGARSignal multiSourceCIGARSignal){
                        for ( Map.Entry<String, Integer> entry : multiSourceCIGARSignal.getSources().entrySet()) {
                            connectedComponentSignals.add(multiSourceCIGARSignal.getReadAlnSignal(entry.getKey(), entry.getValue()));
                        }
                    }
                    else {
                        connectedComponentSignals.add(signal);
                    }
                }
                connectedSignals.add(connectedComponentSignals);
            }
        }
    }

    private void manageAnalyzedSignalsInHoldedReadAln(ReadSignal signal){
        if (signal instanceof MultiSourceCIGARSignal multiSourceSignal) {
            for (String readAlnName : multiSourceSignal.getSources().keySet()) {
                if(onHoldReads.containsKey(readAlnName)){
                    int remainingReadSignals = onHoldReads.get(readAlnName);
                    if(remainingReadSignals == 1){
                        onHoldReads.remove(readAlnName);
                    }
                    else{
                        onHoldReads.put(readAlnName, remainingReadSignals-1);
                    }
                }
            }
        }
        else {
            String readAlnName = signal.getReadAlnName();
            if(onHoldReads.containsKey(readAlnName)){
                int remainingReadSignals = onHoldReads.get(readAlnName);
                if(remainingReadSignals == 1){
                    onHoldReads.remove(readAlnName);
                }
                else{
                    onHoldReads.put(readAlnName, remainingReadSignals-1);
                }
            }
        }
    }

    /**
     * Calculate the gene panel threshold based on sequencing depth.
     * Uses 1% sequencing error rate with 1/4 probability of specific base error.
     * Threshold = error_rate * (1/4) * depth = 0.01 * 0.25 * depth = 0.0025 * depth
     */
    private int calculateGenePanelThreshold(int sequencingDepth) {
        // Sequencing error rate: 1% (0.01)
        return (int) Math.round(0.01 * sequencingDepth);
    }

    private double getSignalThreshold(ReadSignal firstSignal, ReadSignal secondSignal) {
        if (firstSignal.isIndel() && secondSignal.isIndel()){
            return INDEL_SIGNAL_THRESHOLD;
        }
        if(firstSignal.isSoftClip() && secondSignal.isSoftClip()){
            // If both signals are soft-clips, use a different threshold
            return SOFTCLIP_SIGNAL_THRESHOLD;
        }
        if ((firstSignal.isIndel() && secondSignal.isSoftClip())
                || (firstSignal.isSoftClip() && secondSignal.isIndel())){
            return INDEL_SOFTCLIP_SIGNAL_THRESHOLD;
        }
        if (Signal.Source.INSERT_SIZE == firstSignal.getSource() && Signal.Source.INSERT_SIZE == secondSignal.getSource()){
            return INSERT_SIZE_SIGNAL_THRESHOLD;
        }
        if (Signal.Source.STRAND_ORIENTATION == firstSignal.getSource() && Signal.Source.STRAND_ORIENTATION == secondSignal.getSource()){
            return STRAND_ORIENTATION_SIGNAL_THRESHOLD;
        }
        return SOFTCLIP_SIGNAL_THRESHOLD;
    }

    private boolean isDifferentIndelVariation(ReadSignal firstSignal, ReadSignal secondSignal) {
        if (firstSignal.isIndel() && secondSignal.isIndel()){
            return firstSignal.getIndelSignalType() != secondSignal.getIndelSignalType();
        }
        return false;
    }

    private void classifyGermlineComplexSignal(ReadSignal signal) {
        if( !signal.isHandled() ) {
            //The anonymizedRead is null if it comes from the offset before the first pileup position,
            //it is used for classification but is left to be returned by other thread
            if (signal instanceof MultiSourceCIGARSignal multiSourceSignal) {
                for (Map.Entry<String, Integer> entry : multiSourceSignal.getSources().entrySet()){
                    String readAlnId = entry.getKey();
                    AnonymizedRead anonymizedRead = anonymizedReadCache.get(readAlnId);
                    if (anonymizedRead != null) {
                        Signal readSpecificSignal = multiSourceSignal.getReadAlnSignal(readAlnId, entry.getValue());
                        anonymizedRead.addSignalToAnonymize(readSpecificSignal);
                        readSpecificSignal.setIsGermline(true);
                        readSpecificSignal.setHandled(true);
                    }
                }
            }
            else {
                String readAlnId = signal.getReadAlnName();
                AnonymizedRead anonymizedRead = anonymizedReadCache.get(readAlnId);
                if(anonymizedRead != null) anonymizedRead.addSignalToAnonymize((Signal) signal);
            }
            signal.setIsGermline(true);
            signal.setHandled(true);
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

    class SignalCollection{

        private Map<String, ReadSignal> normalSignals;
        private Map<String, ReadSignal> tumorSignals;
        public SignalCollection() {
            init();
        }

        private void init() {
            normalSignals = new LinkedHashMap<>(1_000_000);
            tumorSignals = new LinkedHashMap<>(1_000_000);
        }

        public void add(ReadSignal signal) {
            Map<String, ReadSignal> signals = signal.isFromNormalDataset() ? normalSignals : tumorSignals;
            if (signal.isInterAlignmentSignal()){
                signals.put(signal.getSignalKey(), signal);
                return;
            }
            boolean existsInNormal = normalSignals.containsKey(signal.getSignalKey());
            boolean existsInTumor = tumorSignals.containsKey(signal.getSignalKey());
            if(existsInNormal || existsInTumor) {
                ReadSignal existingSignal = existsInNormal ? normalSignals.get(signal.getSignalKey()) : tumorSignals.get(signal.getSignalKey());
                boolean comeFromSameSourceDataset = signal.isFromNormalDataset() == existingSignal.isFromNormalDataset();
                if(existsInNormal && existsInTumor) {
                    existingSignal = signal.isFromNormalDataset() ? normalSignals.get(signal.getSignalKey()) : tumorSignals.get(signal.getSignalKey());
                    addToExistingSignals(signal, signals, existingSignal);
                    signal.setIsGermline(true);
                }
                else if(comeFromSameSourceDataset) {
                    addToExistingSignals(signal, signals, existingSignal);
                }
                else {
                    signals.put(signal.getSignalKey(), signal);
                }
                if(signal.isFromNormalDataset() || existingSignal.isFromNormalDataset()) {
                    signal.setIsGermline(true);
                    existingSignal.setIsGermline(true);
                }
            }
            else {
                signals.put(signal.getSignalKey(), signal);
            }
        }

        private void addToExistingSignals(ReadSignal signal, Map<String, ReadSignal> signals, ReadSignal existingSignal) {
            if(existingSignal instanceof MultiSourceCIGARSignal existingMultiSourceSignal) {
                existingMultiSourceSignal.addSource(signal.getReadAlnName(), signal.getInReadPosition());
            }
            else {
                MultiSourceCIGARSignal multiSourceCIGARSignal = new MultiSourceCIGARSignal(existingSignal);
                multiSourceCIGARSignal.addSource(signal.getReadAlnName(), signal.getInReadPosition());
                signals.replace(existingSignal.getSignalKey(), multiSourceCIGARSignal);
            }
        }

        public List<ReadSignal> getSignalsList() {
            List<ReadSignal> answer = new ArrayList<>(size());
            Iterator<ReadSignal> normalIt = normalSignals.values().iterator();
            Iterator<ReadSignal> tumorIt = tumorSignals.values().iterator();
            ReadSignal normal = normalIt.hasNext() ? normalIt.next() : null;
            ReadSignal tumor = tumorIt.hasNext() ? tumorIt.next() : null;
            while (normal != null || tumor != null) {
                if (normal == null) {
                    answer.add(tumor);
                    tumor = tumorIt.hasNext() ? tumorIt.next() : null;
                }
                else if (tumor == null) {
                    answer.add(normal);
                    normal = normalIt.hasNext() ? normalIt.next() : null;
                }
                else if (normal.getLocation() <= tumor.getLocation()) {
                    answer.add(normal);
                    normal = normalIt.hasNext() ? normalIt.next() : null;
                }
                else {
                    answer.add(tumor);
                    tumor = tumorIt.hasNext() ? tumorIt.next() : null;
                }
            }
            return answer;
        }

        public int size() {
            return normalSignals.size() + tumorSignals.size();
        }

        public void clear() {
            normalSignals.clear();
            tumorSignals.clear();
        }
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
