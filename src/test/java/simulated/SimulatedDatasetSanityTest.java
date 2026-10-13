package simulated;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInstance;
import simulated.SimulatedDataset.Region;
import simulated.SimulatedDataset.TruthVariant;

import java.util.HashSet;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Checks that the well-formed simulated dataset has the properties that the validation relies on, so a regenerated
 * dataset that lost one of them is detected before its results are interpreted
 * @author Nicolas Gaitan
 */
@Tag(SimulatedTestConfig.SIMULATED_TAG)
@TestInstance(TestInstance.Lifecycle.PER_CLASS)
public class SimulatedDatasetSanityTest {

    public static final int MIN_CONTIGS = 4;
    public static final int MIN_GENOME_SIZE = 50_000;
    public static final int MIN_MEDIAN_DEPTH = 8;
    public static final int MAX_MEDIAN_DEPTH = 12;
    public static final int MAX_DEPTH = 40;
    public static final double MAX_MISMATCH_RATE = 0.01;
    public static final double MIN_MEAN_BASE_QUALITY = 30;
    public static final double MIN_UNIQUELY_MAPPED_FRACTION = 0.9;
    // Reads are sampled as whole fragments, so the edges of the coverage regions are not sharp
    public static final int REGION_EDGE_MARGIN = 100;

    private SimulatedDataset dataset;
    private AlignmentEvidence normal;
    private AlignmentEvidence tumor;

    @BeforeAll
    void loadDataset() throws Exception {
        dataset = new SimulatedDataset(SimulatedTestConfig.getWellformedDir());
        normal = AlignmentEvidence.load(dataset.getNormalFile(), dataset.getReference(), false);
        tumor = AlignmentEvidence.load(dataset.getTumorFile(), dataset.getReference(), false);
    }

    @Test
    void referenceHasEnoughContigsAndBases() {
        assertTrue(dataset.getReference().size() >= MIN_CONTIGS, "Contigs: " + dataset.getReference().size());
        assertTrue(dataset.getGenomeSize() >= MIN_GENOME_SIZE, "Genome size: " + dataset.getGenomeSize());
    }

    @Test
    void medianDepthIsAroundTenInBothSamples() {
        for (AlignmentEvidence sample : new AlignmentEvidence[]{normal, tumor}) {
            int[] depths = sample.getSortedDepths();
            int median = depths[depths.length / 2];
            assertTrue(median >= MIN_MEDIAN_DEPTH && median <= MAX_MEDIAN_DEPTH, "Median depth: " + median);
            assertTrue(depths[depths.length - 1] <= MAX_DEPTH, "Maximum depth: " + depths[depths.length - 1]);
        }
    }

    @Test
    void readsResembleIlluminaSequencing() {
        for (AlignmentEvidence sample : new AlignmentEvidence[]{normal, tumor}) {
            long bases = 0;
            long mismatches = 0;
            long qualitySum = 0;
            int uniquelyMapped = 0;
            int ambiguouslyMapped = 0;
            int primary = 0;
            int duplicates = 0;
            for (SAMRecord record : sample.getRecords()) {
                if (record.isSecondaryOrSupplementary()) continue;
                primary++;
                assertTrue(record.getReadPairedFlag(), "Unpaired read " + record.getReadName());
                if (record.getDuplicateReadFlag()) duplicates++;
                if (record.getMappingQuality() == 60) uniquelyMapped++;
                if (record.getMappingQuality() == 0) ambiguouslyMapped++;
                bases += record.getReadLength();
                mismatches += SequenceUtil.countMismatches(record, dataset.getReference().get(record.getContig()), 0);
                for (byte quality : record.getBaseQualities()) qualitySum += quality;
            }
            double mismatchRate = (double) mismatches / bases;
            // Mismatches include the spiked-in variants, so the sequencing error rate is below this value
            assertTrue(mismatchRate <= MAX_MISMATCH_RATE, "Mismatch rate: " + mismatchRate);
            assertTrue(mismatchRate > 0, "Reads without any sequencing error");
            assertTrue((double) qualitySum / bases >= MIN_MEAN_BASE_QUALITY, "Mean base quality: " + (double) qualitySum / bases);
            assertTrue((double) uniquelyMapped / primary >= MIN_UNIQUELY_MAPPED_FRACTION, "Reads with MAPQ 60: " + uniquelyMapped + " of " + primary);
            assertTrue(ambiguouslyMapped > 0, "No reads with MAPQ 0, the repetitive regions are not represented");
            assertTrue(duplicates > 0, "No duplicate reads");
        }
    }

    @Test
    void coverageRegionsAreCoveredAsDeclared() {
        int checked = 0;
        for (Region region : dataset.getRegions()) {
            if (!SimulatedDataset.REGION_KIND_COVERAGE.equals(region.kind())) continue;
            boolean normalUncovered = region.details().contains("normal=0.0");
            boolean tumorUncovered = region.details().contains("tumor=0.0");
            int normalBases = 0;
            int tumorBases = 0;
            for (int pos = region.start() + 1 + REGION_EDGE_MARGIN; pos <= region.end() - REGION_EDGE_MARGIN; pos++) {
                normalBases += normal.getDepth(region.contig(), pos);
                tumorBases += tumor.getDepth(region.contig(), pos);
            }
            if (normalUncovered) assertEquals(0, normalBases, "Normal coverage in region " + region.name());
            else assertTrue(normalBases > 0, "No normal coverage in region " + region.name());
            if (tumorUncovered) assertEquals(0, tumorBases, "Tumor coverage in region " + region.name());
            else assertTrue(tumorBases > 0, "No tumor coverage in region " + region.name());
            checked++;
        }
        assertTrue(checked >= 5, "Coverage regions checked: " + checked);
    }

    @Test
    void truthTableHasEveryKindOfVariant() {
        Set<String> kinds = new HashSet<>();
        Set<String> tumorVafs = new HashSet<>();
        for (TruthVariant variant : dataset.getTruthVariants()) {
            kinds.add(variant.origin() + ":" + variant.type());
            if (!variant.isGermline()) tumorVafs.add(String.valueOf(variant.expectedTumorVaf()));
        }
        for (String origin : new String[]{SimulatedDataset.ORIGIN_GERMLINE, SimulatedDataset.ORIGIN_SOMATIC}) {
            for (String type : new String[]{"SNV", "INS", "DEL", "SV_DEL", "SV_INS", "SV_INV", "SV_TRA"}) {
                assertTrue(kinds.contains(origin + ":" + type), "No " + origin + " variant of type " + type);
            }
        }
        assertTrue(tumorVafs.size() >= 4, "Different expected VAFs of somatic variants: " + tumorVafs);
    }

    @Test
    void structuralVariantsLeaveAbnormalPairs() {
        int interChromosomal = 0;
        int sameStrand = 0;
        int largeInsert = 0;
        for (SAMRecord record : tumor.getRecords()) {
            if (record.isSecondaryOrSupplementary() || record.getMateUnmappedFlag()) continue;
            if (!record.getReferenceIndex().equals(record.getMateReferenceIndex())) interChromosomal++;
            else if (record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag()) sameStrand++;
            else if (Math.abs(record.getInferredInsertSize()) >= OutputValidator.LARGE_INSERT_SIZE) largeInsert++;
        }
        assertTrue(interChromosomal > 0, "No pairs mapped to different contigs");
        assertTrue(sameStrand > 0, "No pairs with both reads in the same strand");
        assertTrue(largeInsert > 0, "No pairs with large insert size");
    }
}
