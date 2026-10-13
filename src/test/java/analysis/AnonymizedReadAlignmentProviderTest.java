package analysis;

import genomicelements.AnonymizedRead;
import genomicelements.GenomicRegionBaseImpl;
import genomicelements.SomaticSiteVAF;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import testutils.TestAlignments;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Classification of variation over a small tumor-normal pair built by hand, where every read is known:
 * pairs of 50 bp reads tiled every 10 bp, a germline SNV, a somatic SNV, a SNV where only the tumor has coverage,
 * a germline deletion and a somatic deletion
 */
public class AnonymizedReadAlignmentProviderTest {

    private static final String CONTIG = "chr1";
    private static final int CONTIG_LENGTH = 3000;
    private static final int READ_LENGTH = 50;
    private static final int MATE_OFFSET = 150;
    private static final int TILE_STEP = 10;
    // Region covered by both samples and region covered only by the tumor (1-based starts of the first reads)
    private static final int SHARED_FIRST_START = 101;
    private static final int SHARED_LAST_START = 1500;
    private static final int TUMOR_ONLY_FIRST_START = 2001;
    private static final int TUMOR_ONLY_LAST_START = 2500;

    private static final int GERMLINE_SNV = 500;
    private static final int SOMATIC_SNV = 800;
    private static final int TUMOR_ONLY_SNV = 2300;
    // First deleted base of two deletions of 3 bp
    private static final int GERMLINE_DELETION = 1000;
    private static final int SOMATIC_DELETION = 1300;
    private static final int DELETION_LENGTH = 3;
    private static final int MIN_DELETION_FLANK = 8;

    @TempDir
    File tempDir;

    private String reference;
    private SAMFileHeader header;
    private File referenceFile;
    private File normalFile;
    private File tumorFile;
    private List<SAMRecord> normalRecords;
    private List<SAMRecord> tumorRecords;

    @BeforeEach
    void setUp() throws IOException {
        reference = TestAlignments.randomSequence(CONTIG_LENGTH, 17);
        header = TestAlignments.header(new String[]{CONTIG}, new int[]{CONTIG_LENGTH});
        referenceFile = TestAlignments.writeReference(tempDir, "ref.fa", new String[]{CONTIG}, new String[]{reference});
        normalRecords = new ArrayList<>();
        tumorRecords = new ArrayList<>();
        for (int start = SHARED_FIRST_START; start <= SHARED_LAST_START; start += TILE_STEP) {
            addPair(normalRecords, "n" + start, start, new int[]{GERMLINE_SNV}, new int[]{GERMLINE_DELETION});
            addPair(tumorRecords, "t" + start, start, new int[]{GERMLINE_SNV, SOMATIC_SNV}, new int[]{GERMLINE_DELETION, SOMATIC_DELETION});
        }
        for (int start = TUMOR_ONLY_FIRST_START; start <= TUMOR_ONLY_LAST_START; start += TILE_STEP) {
            addPair(tumorRecords, "t" + start, start, new int[]{TUMOR_ONLY_SNV}, new int[0]);
        }
        normalFile = TestAlignments.writeIndexedBam(tempDir, "normal.bam", header, normalRecords);
        tumorFile = TestAlignments.writeIndexedBam(tempDir, "tumor.bam", header, tumorRecords);
    }

    /** Adds a pair carrying the given variants, unless one of its reads ends too close to a deletion */
    private void addPair(List<SAMRecord> records, String name, int start, int[] snvs, int[] deletions) {
        SAMRecord[] pair = new SAMRecord[2];
        int[] starts = {start, start + MATE_OFFSET};
        for (int i = 0; i < 2; i++) {
            int readStart = starts[i];
            String cigar = READ_LENGTH + "M";
            String bases = reference.substring(readStart - 1, readStart - 1 + READ_LENGTH);
            for (int deletion : deletions) {
                int basesBefore = deletion - readStart;
                if (basesBefore >= READ_LENGTH || basesBefore + DELETION_LENGTH <= 0) continue;
                if (basesBefore < MIN_DELETION_FLANK || READ_LENGTH - basesBefore < MIN_DELETION_FLANK) return;
                cigar = basesBefore + "M" + DELETION_LENGTH + "D" + (READ_LENGTH - basesBefore) + "M";
                bases = reference.substring(readStart - 1, deletion - 1)
                        + reference.substring(deletion - 1 + DELETION_LENGTH, readStart - 1 + READ_LENGTH + DELETION_LENGTH);
            }
            pair[i] = TestAlignments.read(header, name, CONTIG, readStart, cigar, bases, i == 0, i == 1, CONTIG,
                    starts[1 - i], i == 0, i == 0 ? MATE_OFFSET + READ_LENGTH : -(MATE_OFFSET + READ_LENGTH));
            for (int snv : snvs) {
                TestAlignments.setBaseAtReferencePosition(pair[i], snv, TestAlignments.otherBase(reference.charAt(snv - 1)));
            }
        }
        records.add(pair[0]);
        records.add(pair[1]);
    }

    private AnonymizedReadAlignmentProvider newProvider() {
        AnonymizedReadAlignmentProvider provider = new AnonymizedReadAlignmentProvider();
        provider.setReadsToExclude(new HashSet<>());
        // Wide thresholds: no pair of this dataset has an abnormal insert size
        provider.setInsertSizeMinThreshold(0);
        provider.setInsertSizeMaxThreshold(10_000);
        provider.setRefSequence(reference.getBytes(StandardCharsets.UTF_8));
        return provider;
    }

    /** Runs the provider over the whole contig and returns the anonymized records of normal [0] and tumor [1] */
    private List<List<SAMRecord>> anonymize(AnonymizedReadAlignmentProvider provider) throws IOException {
        List<List<SAMRecord>> answer = List.of(new ArrayList<>(), new ArrayList<>());
        provider.init(normalFile.getPath(), tumorFile.getPath(), referenceFile.getPath(), new GenomicRegionBaseImpl(CONTIG, 0, 1, CONTIG_LENGTH));
        for (AnonymizedRead read : provider) {
            if (!read.isAnonymized()) read.anonymizeRead();
            answer.get(read.isFromNormalDataset() ? 0 : 1).add(read.getAnonymizedSamRecord());
        }
        provider.close();
        return answer;
    }

    private int countAlt(List<SAMRecord> records, int position) {
        char refBase = reference.charAt(position - 1);
        int answer = 0;
        for (SAMRecord record : records) {
            int readPosition = record.getReadPositionAtReferencePosition(position);
            if (readPosition > 0 && record.getReadString().charAt(readPosition - 1) != refBase) answer++;
        }
        return answer;
    }

    private int countDeletions(List<SAMRecord> records, int firstDeletedBase) {
        int answer = 0;
        for (SAMRecord record : records) {
            if (!record.getCigar().containsOperator(CigarOperator.D)) continue;
            if (record.getAlignmentStart() < firstDeletedBase && record.getAlignmentEnd() > firstDeletedBase) answer++;
        }
        return answer;
    }

    @Test
    void everyReadIsReturnedOnce() throws IOException {
        List<List<SAMRecord>> anonymized = anonymize(newProvider());
        assertEquals(normalRecords.size(), anonymized.get(0).size());
        assertEquals(tumorRecords.size(), anonymized.get(1).size());
    }

    @Test
    void germlineSNVIsRemovedFromBothSamples() throws IOException {
        assertTrue(countAlt(normalRecords, GERMLINE_SNV) > 0);
        assertTrue(countAlt(tumorRecords, GERMLINE_SNV) > 0);
        List<List<SAMRecord>> anonymized = anonymize(newProvider());
        assertEquals(0, countAlt(anonymized.get(0), GERMLINE_SNV));
        assertEquals(0, countAlt(anonymized.get(1), GERMLINE_SNV));
    }

    @Test
    void somaticSNVIsKeptInTheTumor() throws IOException {
        int inputAlt = countAlt(tumorRecords, SOMATIC_SNV);
        assertTrue(inputAlt > 0);
        List<List<SAMRecord>> anonymized = anonymize(newProvider());
        assertEquals(inputAlt, countAlt(anonymized.get(1), SOMATIC_SNV));
        assertEquals(0, countAlt(anonymized.get(0), SOMATIC_SNV));
    }

    @Test
    void variantWhereTheNormalHasNoCoverageIsRemoved() throws IOException {
        assertTrue(countAlt(tumorRecords, TUMOR_ONLY_SNV) > 0);
        List<List<SAMRecord>> anonymized = anonymize(newProvider());
        assertEquals(0, countAlt(anonymized.get(1), TUMOR_ONLY_SNV));
    }

    @Test
    void germlineDeletionIsRemovedFromBothSamplesAndReadsKeepTheirLength() throws IOException {
        assertTrue(countDeletions(normalRecords, GERMLINE_DELETION) > 0);
        assertTrue(countDeletions(tumorRecords, GERMLINE_DELETION) > 0);
        List<List<SAMRecord>> anonymized = anonymize(newProvider());
        assertEquals(0, countDeletions(anonymized.get(0), GERMLINE_DELETION));
        assertEquals(0, countDeletions(anonymized.get(1), GERMLINE_DELETION));
        for (SAMRecord record : anonymized.get(0)) {
            assertEquals(READ_LENGTH, record.getReadLength(), record.getReadName());
            assertEquals(0, record.getIntegerAttribute("NM") == null ? 0 : record.getIntegerAttribute("NM"), record.getReadName());
        }
    }

    @Test
    void somaticDeletionIsKeptInTheTumor() throws IOException {
        int inputDeletions = countDeletions(tumorRecords, SOMATIC_DELETION);
        assertTrue(inputDeletions > 0);
        List<List<SAMRecord>> anonymized = anonymize(newProvider());
        assertEquals(inputDeletions, countDeletions(anonymized.get(1), SOMATIC_DELETION));
    }

    @Test
    void genePanelModeAlsoRemovesGermlineAndRecordsSomaticSites() throws IOException {
        AnonymizedReadAlignmentProvider provider = newProvider();
        provider.setSampleType(GenomeAnonymizer.SAMPLE_TYPE_GENE_PANEL);
        provider.setMinDepthForVAFCorrection(1);
        int inputAlt = countAlt(tumorRecords, SOMATIC_SNV);
        List<List<SAMRecord>> anonymized = anonymize(provider);
        assertEquals(0, countAlt(anonymized.get(0), GERMLINE_SNV));
        assertEquals(0, countAlt(anonymized.get(1), GERMLINE_SNV));
        assertEquals(0, countAlt(anonymized.get(1), TUMOR_ONLY_SNV));
        assertEquals(inputAlt, countAlt(anonymized.get(1), SOMATIC_SNV));
        Map<String, TreeMap<Integer, SomaticSiteVAF>> somaticSites = provider.getSomaticSiteVAFs();
        assertNotNull(somaticSites.get(CONTIG));
        SomaticSiteVAF site = somaticSites.get(CONTIG).get(SOMATIC_SNV);
        assertNotNull(site, "The somatic SNV was not recorded for VAF correction");
        assertEquals((byte) TestAlignments.otherBase(reference.charAt(SOMATIC_SNV - 1)), site.altAllele());
        assertEquals(1.0f, site.tumorVAF(), 1e-6);
        assertEquals(0.0f, site.normalVAF(), 1e-6);
        assertEquals(TestAlignments.DEFAULT_BASE_QUALITY, site.medianTumorAltQuality());
        // Germline sites are not candidates for VAF correction
        assertTrue(!somaticSites.get(CONTIG).containsKey(GERMLINE_SNV));
    }
}
