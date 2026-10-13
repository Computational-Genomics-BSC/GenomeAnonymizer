package simulated;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import simulated.SimulatedDataset.Region;
import simulated.SimulatedDataset.TruthVariant;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import static simulated.AlignmentEvidence.DELETION;
import static simulated.AlignmentEvidence.INSERTION;

/**
 * Compares the input and the output of a GenomeAnonymizer run on a simulated dataset. It produces three tables:
 * the findings of a scan of the whole genome that does not use the truth (germline leaks and artifacts), the fate
 * of every spiked-in small variant with its VAFs before and after, and the fate of the signals of every spiked-in
 * structural variant.
 *
 * The expectations follow the rule of the WGS mode: an allele seen in at least one eligible read of the normal
 * sample, or located where one of the two samples has no coverage, must not survive in any output. Everything else
 * seen in the tumor is expected to be kept untouched.
 * @author Nicolas Gaitan
 */
public class OutputValidator {

    public static final String CATEGORY_LEAK = "LEAK";
    public static final String CATEGORY_ARTIFACT = "ARTIFACT";

    public static final String REASON_ALT_IN_NORMAL = "ALT_IN_NORMAL";
    public static final String REASON_NEAR_NORMAL_INDEL = "NEAR_NORMAL_INDEL";
    public static final String REASON_UNCOVERED_IN_NORMAL = "UNCOVERED_IN_NORMAL";
    public static final String REASON_UNCOVERED_IN_TUMOR = "UNCOVERED_IN_TUMOR";
    public static final String REASON_NOT_IN_INPUT = "NOT_IN_INPUT";

    public static final String EXPECT_REMOVE_SEEN_IN_NORMAL = "MUST_REMOVE_SEEN_IN_NORMAL";
    public static final String EXPECT_REMOVE_UNCOVERED = "MUST_REMOVE_UNCOVERED";
    public static final String EXPECT_KEEP = "EXPECT_KEEP";
    public static final String EXPECT_NOT_SAMPLED = "NOT_SAMPLED";
    public static final String EXPECT_UNDETERMINED = "UNDETERMINED";

    public static final String VERDICT_OK = "OK";
    public static final String VERDICT_LEAK = "LEAK";
    public static final String VERDICT_SOMATIC_LOST = "SOMATIC_LOST";
    public static final String VERDICT_SOMATIC_REDUCED = "SOMATIC_REDUCED";
    public static final String VERDICT_ALT_INCREASED = "ALT_INCREASED";
    public static final String VERDICT_ARTIFACT = "ARTIFACT";
    // Germline variant that by chance was not sampled in the normal, so no method can tell it apart from a somatic one
    public static final String VERDICT_GERMLINE_UNSEEN_KEPT = "GERMLINE_UNSEEN_IN_NORMAL_KEPT";
    public static final String VERDICT_RESIDUAL_SIGNAL = "RESIDUAL_SIGNAL";

    // Aligners may place the same indel some bases apart in different reads
    public static final int INDEL_POSITION_TOLERANCE = 15;
    // Two indels of the same type closer than this in position and length are the same event for the anonymizer
    public static final int NEAR_INDEL_DISTANCE = 5;

    // Abnormal pair definitions, derived from the simulated library (insert size 420 +- 70, reads of 150 bp)
    public static final int LARGE_INSERT_SIZE = 700;
    public static final int SMALL_INSERT_SIZE = 300;
    public static final int PAIR_SIGNAL_WINDOW = 500;
    public static final int CLIP_SIGNAL_WINDOW = 20;
    public static final int SUPPLEMENTARY_SIGNAL_WINDOW = 150;
    public static final int MIN_CLIP_LENGTH = 10;

    public static final String NOT_AVAILABLE = "NA";

    private final SimulatedDataset dataset;
    // Inputs as GenomeAnonymizer sees them (eligible reads) and complete
    private final AlignmentEvidence normalIn;
    private final AlignmentEvidence tumorIn;
    private final AlignmentEvidence normalInAll;
    private final AlignmentEvidence tumorInAll;
    private final AlignmentEvidence normalOut;
    private final AlignmentEvidence tumorOut;

    public OutputValidator(SimulatedDataset dataset, AlignmentEvidence normalIn, AlignmentEvidence tumorIn,
                           AlignmentEvidence normalInAll, AlignmentEvidence tumorInAll,
                           AlignmentEvidence normalOut, AlignmentEvidence tumorOut) {
        this.dataset = dataset;
        this.normalIn = normalIn;
        this.tumorIn = tumorIn;
        this.normalInAll = normalInAll;
        this.tumorInAll = tumorInAll;
        this.normalOut = normalOut;
        this.tumorOut = tumorOut;
    }

    public static String formatVaf(int alt, int depth) {
        if (depth == 0) return NOT_AVAILABLE;
        return String.format(Locale.ROOT, "%.3f", (double) alt / depth);
    }

    // ---------------------------------------------------------------------------------------------------------
    // Genome wide scan, independent of the truth table
    // ---------------------------------------------------------------------------------------------------------

    /**
     * Looks at every non reference allele and every indel present in the outputs and reports those that should not
     * be there: alleles that the normal input also has, alleles where one of the samples has no coverage, and
     * alleles that the input of the same sample did not have (artifacts created by the anonymization)
     */
    public Table scanGenome() {
        Table findings = new Table("key", "category", "reason", "contig", "pos", "kind", "ref", "alt", "length",
                "normal_in_alt", "normal_in_depth", "tumor_in_alt", "tumor_in_depth", "normal_out_alt", "tumor_out_alt", "truth_id");
        for (Map.Entry<String, byte[]> entry : dataset.getReference().entrySet()) {
            String contig = entry.getKey();
            byte[] sequence = entry.getValue();
            for (int pos = 1; pos <= sequence.length; pos++) {
                char ref = (char) sequence[pos - 1];
                // Positions with ambiguous reference bases cannot be evaluated
                if (AlignmentEvidence.BASES.indexOf(ref) < 0 || ref == 'N') continue;
                for (char alt : new char[]{'A', 'C', 'G', 'T'}) {
                    if (alt == ref) continue;
                    int normalOutAlt = normalOut.getBaseCount(contig, pos, alt);
                    int tumorOutAlt = tumorOut.getBaseCount(contig, pos, alt);
                    if (normalOutAlt + tumorOutAlt == 0) continue;
                    int normalInAlt = normalIn.getBaseCount(contig, pos, alt);
                    int tumorInAlt = tumorIn.getBaseCount(contig, pos, alt);
                    int normalInDepth = normalIn.getDepth(contig, pos);
                    int tumorInDepth = tumorIn.getDepth(contig, pos);
                    boolean createdInNormal = normalOutAlt > 0 && normalInAll.getBaseCount(contig, pos, alt) == 0;
                    boolean createdInTumor = tumorOutAlt > 0 && tumorInAll.getBaseCount(contig, pos, alt) == 0;
                    String category = CATEGORY_LEAK;
                    String reason;
                    if (normalInAlt > 0) reason = REASON_ALT_IN_NORMAL;
                    else if (normalInDepth == 0) reason = REASON_UNCOVERED_IN_NORMAL;
                    else if (tumorInDepth == 0) reason = REASON_UNCOVERED_IN_TUMOR;
                    else if (createdInNormal || createdInTumor) {
                        category = CATEGORY_ARTIFACT;
                        reason = REASON_NOT_IN_INPUT;
                    }
                    else continue;
                    String key = contig + ":" + pos + ":SNV:" + ref + ">" + alt;
                    findings.addRow(key, category, reason, contig, pos, "SNV", ref, alt, 1, normalInAlt, normalInDepth,
                            tumorInAlt, tumorInDepth, normalOutAlt, tumorOutAlt, findTruthId(contig, pos, "SNV", String.valueOf(alt), 1));
                }
            }
        }
        Set<String> outputIndels = new TreeSet<>(normalOut.getIndelCounts().keySet());
        outputIndels.addAll(tumorOut.getIndelCounts().keySet());
        for (String indelKey : outputIndels) {
            String[] elements = indelKey.split(":");
            String contig = elements[0];
            int anchor = Integer.parseInt(elements[1]);
            char type = elements[2].charAt(0);
            int length = Integer.parseInt(elements[3]);
            int normalOutAlt = normalOut.getIndelCount(contig, anchor, type, length);
            int tumorOutAlt = tumorOut.getIndelCount(contig, anchor, type, length);
            int normalInAlt = normalIn.getIndelCount(contig, anchor, type, length);
            int tumorInAlt = tumorIn.getIndelCount(contig, anchor, type, length);
            int normalInDepth = normalIn.getDepth(contig, anchor);
            int tumorInDepth = tumorIn.getDepth(contig, anchor);
            boolean createdInNormal = normalOutAlt > 0 && normalInAll.getIndelCount(contig, anchor, type, length) == 0;
            boolean createdInTumor = tumorOutAlt > 0 && tumorInAll.getIndelCount(contig, anchor, type, length) == 0;
            String category = CATEGORY_LEAK;
            String reason;
            if (normalInAlt > 0) reason = REASON_ALT_IN_NORMAL;
            else if (hasNearIndel(normalIn, contig, anchor, type, length)) reason = REASON_NEAR_NORMAL_INDEL;
            else if (normalInDepth == 0) reason = REASON_UNCOVERED_IN_NORMAL;
            else if (tumorInDepth == 0) reason = REASON_UNCOVERED_IN_TUMOR;
            else if (createdInNormal || createdInTumor) {
                category = CATEGORY_ARTIFACT;
                reason = REASON_NOT_IN_INPUT;
            }
            else continue;
            String kind = type == INSERTION ? "INS" : "DEL";
            String key = contig + ":" + anchor + ":" + kind + ":" + length;
            findings.addRow(key, category, reason, contig, anchor, kind, ".", ".", length, normalInAlt, normalInDepth,
                    tumorInAlt, tumorInDepth, normalOutAlt, tumorOutAlt, findTruthId(contig, anchor, kind, null, length));
        }
        return findings;
    }

    private boolean hasNearIndel(AlignmentEvidence evidence, String contig, int anchor, char type, int length) {
        for (int otherLength = Math.max(1, length - NEAR_INDEL_DISTANCE); otherLength <= length + NEAR_INDEL_DISTANCE; otherLength++) {
            if (evidence.getIndelCount(contig, anchor - NEAR_INDEL_DISTANCE, anchor + NEAR_INDEL_DISTANCE, type, otherLength) > 0) {
                return true;
            }
        }
        return false;
    }

    private String findTruthId(String contig, int pos, String kind, String altBase, int length) {
        for (TruthVariant variant : dataset.getTruthVariants()) {
            if (variant.isStructural() || !variant.type().equals(kind) || !variant.contig().equals(contig)) continue;
            if (variant.isSNV()) {
                if (variant.pos() == pos && variant.alt().equals(altBase)) return variant.id();
            }
            else if (variant.length() == length) {
                int[] window = getIndelWindow(variant);
                if (pos >= window[0] && pos <= window[1]) return variant.id();
            }
        }
        return ".";
    }

    // ---------------------------------------------------------------------------------------------------------
    // Spiked-in SNVs and indels
    // ---------------------------------------------------------------------------------------------------------

    /**
     * Interval of anchor positions where the indel of a truth variant can be found. Inside tandem repeats the
     * aligner can place it anywhere along the repeat
     */
    private int[] getIndelWindow(TruthVariant variant) {
        int first = variant.pos() - INDEL_POSITION_TOLERANCE;
        int last = variant.pos() + INDEL_POSITION_TOLERANCE;
        for (Region region : dataset.getRegions()) {
            if (!SimulatedDataset.REGION_KIND_REFERENCE.equals(region.kind())) continue;
            if (region.name().startsWith("segdup") || region.name().equals("n_block")) continue;
            if (region.overlaps(variant.contig(), variant.pos(), variant.end())) {
                first = Math.min(first, region.start() + 1 - INDEL_POSITION_TOLERANCE);
                last = Math.max(last, region.end() + INDEL_POSITION_TOLERANCE);
            }
        }
        return new int[]{first, last};
    }

    private int[] countAllele(AlignmentEvidence evidence, TruthVariant variant) {
        int alt;
        if (variant.isSNV()) {
            alt = evidence.getBaseCount(variant.contig(), variant.pos(), variant.alt().charAt(0));
        }
        else {
            int[] window = getIndelWindow(variant);
            alt = evidence.getIndelCount(variant.contig(), window[0], window[1], variant.isInsertion() ? INSERTION : DELETION, variant.length());
        }
        return new int[]{alt, evidence.getDepth(variant.contig(), variant.pos())};
    }

    public Table evaluateSmallVariants() {
        Table report = new Table("id", "type", "origin", "level", "zygosity", "contig", "pos", "ref", "alt", "regions",
                "expected_normal_vaf", "expected_tumor_vaf",
                "normal_in_alt", "normal_in_depth", "normal_in_vaf", "tumor_in_alt", "tumor_in_depth", "tumor_in_vaf",
                "normal_out_alt", "normal_out_depth", "normal_out_vaf", "tumor_out_alt", "tumor_out_depth", "tumor_out_vaf",
                "expectation", "verdict");
        for (TruthVariant variant : dataset.getTruthVariants()) {
            if (variant.isStructural()) continue;
            int[] nIn = countAllele(normalIn, variant);
            int[] tIn = countAllele(tumorIn, variant);
            int[] nOut = countAllele(normalOut, variant);
            int[] tOut = countAllele(tumorOut, variant);
            String expectation;
            if (nIn[0] > 0) expectation = EXPECT_REMOVE_SEEN_IN_NORMAL;
            else if (tIn[0] == 0) expectation = EXPECT_NOT_SAMPLED;
            else if (nIn[1] == 0 || tIn[1] == 0) expectation = EXPECT_REMOVE_UNCOVERED;
            else expectation = EXPECT_KEEP;
            String verdict = VERDICT_OK;
            if (EXPECT_REMOVE_SEEN_IN_NORMAL.equals(expectation) || EXPECT_REMOVE_UNCOVERED.equals(expectation)) {
                if (nOut[0] + tOut[0] > 0) verdict = VERDICT_LEAK;
            }
            else if (EXPECT_NOT_SAMPLED.equals(expectation)) {
                if (nOut[0] + tOut[0] > 0) verdict = VERDICT_ARTIFACT;
            }
            else {
                if (nOut[0] > 0) verdict = VERDICT_ARTIFACT;
                else if (tOut[0] == 0) verdict = VERDICT_SOMATIC_LOST;
                else if (tOut[0] < tIn[0]) verdict = VERDICT_SOMATIC_REDUCED;
                else if (tOut[0] > tIn[0]) verdict = VERDICT_ALT_INCREASED;
                else if (variant.isGermline()) verdict = VERDICT_GERMLINE_UNSEEN_KEPT;
            }
            report.addRow(variant.id(), variant.type(), variant.origin(), variant.level(), variant.zygosity(), variant.contig(),
                    variant.pos(), variant.ref(), variant.alt(), variant.regions(),
                    String.format(Locale.ROOT, "%.3f", variant.expectedNormalVaf()), String.format(Locale.ROOT, "%.3f", variant.expectedTumorVaf()),
                    nIn[0], nIn[1], formatVaf(nIn[0], nIn[1]), tIn[0], tIn[1], formatVaf(tIn[0], tIn[1]),
                    nOut[0], nOut[1], formatVaf(nOut[0], nOut[1]), tOut[0], tOut[1], formatVaf(tOut[0], tOut[1]),
                    expectation, verdict);
        }
        return report;
    }

    // ---------------------------------------------------------------------------------------------------------
    // Spiked-in structural variants
    // ---------------------------------------------------------------------------------------------------------

    private List<Object[]> getBreakpoints(TruthVariant variant) {
        List<Object[]> breakpoints = new ArrayList<>();
        breakpoints.add(new Object[]{variant.contig(), variant.pos()});
        if ("SV_TRA".equals(variant.type())) breakpoints.add(new Object[]{variant.mateContig(), variant.matePos()});
        else if (variant.end() != variant.pos()) breakpoints.add(new Object[]{variant.contig(), variant.end()});
        return breakpoints;
    }

    private boolean isNear(SAMRecord record, List<Object[]> breakpoints, int window) {
        for (Object[] breakpoint : breakpoints) {
            int pos = (Integer) breakpoint[1];
            if (record.getContig().equals(breakpoint[0]) && record.getAlignmentStart() <= pos + window
                    && record.getAlignmentEnd() >= pos - window) return true;
        }
        return false;
    }

    private boolean hasPairSignal(SAMRecord record, String type) {
        if (record.isSecondaryOrSupplementary() || !record.getReadPairedFlag() || record.getMateUnmappedFlag()) return false;
        boolean sameContig = record.getReferenceIndex().equals(record.getMateReferenceIndex());
        if ("SV_TRA".equals(type)) return !sameContig;
        if (!sameContig) return false;
        int insertSize = Math.abs(record.getInferredInsertSize());
        boolean sameStrand = record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag();
        if ("SV_DEL".equals(type)) return !sameStrand && insertSize >= LARGE_INSERT_SIZE;
        if ("SV_INS".equals(type)) return !sameStrand && insertSize > 0 && insertSize <= SMALL_INSERT_SIZE;
        if ("SV_INV".equals(type)) return sameStrand;
        if ("SV_DUP".equals(type)) {
            // Everted pair: the read on the forward strand is the rightmost one
            if (sameStrand) return false;
            boolean forwardIsRightmost = record.getReadNegativeStrandFlag()
                    ? record.getMateAlignmentStart() > record.getAlignmentStart()
                    : record.getAlignmentStart() > record.getMateAlignmentStart();
            return forwardIsRightmost;
        }
        return false;
    }

    private boolean hasClipSignal(SAMRecord record, List<Object[]> breakpoints) {
        List<CigarElement> cigarElements = record.getCigar().getCigarElements();
        if (cigarElements.isEmpty()) return false;
        CigarElement first = cigarElements.get(0);
        CigarElement last = cigarElements.get(cigarElements.size() - 1);
        for (Object[] breakpoint : breakpoints) {
            if (!record.getContig().equals(breakpoint[0])) continue;
            int pos = (Integer) breakpoint[1];
            if (isClip(first) && Math.abs(record.getAlignmentStart() - pos) <= CLIP_SIGNAL_WINDOW) return true;
            if (isClip(last) && Math.abs(record.getAlignmentEnd() - pos) <= CLIP_SIGNAL_WINDOW) return true;
        }
        return false;
    }

    private boolean isClip(CigarElement cigarElement) {
        CigarOperator op = cigarElement.getOperator();
        return (CigarOperator.S == op || CigarOperator.H == op) && cigarElement.getLength() >= MIN_CLIP_LENGTH;
    }

    /** @return number of reads with abnormal pair signal, with clipping signal and supplementary alignments */
    private int[] countStructuralSignals(AlignmentEvidence evidence, TruthVariant variant) {
        List<Object[]> breakpoints = getBreakpoints(variant);
        int[] answer = new int[3];
        for (SAMRecord record : evidence.getRecords()) {
            if (isNear(record, breakpoints, PAIR_SIGNAL_WINDOW) && hasPairSignal(record, variant.type())) answer[0]++;
            if (hasClipSignal(record, breakpoints)) answer[1]++;
            if (record.getSupplementaryAlignmentFlag() && isNear(record, breakpoints, SUPPLEMENTARY_SIGNAL_WINDOW)) answer[2]++;
        }
        return answer;
    }

    public Table evaluateStructuralVariants() {
        Table report = new Table("id", "type", "origin", "level", "zygosity", "contig", "pos", "end", "mate_contig", "mate_pos",
                "normal_in_pairs", "normal_out_pairs", "tumor_in_pairs", "tumor_out_pairs",
                "normal_in_clips", "normal_out_clips", "tumor_in_clips", "tumor_out_clips",
                "normal_in_supplementary", "normal_out_supplementary", "tumor_in_supplementary", "tumor_out_supplementary",
                "expectation", "verdict");
        for (TruthVariant variant : dataset.getTruthVariants()) {
            if (!variant.isStructural()) continue;
            int[] nIn = countStructuralSignals(normalIn, variant);
            int[] tIn = countStructuralSignals(tumorIn, variant);
            int[] nOut = countStructuralSignals(normalOut, variant);
            int[] tOut = countStructuralSignals(tumorOut, variant);
            int normalInSignal = nIn[0] + nIn[1];
            int tumorInSignal = tIn[0] + tIn[1];
            String expectation;
            if (normalInSignal > 0 && variant.isGermline()) expectation = EXPECT_REMOVE_SEEN_IN_NORMAL;
            else if (normalInSignal == 0 && tumorInSignal > 0 && !variant.isGermline()) expectation = EXPECT_KEEP;
            else if (normalInSignal + tumorInSignal == 0) expectation = EXPECT_NOT_SAMPLED;
            else expectation = EXPECT_UNDETERMINED;
            String verdict = VERDICT_OK;
            if (EXPECT_REMOVE_SEEN_IN_NORMAL.equals(expectation)) {
                if (nOut[0] + nOut[1] + nOut[2] + tOut[0] + tOut[1] + tOut[2] > 0) verdict = VERDICT_RESIDUAL_SIGNAL;
            }
            else if (EXPECT_KEEP.equals(expectation)) {
                if (tOut[0] + tOut[1] == 0) verdict = VERDICT_SOMATIC_LOST;
                else if (tOut[0] < tIn[0] || tOut[1] < tIn[1]) verdict = VERDICT_SOMATIC_REDUCED;
            }
            report.addRow(variant.id(), variant.type(), variant.origin(), variant.level(), variant.zygosity(), variant.contig(),
                    variant.pos(), variant.end(), variant.mateContig(), variant.matePos() < 0 ? "." : variant.matePos(),
                    nIn[0], nOut[0], tIn[0], tOut[0], nIn[1], nOut[1], tIn[1], tOut[1], nIn[2], nOut[2], tIn[2], tOut[2],
                    expectation, verdict);
        }
        return report;
    }

    // ---------------------------------------------------------------------------------------------------------
    // Comparison against a baseline
    // ---------------------------------------------------------------------------------------------------------

    /**
     * Lists every difference between two tables in the given columns, including rows present in only one of them
     */
    public static Table compareToBaseline(Table current, Table baseline, String... columns) {
        Table differences = new Table("key", "id", "column", "baseline", "current");
        Set<String> keys = new java.util.LinkedHashSet<>(baseline.getRows().keySet());
        keys.addAll(current.getRows().keySet());
        for (String key : keys) {
            Map<String, String> currentRow = current.getRow(key);
            Map<String, String> baselineRow = baseline.getRow(key);
            if (currentRow == null || baselineRow == null) {
                differences.addRow(key + "|row", key, "row", baselineRow == null ? "<absent>" : "<present>",
                        currentRow == null ? "<absent>" : "<present>");
                continue;
            }
            for (String column : columns) {
                String currentValue = currentRow.get(column);
                String baselineValue = baselineRow.get(column);
                if (currentValue == null || !currentValue.equals(baselineValue)) {
                    differences.addRow(key + "|" + column, key, column, baselineValue, currentValue);
                }
            }
        }
        return differences;
    }

    /** Keys of the rows of a table that have a given value in a column */
    public static Set<String> keysWith(Table table, String column, String value) {
        Set<String> answer = new HashSet<>();
        for (Map.Entry<String, Map<String, String>> entry : table.getRows().entrySet()) {
            if (value.equals(entry.getValue().get(column))) answer.add(entry.getKey());
        }
        return answer;
    }
}
