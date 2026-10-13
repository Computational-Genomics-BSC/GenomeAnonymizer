package simulated;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Everything the validation needs to know about one alignment file: base counts per reference position, indels
 * found in the CIGARs and the records themselves. The simulated files are small, so all is kept in memory
 * @author Nicolas Gaitan
 */
public class AlignmentEvidence {

    public static final int BASE_COUNT = 5;
    public static final String BASES = "ACGTN";
    public static final char INSERTION = 'I';
    public static final char DELETION = 'D';

    private final SAMFileHeader header;
    // contig -> [0-based position][base index]
    private final Map<String, int[][]> baseCounts = new HashMap<>();
    // contig:anchor:type:length -> number of reads, the anchor is the 1-based position of the base before the indel
    private final Map<String, Integer> indelCounts = new HashMap<>();
    private final List<SAMRecord> records = new ArrayList<>();

    private AlignmentEvidence(SAMFileHeader header, Map<String, byte[]> reference) {
        this.header = header;
        for (Map.Entry<String, byte[]> entry : reference.entrySet()) {
            baseCounts.put(entry.getKey(), new int[entry.getValue().length][BASE_COUNT]);
        }
    }

    /**
     * @param onlyEligibleReads if true, only the reads that GenomeAnonymizer takes into account by default are
     *                          loaded: mapped, with the mate mapped and not marked as duplicates
     */
    public static AlignmentEvidence load(File file, Map<String, byte[]> reference, boolean onlyEligibleReads) throws IOException {
        SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        try (SamReader reader = factory.open(file)) {
            AlignmentEvidence answer = new AlignmentEvidence(reader.getFileHeader(), reference);
            for (SAMRecord record : reader) {
                if (record.getReadUnmappedFlag()) continue;
                if (onlyEligibleReads && !isEligible(record)) continue;
                answer.addRecord(record);
            }
            return answer;
        }
    }

    public static boolean isEligible(SAMRecord record) {
        if (record.getReadUnmappedFlag() || record.getDuplicateReadFlag()) return false;
        return record.getReadPairedFlag() && !record.getMateUnmappedFlag();
    }

    private void addRecord(SAMRecord record) {
        records.add(record);
        int[][] contigCounts = baseCounts.get(record.getContig());
        byte[] bases = record.getReadBases();
        int refPos = record.getAlignmentStart();
        int readPos = 0;
        for (CigarElement cigarElement : record.getCigar().getCigarElements()) {
            CigarOperator op = cigarElement.getOperator();
            int length = cigarElement.getLength();
            if (CigarOperator.M == op || CigarOperator.EQ == op || CigarOperator.X == op) {
                for (int i = 0; i < length; i++) {
                    // Bases aligned outside the reference are ignored, they are reported by the record validation
                    if (refPos + i - 1 < contigCounts.length) contigCounts[refPos + i - 1][baseIndex(bases[readPos + i])]++;
                }
            }
            else if (CigarOperator.I == op) {
                indelCounts.merge(indelKey(record.getContig(), refPos - 1, INSERTION, length), 1, Integer::sum);
            }
            else if (CigarOperator.D == op) {
                indelCounts.merge(indelKey(record.getContig(), refPos - 1, DELETION, length), 1, Integer::sum);
            }
            if (op.consumesReferenceBases()) refPos += length;
            if (op.consumesReadBases()) readPos += length;
        }
    }

    public static int baseIndex(byte base) {
        int answer = BASES.indexOf(Character.toUpperCase((char) base));
        return answer < 0 ? BASE_COUNT - 1 : answer;
    }

    public static String indelKey(String contig, int anchor, char type, int length) {
        return contig + ":" + anchor + ":" + type + ":" + length;
    }

    public SAMFileHeader getHeader() {
        return header;
    }

    public List<SAMRecord> getRecords() {
        return records;
    }

    public Map<String, Integer> getIndelCounts() {
        return indelCounts;
    }

    /** @param pos 1-based position */
    public int getDepth(String contig, int pos) {
        int[][] contigCounts = baseCounts.get(contig);
        if (pos < 1 || pos > contigCounts.length) return 0;
        int answer = 0;
        for (int count : contigCounts[pos - 1]) answer += count;
        return answer;
    }

    /** @param pos 1-based position */
    public int getBaseCount(String contig, int pos, char base) {
        int[][] contigCounts = baseCounts.get(contig);
        if (pos < 1 || pos > contigCounts.length) return 0;
        return contigCounts[pos - 1][baseIndex((byte) base)];
    }

    public int getIndelCount(String contig, int anchor, char type, int length) {
        return indelCounts.getOrDefault(indelKey(contig, anchor, type, length), 0);
    }

    /** Number of reads with an indel of the given type and length anchored in the interval [first, last] */
    public int getIndelCount(String contig, int first, int last, char type, int length) {
        int answer = 0;
        for (int anchor = first; anchor <= last; anchor++) {
            answer += getIndelCount(contig, anchor, type, length);
        }
        return answer;
    }

    /** Sorted depths of every position of the genome */
    public int[] getSortedDepths() {
        List<Integer> depths = new ArrayList<>();
        for (Map.Entry<String, int[][]> entry : baseCounts.entrySet()) {
            for (int pos = 1; pos <= entry.getValue().length; pos++) {
                depths.add(getDepth(entry.getKey(), pos));
            }
        }
        int[] answer = new int[depths.size()];
        for (int i = 0; i < answer.length; i++) answer[i] = depths.get(i);
        java.util.Arrays.sort(answer);
        return answer;
    }
}
