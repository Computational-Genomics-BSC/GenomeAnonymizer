package simulated;

import htsjdk.samtools.reference.FastaSequenceIndexEntry;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * A simulated tumor-normal pair with its reference, the truth table of spiked-in variants and its special regions,
 * as written by src/test/simulation/build_datasets.sh
 * @author Nicolas Gaitan
 */
public class SimulatedDataset {

    public static final String NORMAL_FILE = "normal.bam";
    public static final String TUMOR_FILE = "tumor.bam";
    public static final String REFERENCE_FILE = "reference.fa";
    public static final String TRUTH_FILE = "truth_variants.tsv";
    public static final String REGIONS_FILE = "regions.bed";

    public static final String ORIGIN_GERMLINE = "germline";
    public static final String ORIGIN_SOMATIC = "somatic";
    public static final String REGION_KIND_COVERAGE = "coverage";
    public static final String REGION_KIND_REFERENCE = "reference";

    /**
     * pos and end are 1-based and inclusive. For indels pos is the base before the event, as in VCF. For structural
     * variants pos and end are the last bases before each breakpoint
     */
    public record TruthVariant(String id, String type, String origin, String level, String zygosity, String contig,
                               int pos, int end, String ref, String alt, int length, String mateContig, int matePos,
                               double expectedNormalVaf, double expectedTumorVaf, String regions) {
        public boolean isGermline() {
            return ORIGIN_GERMLINE.equals(origin);
        }

        public boolean isStructural() {
            return type.startsWith("SV_");
        }

        public boolean isSNV() {
            return "SNV".equals(type);
        }

        public boolean isInsertion() {
            return "INS".equals(type);
        }
    }

    /** start is 0-based and end is exclusive, as in BED */
    public record Region(String contig, int start, int end, String name, String kind, String details) {
        public boolean overlaps(String otherContig, int first, int last) {
            return contig.equals(otherContig) && first <= end && last >= start + 1;
        }
    }

    private final File directory;
    // Reference sequences in upper case, 0-based arrays
    private final Map<String, byte[]> reference = new LinkedHashMap<>();
    private final List<TruthVariant> truthVariants = new ArrayList<>();
    private final List<Region> regions = new ArrayList<>();

    public SimulatedDataset(File directory) throws IOException {
        this.directory = directory;
        loadReference();
        loadTruth();
        loadRegions();
    }

    private void loadReference() throws IOException {
        try (IndexedFastaSequenceFile fasta = new IndexedFastaSequenceFile(getReferenceFile())) {
            for (FastaSequenceIndexEntry entry : fasta.getIndex()) {
                byte[] bases = fasta.getSequence(entry.getContig()).getBases();
                for (int i = 0; i < bases.length; i++) {
                    bases[i] = (byte) Character.toUpperCase(bases[i]);
                }
                reference.put(entry.getContig(), bases);
            }
        }
    }

    private void loadTruth() throws IOException {
        Table table = Table.read(new File(directory, TRUTH_FILE));
        for (Map<String, String> row : table.getRows().values()) {
            String matePos = row.get("mate_pos");
            truthVariants.add(new TruthVariant(row.get("id"), row.get("type"), row.get("origin"), row.get("level"),
                    row.get("zygosity"), row.get("contig"), Integer.parseInt(row.get("pos")), Integer.parseInt(row.get("end")),
                    row.get("ref"), row.get("alt"), Integer.parseInt(row.get("length")), row.get("mate_contig"),
                    ".".equals(matePos) ? -1 : Integer.parseInt(matePos), Double.parseDouble(row.get("expected_normal_vaf")),
                    Double.parseDouble(row.get("expected_tumor_vaf")), row.get("regions")));
        }
    }

    private void loadRegions() throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(new File(directory, REGIONS_FILE)))) {
            String line = reader.readLine();
            while (line != null) {
                String[] elements = line.split("\t");
                regions.add(new Region(elements[0], Integer.parseInt(elements[1]), Integer.parseInt(elements[2]),
                        elements[3], elements[4], elements[5]));
                line = reader.readLine();
            }
        }
    }

    public File getNormalFile() {
        return new File(directory, NORMAL_FILE);
    }

    public File getTumorFile() {
        return new File(directory, TUMOR_FILE);
    }

    public File getReferenceFile() {
        return new File(directory, REFERENCE_FILE);
    }

    public Map<String, byte[]> getReference() {
        return reference;
    }

    public List<TruthVariant> getTruthVariants() {
        return truthVariants;
    }

    public List<Region> getRegions() {
        return regions;
    }

    public Region getRegion(String name) {
        for (Region region : regions) {
            if (region.name().equals(name)) return region;
        }
        throw new IllegalArgumentException("Region not found in dataset: " + name);
    }

    public long getGenomeSize() {
        long answer = 0;
        for (byte[] sequence : reference.values()) answer += sequence.length;
        return answer;
    }
}
