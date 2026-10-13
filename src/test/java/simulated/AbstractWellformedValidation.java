package simulated;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.SAMValidationError;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.SequenceUtil;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInstance;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Pattern;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Validation of a complete run of the GenomeAnonymizer JAR over the well-formed simulated tumor-normal pair.
 *
 * The run is launched once per subclass (one per sample type mode). The checks are of three kinds:
 *  - Hard requirements of any output: valid files, identifiers anonymized, reproducibility.
 *  - Privacy: germline alleles found in the outputs and artifacts. Any finding that the baseline does not have
 *    makes the test fail.
 *  - Recorded behaviour: the fate and the VAFs of every spiked-in variant are written to the reports and compared
 *    with the baseline, which is the result of the build that was accepted last. Differences are listed in
 *    baseline_differences.tsv and only make the test fail with -PstrictBaseline. They are accepted with
 *    -PupdateBaseline.
 * @author Nicolas Gaitan
 */
@Tag(SimulatedTestConfig.SIMULATED_TAG)
@TestInstance(TestInstance.Lifecycle.PER_CLASS)
public abstract class AbstractWellformedValidation {

    public static final String GENOME_SCAN_FILE = "genome_scan.tsv";
    public static final String SMALL_VARIANTS_FILE = "small_variants.tsv";
    public static final String STRUCTURAL_VARIANTS_FILE = "structural_variants.tsv";
    public static final String BASELINE_DIFFERENCES_FILE = "baseline_differences.tsv";
    public static final String SUMMARY_FILE = "summary.txt";
    public static final String KNOWN_FAILURES_FILE = "known_failures.tsv";

    public static final String[] SMALL_VARIANT_TRACKED_COLUMNS = {"normal_in_alt", "normal_in_depth", "tumor_in_alt",
            "tumor_in_depth", "normal_out_alt", "normal_out_depth", "tumor_out_alt", "tumor_out_depth", "expectation", "verdict"};
    public static final String[] STRUCTURAL_VARIANT_TRACKED_COLUMNS = {"normal_in_pairs", "normal_out_pairs", "tumor_in_pairs",
            "tumor_out_pairs", "normal_in_clips", "normal_out_clips", "tumor_in_clips", "tumor_out_clips",
            "normal_out_supplementary", "tumor_out_supplementary", "expectation", "verdict"};
    public static final String[] GENOME_SCAN_TRACKED_COLUMNS = {"category", "reason", "normal_out_alt", "tumor_out_alt"};

    // The whole simulated validation has to stay well below this, so a single run is given a fraction of it
    public static final double MAX_RUN_SECONDS = 300;
    public static final double MIN_KEPT_ALIGNMENTS_FRACTION = 0.9;
    private static final Pattern ANONYMIZED_READ_NAME = Pattern.compile("[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}");

    protected SimulatedDataset dataset;
    protected AnonymizerRunner.RunResult run;
    protected File normalOutputFile;
    protected File tumorOutputFile;
    protected AlignmentEvidence normalInAll;
    protected AlignmentEvidence tumorInAll;
    protected AlignmentEvidence normalOut;
    protected AlignmentEvidence tumorOut;

    protected Table genomeScan;
    protected Table smallVariants;
    protected Table structuralVariants;
    protected Table baselineGenomeScan;
    protected Table baselineSmallVariants;
    protected Table baselineStructuralVariants;
    protected Table baselineDifferences;
    // Requirement id -> number of problems found in this run, and in the run of the baseline
    protected Map<String, Integer> requirementProblems = new TreeMap<>();
    protected Map<String, Integer> knownFailures = new TreeMap<>();

    /** Name of the mode, used for the directories of outputs, reports and baselines */
    protected abstract String getMode();

    /** Arguments added to the common ones (inputs, reference, output prefix, threads and seed) */
    protected abstract String[] getExtraArgs();

    protected File getModeOutputDir() {
        return new File(SimulatedTestConfig.getOutputDir(), getMode());
    }

    protected File getModeReportDir() {
        return new File(SimulatedTestConfig.getReportDir(), getMode());
    }

    protected File getModeBaselineDir() {
        return new File(SimulatedTestConfig.getBaselineDir(), getMode());
    }

    @BeforeAll
    void runAnonymizerAndEvaluate() throws Exception {
        dataset = new SimulatedDataset(SimulatedTestConfig.getWellformedDir());
        String outputPrefix = new File(getModeOutputDir(), "simulated").getPath();
        run = AnonymizerRunner.run(dataset.getNormalFile(), dataset.getTumorFile(), dataset.getReferenceFile(), outputPrefix, getExtraArgs());
        assertEquals(0, run.exitCode(), "GenomeAnonymizer failed on the well-formed dataset. Log: " + run.logFile());
        // With VAF correction the final product is the corrected pair of files, when they are written
        normalOutputFile = run.normalFixVAFOutput().exists() ? run.normalFixVAFOutput() : run.normalOutput();
        tumorOutputFile = run.tumorFixVAFOutput().exists() ? run.tumorFixVAFOutput() : run.tumorOutput();
        assertTrue(normalOutputFile.exists(), "Missing output " + normalOutputFile);
        assertTrue(tumorOutputFile.exists(), "Missing output " + tumorOutputFile);

        Map<String, byte[]> reference = dataset.getReference();
        AlignmentEvidence normalIn = AlignmentEvidence.load(dataset.getNormalFile(), reference, true);
        AlignmentEvidence tumorIn = AlignmentEvidence.load(dataset.getTumorFile(), reference, true);
        normalInAll = AlignmentEvidence.load(dataset.getNormalFile(), reference, false);
        tumorInAll = AlignmentEvidence.load(dataset.getTumorFile(), reference, false);
        normalOut = AlignmentEvidence.load(normalOutputFile, reference, false);
        tumorOut = AlignmentEvidence.load(tumorOutputFile, reference, false);
        OutputValidator validator = new OutputValidator(dataset, normalIn, tumorIn, normalInAll, tumorInAll, normalOut, tumorOut);
        genomeScan = validator.scanGenome();
        smallVariants = validator.evaluateSmallVariants();
        structuralVariants = validator.evaluateStructuralVariants();

        if (SimulatedTestConfig.updateBaseline()) {
            genomeScan.write(new File(getModeBaselineDir(), GENOME_SCAN_FILE));
            smallVariants.write(new File(getModeBaselineDir(), SMALL_VARIANTS_FILE));
            structuralVariants.write(new File(getModeBaselineDir(), STRUCTURAL_VARIANTS_FILE));
            System.out.println("[" + getMode() + "] Baseline updated in " + getModeBaselineDir());
        }
        File baselineFile = new File(getModeBaselineDir(), SMALL_VARIANTS_FILE);
        assertTrue(baselineFile.exists(), "There is no baseline for mode " + getMode() + " in " + getModeBaselineDir()
                + ". Create it with ./gradlew simulatedTest -PupdateBaseline");
        knownFailures = readKnownFailures();
        baselineGenomeScan = Table.read(new File(getModeBaselineDir(), GENOME_SCAN_FILE));
        baselineSmallVariants = Table.read(baselineFile);
        baselineStructuralVariants = Table.read(new File(getModeBaselineDir(), STRUCTURAL_VARIANTS_FILE));
        baselineDifferences = new Table("key", "table", "id", "column", "baseline", "current");
        addDifferences("genome_scan", OutputValidator.compareToBaseline(genomeScan, baselineGenomeScan, GENOME_SCAN_TRACKED_COLUMNS));
        addDifferences("small_variants", OutputValidator.compareToBaseline(smallVariants, baselineSmallVariants, SMALL_VARIANT_TRACKED_COLUMNS));
        addDifferences("structural_variants", OutputValidator.compareToBaseline(structuralVariants, baselineStructuralVariants, STRUCTURAL_VARIANT_TRACKED_COLUMNS));
        writeReports();
    }

    private void addDifferences(String tableName, Table differences) {
        for (Map<String, String> row : differences.getRows().values()) {
            baselineDifferences.addRow(tableName + "|" + row.get("key"), tableName, row.get("id"), row.get("column"),
                    row.get("baseline"), row.get("current"));
        }
    }

    // ---------------------------------------------------------------------------------------------------------
    // Requirements of any output
    // ---------------------------------------------------------------------------------------------------------

    /**
     * Registers the result of a requirement. A requirement with problems makes the test fail, unless the baseline
     * already had at least as many problems for it: those are known failures of the accepted build, and they are
     * printed as warnings and listed in the summary.
     */
    private void checkRequirement(String id, List<String> problems) {
        requirementProblems.put(id, problems.size());
        if (problems.isEmpty()) return;
        String message = "Requirement " + id + " has " + problems.size() + " problems in mode " + getMode() + ", first ones:\n"
                + String.join("\n", problems.subList(0, Math.min(10, problems.size())));
        int knownProblems = knownFailures.getOrDefault(id, 0);
        if (SimulatedTestConfig.updateBaseline() || problems.size() <= knownProblems) {
            System.out.println("[" + getMode() + "] KNOWN FAILURE: " + message);
            return;
        }
        if (knownProblems > 0) message = "(the baseline had " + knownProblems + " problems) " + message;
        fail(message);
    }

    private Map<String, Integer> readKnownFailures() throws IOException {
        Map<String, Integer> answer = new TreeMap<>();
        File file = new File(getModeBaselineDir(), KNOWN_FAILURES_FILE);
        if (!file.exists()) return answer;
        for (Map<String, String> row : Table.read(file).getRows().values()) {
            answer.put(row.get("requirement"), Integer.parseInt(row.get("problems")));
        }
        return answer;
    }

    @Test
    void runFinishesInTime() {
        List<String> problems = new ArrayList<>();
        if (run.seconds() >= MAX_RUN_SECONDS) problems.add("The run took " + run.seconds() + " seconds, the limit is " + MAX_RUN_SECONDS);
        checkRequirement("run_finishes_in_time", problems);
    }

    @Test
    void outputsAreSortedAndIndexed() throws IOException {
        List<String> problems = new ArrayList<>();
        for (File file : new File[]{normalOutputFile, tumorOutputFile}) {
            SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
            try (SamReader reader = factory.open(file)) {
                SAMFileHeader header = reader.getFileHeader();
                if (SAMFileHeader.SortOrder.coordinate != header.getSortOrder()) problems.add("Sort order " + header.getSortOrder() + " in " + file.getName());
                if (!reader.hasIndex()) problems.add("Output without index: " + file.getName());
                List<String> contigs = new ArrayList<>();
                for (SAMSequenceRecord sequence : header.getSequenceDictionary().getSequences()) contigs.add(sequence.getSequenceName());
                if (!contigs.equals(new ArrayList<>(dataset.getReference().keySet()))) problems.add("Contigs " + contigs + " in " + file.getName());
                int previousIdx = -1;
                int previousStart = -1;
                for (SAMRecord record : reader) {
                    int idx = record.getReferenceIndex();
                    if (idx < previousIdx || (idx == previousIdx && record.getAlignmentStart() < previousStart)) {
                        problems.add("Record out of order in " + file.getName() + ": " + record.getSAMString().trim());
                    }
                    previousIdx = idx;
                    previousStart = record.getAlignmentStart();
                }
            }
        }
        checkRequirement("outputs_sorted_and_indexed", problems);
    }

    @Test
    void recordsAreValid() {
        List<String> problems = new ArrayList<>();
        for (AlignmentEvidence output : new AlignmentEvidence[]{normalOut, tumorOut}) {
            for (SAMRecord record : output.getRecords()) {
                List<SAMValidationError> errors = record.isValid();
                if (errors != null) {
                    for (SAMValidationError error : errors) problems.add(error.getMessage() + " in " + record.getSAMString().trim());
                }
                if (record.getReadLength() != record.getCigar().getReadLength()) {
                    problems.add("CIGAR and sequence lengths differ in " + record.getSAMString().trim());
                }
                if (record.getReadLength() != record.getBaseQualities().length) {
                    problems.add("Sequence and qualities lengths differ in " + record.getSAMString().trim());
                }
                if (record.getAlignmentEnd() > dataset.getReference().get(record.getContig()).length) {
                    problems.add("Record aligned beyond the end of the contig: " + record.getSAMString().trim());
                }
            }
        }
        checkRequirement("records_valid", problems);
    }

    @Test
    void editDistanceTagsMatchReference() {
        List<String> problems = new ArrayList<>();
        for (AlignmentEvidence output : new AlignmentEvidence[]{normalOut, tumorOut}) {
            for (SAMRecord record : output.getRecords()) {
                Integer editDistance = record.getIntegerAttribute("NM");
                if (editDistance == null) continue;
                byte[] contigBases = dataset.getReference().get(record.getContig());
                int expected = SequenceUtil.calculateSamNmTag(record, contigBases);
                if (expected != editDistance) problems.add("NM is " + editDistance + " and should be " + expected + " in " + record.getSAMString().trim());
            }
        }
        checkRequirement("edit_distance_tags_match_reference", problems);
    }

    @Test
    void matesAreConsistent() {
        List<String> problems = new ArrayList<>();
        for (AlignmentEvidence output : new AlignmentEvidence[]{normalOut, tumorOut}) {
            Map<String, SAMRecord> primaryRecords = new HashMap<>();
            for (SAMRecord record : output.getRecords()) {
                if (record.isSecondaryOrSupplementary() || !record.getReadPairedFlag()) continue;
                primaryRecords.put(record.getReadName() + "/" + (record.getFirstOfPairFlag() ? 1 : 2), record);
            }
            for (SAMRecord record : output.getRecords()) {
                if (record.isSecondaryOrSupplementary() || !record.getReadPairedFlag() || record.getMateUnmappedFlag()) continue;
                SAMRecord mate = primaryRecords.get(record.getReadName() + "/" + (record.getFirstOfPairFlag() ? 2 : 1));
                if (mate == null) {
                    problems.add("Mate not found for " + record.getSAMString().trim());
                    continue;
                }
                if (!record.getMateReferenceName().equals(mate.getContig()) || record.getMateAlignmentStart() != mate.getAlignmentStart()) {
                    problems.add("Mate position is " + mate.getContig() + ":" + mate.getAlignmentStart() + " for " + record.getSAMString().trim());
                }
                else if (record.getMateNegativeStrandFlag() != mate.getReadNegativeStrandFlag()) {
                    problems.add("Mate strand flag is wrong in " + record.getSAMString().trim());
                }
            }
        }
        checkRequirement("mates_consistent", problems);
    }

    @Test
    void readNamesAreAnonymized() {
        List<String> problems = new ArrayList<>();
        AlignmentEvidence[][] samples = {{normalInAll, normalOut}, {tumorInAll, tumorOut}};
        for (AlignmentEvidence[] sample : samples) {
            Set<String> originalNames = new HashSet<>();
            for (SAMRecord record : sample[0].getRecords()) originalNames.add(record.getReadName());
            for (SAMRecord record : sample[1].getRecords()) {
                String name = record.getReadName();
                if (originalNames.contains(name)) problems.add("Original read name in the output: " + name);
                else if (!ANONYMIZED_READ_NAME.matcher(name).matches()) problems.add("Unexpected read name in the output: " + name);
            }
        }
        checkRequirement("read_names_anonymized", problems);
    }

    private Set<String> getOriginalIdentifiers(AlignmentEvidence input) {
        Set<String> answer = new TreeSet<>();
        for (SAMReadGroupRecord readGroup : input.getHeader().getReadGroups()) {
            answer.add(readGroup.getId());
            if (readGroup.getSample() != null) answer.add(readGroup.getSample());
            if (readGroup.getLibrary() != null) answer.add(readGroup.getLibrary());
            if (readGroup.getPlatformUnit() != null) answer.add(readGroup.getPlatformUnit());
        }
        return answer;
    }

    @Test
    void readGroupsAreAnonymized() {
        List<String> problems = new ArrayList<>();
        AlignmentEvidence[][] samples = {{normalInAll, normalOut}, {tumorInAll, tumorOut}};
        for (AlignmentEvidence[] sample : samples) {
            List<SAMReadGroupRecord> outputReadGroups = sample[1].getHeader().getReadGroups();
            if (outputReadGroups.size() != 1) {
                problems.add("Expected one read group in the output and found " + outputReadGroups.size());
                continue;
            }
            SAMReadGroupRecord outputReadGroup = outputReadGroups.get(0);
            String readGroupLine = outputReadGroup.getSAMString();
            for (String original : getOriginalIdentifiers(sample[0])) {
                if (readGroupLine.contains(original)) problems.add("Original identifier '" + original + "' in the read group of the output: " + readGroupLine);
            }
            for (SAMRecord record : sample[1].getRecords()) {
                if (!outputReadGroup.getId().equals(record.getStringAttribute("RG"))) problems.add("Read group " + record.getStringAttribute("RG") + " in " + record.getReadName());
            }
        }
        checkRequirement("read_groups_anonymized", problems);
    }

    @Test
    void headerDoesNotContainOriginalIdentifiers() {
        List<String> problems = new ArrayList<>();
        AlignmentEvidence[][] samples = {{normalInAll, normalOut}, {tumorInAll, tumorOut}};
        for (AlignmentEvidence[] sample : samples) {
            StringWriter headerText = new StringWriter();
            new SAMTextHeaderCodec().encode(headerText, sample[1].getHeader());
            for (String line : headerText.toString().split("\n")) {
                for (String original : getOriginalIdentifiers(sample[0])) {
                    if (line.contains(original)) problems.add("Original identifier '" + original + "' in the header of the output: " + line);
                }
            }
        }
        checkRequirement("header_without_original_identifiers", problems);
    }

    @Test
    void temporaryTagsAndFilesAreRemoved() {
        List<String> problems = new ArrayList<>();
        for (AlignmentEvidence output : new AlignmentEvidence[]{normalOut, tumorOut}) {
            for (SAMRecord record : output.getRecords()) {
                if (record.getAttribute("op") != null) problems.add("Temporary tag op left in " + record.getReadName());
                if (record.getAttribute("MC") != null) problems.add("Tag MC left in " + record.getReadName());
            }
        }
        String prefixName = new File(run.outputPrefix()).getName();
        File[] partitionFiles = getModeOutputDir().listFiles((dir, name) -> name.startsWith(prefixName + "_"));
        if (partitionFiles != null) {
            for (File file : partitionFiles) problems.add("Partition file left: " + file.getName());
        }
        checkRequirement("temporary_tags_and_files_removed", problems);
    }

    @Test
    void outputsKeepMostReadsAndDoNotGainAny() {
        List<String> problems = new ArrayList<>();
        AlignmentEvidence[][] samples = {{normalInAll, normalOut}, {tumorInAll, tumorOut}};
        for (AlignmentEvidence[] sample : samples) {
            int in = sample[0].getRecords().size();
            int out = sample[1].getRecords().size();
            if (out > in) problems.add("The output has " + out + " alignments and the input " + in);
            // Only duplicates, reads with unmapped mates and discarded pairs can be lost
            if (out < MIN_KEPT_ALIGNMENTS_FRACTION * in) problems.add("The output keeps only " + out + " of " + in + " alignments");
        }
        checkRequirement("outputs_keep_most_reads", problems);
    }

    @Test
    void runIsReproducible() throws Exception {
        List<String> problems = new ArrayList<>();
        String outputPrefix = new File(getModeOutputDir(), "repeat").getPath();
        AnonymizerRunner.RunResult second = AnonymizerRunner.run(dataset.getNormalFile(), dataset.getTumorFile(),
                dataset.getReferenceFile(), outputPrefix, getExtraArgs());
        assertEquals(0, second.exitCode(), "Second run failed. Log: " + second.logFile());
        File secondNormal = second.normalFixVAFOutput().exists() ? second.normalFixVAFOutput() : second.normalOutput();
        File secondTumor = second.tumorFixVAFOutput().exists() ? second.tumorFixVAFOutput() : second.tumorOutput();
        compareRecords(normalOut, AlignmentEvidence.load(secondNormal, dataset.getReference(), false), "normal", problems);
        compareRecords(tumorOut, AlignmentEvidence.load(secondTumor, dataset.getReference(), false), "tumor", problems);
        checkRequirement("run_reproducible_with_same_seed", problems);
    }

    private void compareRecords(AlignmentEvidence first, AlignmentEvidence second, String sample, List<String> problems) {
        if (first.getRecords().size() != second.getRecords().size()) {
            problems.add("Two runs wrote " + first.getRecords().size() + " and " + second.getRecords().size() + " alignments for the " + sample);
            return;
        }
        // Records sharing a start position are not written in a defined order, so the comparison is done on sorted text
        List<String> firstRecords = new ArrayList<>();
        List<String> secondRecords = new ArrayList<>();
        for (SAMRecord record : first.getRecords()) firstRecords.add(record.getSAMString());
        for (SAMRecord record : second.getRecords()) secondRecords.add(record.getSAMString());
        Collections.sort(firstRecords);
        Collections.sort(secondRecords);
        for (int i = 0; i < firstRecords.size(); i++) {
            if (!firstRecords.get(i).equals(secondRecords.get(i))) {
                problems.add("Alignment of the " + sample + " differs between two runs: " + firstRecords.get(i).trim());
            }
        }
    }

    // ---------------------------------------------------------------------------------------------------------
    // Privacy
    // ---------------------------------------------------------------------------------------------------------

    @Test
    void noNewGermlineLeaksInGenomeScan() {
        assertNoNewFindings(OutputValidator.CATEGORY_LEAK, "germline alleles leaked to the output");
    }

    @Test
    void noNewArtifactsInGenomeScan() {
        assertNoNewFindings(OutputValidator.CATEGORY_ARTIFACT, "alleles in the output that the input does not have");
    }

    private void assertNoNewFindings(String category, String description) {
        Set<String> current = OutputValidator.keysWith(genomeScan, "category", category);
        Set<String> known = OutputValidator.keysWith(baselineGenomeScan, "category", category);
        Set<String> newFindings = new TreeSet<>(current);
        newFindings.removeAll(known);
        StringBuilder message = new StringBuilder();
        for (String key : newFindings) {
            message.append(String.join("\t", genomeScan.getRow(key).values())).append("\n");
        }
        assertTrue(newFindings.isEmpty(), newFindings.size() + " new " + description + " (not in the baseline). See "
                + new File(getModeReportDir(), GENOME_SCAN_FILE) + "\n" + String.join("\t", genomeScan.getHeader()) + "\n" + message);
    }

    @Test
    void noNewLeaksOfSpikedGermlineVariants() {
        Set<String> current = OutputValidator.keysWith(smallVariants, "verdict", OutputValidator.VERDICT_LEAK);
        Set<String> known = OutputValidator.keysWith(baselineSmallVariants, "verdict", OutputValidator.VERDICT_LEAK);
        Set<String> newLeaks = new TreeSet<>(current);
        newLeaks.removeAll(known);
        assertTrue(newLeaks.isEmpty(), "Spiked-in variants that must be removed and are still in the output (not in the baseline): "
                + newLeaks + ". See " + new File(getModeReportDir(), SMALL_VARIANTS_FILE));
        Set<String> currentSV = OutputValidator.keysWith(structuralVariants, "verdict", OutputValidator.VERDICT_RESIDUAL_SIGNAL);
        Set<String> knownSV = OutputValidator.keysWith(baselineStructuralVariants, "verdict", OutputValidator.VERDICT_RESIDUAL_SIGNAL);
        Set<String> newResidualSignals = new TreeSet<>(currentSV);
        newResidualSignals.removeAll(knownSV);
        assertTrue(newResidualSignals.isEmpty(), "Spiked-in germline structural variants with signal left in the output (not in the baseline): "
                + newResidualSignals + ". See " + new File(getModeReportDir(), STRUCTURAL_VARIANTS_FILE));
    }

    // ---------------------------------------------------------------------------------------------------------
    // Recorded behaviour
    // ---------------------------------------------------------------------------------------------------------

    @Test
    void resultsMatchBaseline() {
        if (baselineDifferences.size() == 0) return;
        String message = baselineDifferences.size() + " differences against the baseline of mode " + getMode() + ", listed in "
                + new File(getModeReportDir(), BASELINE_DIFFERENCES_FILE)
                + ". If they are the intended effect of the change, accept them with ./gradlew simulatedTest -PupdateBaseline";
        if (SimulatedTestConfig.strictBaseline()) fail(message);
        System.out.println("[" + getMode() + "] WARNING: " + message);
    }

    // ---------------------------------------------------------------------------------------------------------
    // Reports
    // ---------------------------------------------------------------------------------------------------------

    private void writeReports() throws IOException {
        File reportDir = getModeReportDir();
        genomeScan.write(new File(reportDir, GENOME_SCAN_FILE));
        structuralVariants.write(new File(reportDir, STRUCTURAL_VARIANTS_FILE));
        baselineDifferences.write(new File(reportDir, BASELINE_DIFFERENCES_FILE));
        // The report of small variants also shows the VAFs of the baseline next to the ones of this run
        List<String> header = new ArrayList<>(smallVariants.getHeader());
        header.add("baseline_normal_out_vaf");
        header.add("baseline_tumor_out_vaf");
        header.add("baseline_verdict");
        Table extended = new Table(header.toArray(new String[0]));
        for (Map.Entry<String, Map<String, String>> entry : smallVariants.getRows().entrySet()) {
            List<Object> values = new ArrayList<>(entry.getValue().values());
            Map<String, String> baselineRow = baselineSmallVariants.getRow(entry.getKey());
            values.add(baselineRow == null ? OutputValidator.NOT_AVAILABLE : baselineRow.get("normal_out_vaf"));
            values.add(baselineRow == null ? OutputValidator.NOT_AVAILABLE : baselineRow.get("tumor_out_vaf"));
            values.add(baselineRow == null ? OutputValidator.NOT_AVAILABLE : baselineRow.get("verdict"));
            extended.addRow(values.toArray());
        }
        extended.write(new File(reportDir, SMALL_VARIANTS_FILE));
    }

    private double[] meanAbsoluteVafChange(String firstColumn, Table firstTable, String secondColumn, Table secondTable) {
        double sum = 0;
        int n = 0;
        for (Map.Entry<String, Map<String, String>> entry : firstTable.getRows().entrySet()) {
            Map<String, String> secondRow = secondTable.getRow(entry.getKey());
            if (secondRow == null) continue;
            String first = entry.getValue().get(firstColumn);
            String second = secondRow.get(secondColumn);
            if (OutputValidator.NOT_AVAILABLE.equals(first) || OutputValidator.NOT_AVAILABLE.equals(second)) continue;
            if (!OutputValidator.EXPECT_KEEP.equals(entry.getValue().get("expectation"))) continue;
            sum += Math.abs(Double.parseDouble(first) - Double.parseDouble(second));
            n++;
        }
        return new double[]{n == 0 ? 0 : sum / n, n};
    }

    private Map<String, Integer> countBy(Table table, String... columns) {
        Map<String, Integer> answer = new TreeMap<>();
        for (Map<String, String> row : table.getRows().values()) {
            StringBuilder key = new StringBuilder();
            for (String column : columns) {
                if (key.length() > 0) key.append(" / ");
                key.append(row.get(column));
            }
            answer.merge(key.toString(), 1, Integer::sum);
        }
        return answer;
    }

    @AfterAll
    void writeSummary() throws IOException {
        if (smallVariants == null || baselineSmallVariants == null) return;
        if (SimulatedTestConfig.updateBaseline()) {
            Table failures = new Table("requirement", "problems");
            for (Map.Entry<String, Integer> entry : requirementProblems.entrySet()) {
                if (entry.getValue() > 0) failures.addRow(entry.getKey(), entry.getValue());
            }
            failures.write(new File(getModeBaselineDir(), KNOWN_FAILURES_FILE));
            knownFailures = readKnownFailures();
        }
        StringBuilder summary = new StringBuilder();
        summary.append("Simulated validation, mode ").append(getMode()).append("\n");
        summary.append("Arguments: ").append(String.join(" ", getExtraArgs())).append("\n");
        summary.append(String.format(Locale.ROOT, "Run time: %.1f seconds with %d threads%n", run.seconds(), SimulatedTestConfig.getThreads()));
        summary.append("Evaluated outputs: ").append(normalOutputFile.getName()).append(", ").append(tumorOutputFile.getName()).append("\n");
        summary.append("Alignments normal: ").append(normalInAll.getRecords().size()).append(" in, ").append(normalOut.getRecords().size()).append(" out\n");
        summary.append("Alignments tumor: ").append(tumorInAll.getRecords().size()).append(" in, ").append(tumorOut.getRecords().size()).append(" out\n");
        summary.append("\nRequirements with problems (number of problems, number in the baseline):\n");
        boolean anyFailure = false;
        for (Map.Entry<String, Integer> entry : requirementProblems.entrySet()) {
            int known = knownFailures.getOrDefault(entry.getKey(), 0);
            if (entry.getValue() == 0 && known == 0) continue;
            anyFailure = true;
            summary.append("  ").append(entry.getKey()).append("\t").append(entry.getValue()).append("\t").append(known);
            if (entry.getValue() == 0) summary.append("\tFIXED, accept it with -PupdateBaseline");
            summary.append("\n");
        }
        if (!anyFailure) summary.append("  none\n");
        summary.append("\nGenome scan findings (category / reason / kind):\n");
        appendCounts(summary, countBy(genomeScan, "category", "reason", "kind"));
        summary.append("\nSpiked-in small variants (origin / expectation / verdict):\n");
        appendCounts(summary, countBy(smallVariants, "origin", "expectation", "verdict"));
        summary.append("\nSpiked-in structural variants (id / expectation / verdict):\n");
        appendCounts(summary, countBy(structuralVariants, "id", "expectation", "verdict"));
        double[] inputChange = meanAbsoluteVafChange("tumor_in_vaf", smallVariants, "tumor_out_vaf", smallVariants);
        double[] baselineChange = meanAbsoluteVafChange("tumor_out_vaf", smallVariants, "tumor_out_vaf", baselineSmallVariants);
        summary.append(String.format(Locale.ROOT, "%nTumor VAF of variants expected to be kept: mean absolute change of %.4f from the input"
                + " and of %.4f from the baseline (%d variants)%n", inputChange[0], baselineChange[0], (int) inputChange[1]));
        summary.append("Differences against the baseline: ").append(baselineDifferences.size()).append("\n");
        File summaryFile = new File(getModeReportDir(), SUMMARY_FILE);
        try (PrintWriter writer = new PrintWriter(summaryFile)) {
            writer.print(summary);
        }
        System.out.println("\n" + summary + "Reports: " + getModeReportDir() + "\n");
    }

    private void appendCounts(StringBuilder summary, Map<String, Integer> counts) {
        if (counts.isEmpty()) summary.append("  none\n");
        for (Map.Entry<String, Integer> entry : counts.entrySet()) {
            summary.append("  ").append(entry.getValue()).append("\t").append(entry.getKey()).append("\n");
        }
    }
}
