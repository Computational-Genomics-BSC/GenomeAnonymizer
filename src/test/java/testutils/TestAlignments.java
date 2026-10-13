package testutils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Helper functions to build small references, read alignments and indexed BAM files for the unit tests
 * @author Nicolas Gaitan
 */
public class TestAlignments {

    public static final byte DEFAULT_BASE_QUALITY = 30;
    public static final int DEFAULT_MAPPING_QUALITY = 60;
    public static final String DEFAULT_READ_GROUP = "test";
    private static final String BASES = "ACGT";

    public static String randomSequence(int length, long seed) {
        Random random = new Random(seed);
        StringBuilder builder = new StringBuilder(length);
        for (int i = 0; i < length; i++) {
            builder.append(BASES.charAt(random.nextInt(BASES.length())));
        }
        return builder.toString();
    }

    /** Returns a base different from the given one, always the same for the same input */
    public static char otherBase(char base) {
        return BASES.charAt((BASES.indexOf(base) + 1) % BASES.length());
    }

    public static SAMFileHeader header(String[] contigNames, int[] contigLengths) {
        SAMFileHeader header = new SAMFileHeader();
        for (int i = 0; i < contigNames.length; i++) {
            header.addSequence(new SAMSequenceRecord(contigNames[i], contigLengths[i]));
        }
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        // The readers of the anonymizer take the platform from the first read group
        SAMReadGroupRecord readGroup = new SAMReadGroupRecord(DEFAULT_READ_GROUP);
        readGroup.setSample(DEFAULT_READ_GROUP);
        readGroup.setPlatform("ILLUMINA");
        header.addReadGroup(readGroup);
        return header;
    }

    /**
     * Builds a mapped read of a pair with constant base qualities
     * @param start 1-based alignment start
     * @param mateStart 1-based alignment start of the mate
     */
    public static SAMRecord read(SAMFileHeader header, String name, String contig, int start, String cigar, String bases,
                                 boolean firstOfPair, boolean negativeStrand, String mateContig, int mateStart,
                                 boolean mateNegativeStrand, int insertSize) {
        SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName(contig);
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setReadString(bases);
        byte[] qualities = new byte[bases.length()];
        Arrays.fill(qualities, DEFAULT_BASE_QUALITY);
        record.setBaseQualities(qualities);
        record.setMappingQuality(DEFAULT_MAPPING_QUALITY);
        record.setReadPairedFlag(true);
        record.setProperPairFlag(true);
        record.setFirstOfPairFlag(firstOfPair);
        record.setSecondOfPairFlag(!firstOfPair);
        record.setReadNegativeStrandFlag(negativeStrand);
        record.setMateReferenceName(mateContig);
        record.setMateAlignmentStart(mateStart);
        record.setMateNegativeStrandFlag(mateNegativeStrand);
        record.setInferredInsertSize(insertSize);
        record.setAttribute("RG", DEFAULT_READ_GROUP);
        return record;
    }

    /**
     * Builds the two reads of a proper pair without differences against the reference
     * @param reference sequence of the contig
     * @param start 1-based alignment start of the first read
     * @param mateStart 1-based alignment start of the second read
     */
    public static SAMRecord[] properPair(SAMFileHeader header, String name, String contig, String reference, int start,
                                         int mateStart, int readLength) {
        int insertSize = mateStart + readLength - start;
        String cigar = readLength + "M";
        SAMRecord first = read(header, name, contig, start, cigar, reference.substring(start - 1, start - 1 + readLength),
                true, false, contig, mateStart, true, insertSize);
        SAMRecord second = read(header, name, contig, mateStart, cigar, reference.substring(mateStart - 1, mateStart - 1 + readLength),
                false, true, contig, start, false, -insertSize);
        return new SAMRecord[]{first, second};
    }

    /** Replaces the base of a read that is aligned to the given 1-based reference position, if the read covers it */
    public static void setBaseAtReferencePosition(SAMRecord record, int referencePosition, char base) {
        int readPosition = record.getReadPositionAtReferencePosition(referencePosition);
        if (readPosition == 0) return;
        byte[] bases = record.getReadBases();
        bases[readPosition - 1] = (byte) base;
        record.setReadBases(bases);
    }

    public static File writeIndexedBam(File directory, String name, SAMFileHeader header, List<SAMRecord> records) {
        File answer = new File(directory, name);
        SAMFileWriterFactory factory = new SAMFileWriterFactory();
        factory.setCreateIndex(true);
        try (SAMFileWriter writer = factory.makeWriter(header, false, answer, null)) {
            for (SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }
        return answer;
    }

    /** Writes a FASTA file with one line per sequence and its .fai index */
    public static File writeReference(File directory, String name, String[] contigNames, String[] sequences) throws IOException {
        File answer = new File(directory, name);
        long offset = 0;
        try (PrintWriter fasta = new PrintWriter(answer);
             PrintWriter index = new PrintWriter(new File(directory, name + ".fai"))) {
            for (int i = 0; i < contigNames.length; i++) {
                String headerLine = ">" + contigNames[i];
                fasta.print(headerLine + "\n" + sequences[i] + "\n");
                offset += headerLine.length() + 1;
                index.print(contigNames[i] + "\t" + sequences[i].length() + "\t" + offset + "\t" + sequences[i].length()
                        + "\t" + (sequences[i].length() + 1) + "\n");
                offset += sequences[i].length() + 1;
            }
        }
        return answer;
    }
}
