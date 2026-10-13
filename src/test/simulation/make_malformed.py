#!/usr/bin/env python3
"""
Builds the malformed datasets of the simulated validation tests from the well-formed simulated tumor-normal pair.

Every case is a directory with normal.bam, tumor.bam and reference.fa where exactly one thing is broken, so the
failure of GenomeAnonymizer can be attributed to it. To keep the cases small, only the reads of one contig are used.
The expected behaviour of each case is not decided here, it is declared in cases.tsv next to the datasets.
"""
import argparse
import os
import random
import shutil

import pysam

SUBSET_CONTIG = "sim4"
SEED = 7


def subset_records(path):
    with pysam.AlignmentFile(path, "rb") as bam:
        header = bam.header.to_dict()
        records = [record for record in bam.fetch(SUBSET_CONTIG) if record.next_reference_name == SUBSET_CONTIG]
    return header, records


def write_bam(path, header, records, index=True, mode="wb"):
    with pysam.AlignmentFile(path, mode, header=header) as out:
        for record in records:
            out.write(record)
    if index:
        pysam.index(path)


def copy_reference(source_dir, case_dir, with_index=True):
    shutil.copy(os.path.join(source_dir, "reference.fa"), os.path.join(case_dir, "reference.fa"))
    if with_index:
        shutil.copy(os.path.join(source_dir, "reference.fa.fai"), os.path.join(case_dir, "reference.fa.fai"))


def new_case(out_dir, name, source_dir, tumor, reference_index=True):
    """Creates the directory of a case with a correct tumor BAM and reference, and returns its path."""
    case_dir = os.path.join(out_dir, name)
    os.makedirs(case_dir)
    copy_reference(source_dir, case_dir, reference_index)
    write_bam(os.path.join(case_dir, "tumor.bam"), tumor[0], tumor[1])
    return case_dir


def primary_proper_records(records):
    return [i for i, record in enumerate(records)
            if not record.is_supplementary and not record.is_secondary and not record.is_duplicate
            and record.cigarstring == str(record.query_length) + "M"]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--indir", required=True, help="Directory of the well-formed dataset")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory of the malformed datasets")
    args = parser.parse_args()
    rng = random.Random(SEED)
    os.makedirs(args.outdir, exist_ok=True)
    normal = subset_records(os.path.join(args.indir, "normal.bam"))
    tumor = subset_records(os.path.join(args.indir, "tumor.bam"))
    header, records = normal

    # Control: the same subset without any defect, to prove that the failures come from the defects
    case_dir = new_case(args.outdir, "control_wellformed", args.indir, tumor)
    write_bam(os.path.join(case_dir, "normal.bam"), header, records)

    # Normal BAM is not sorted by coordinate (and therefore cannot be indexed)
    case_dir = new_case(args.outdir, "unsorted_normal", args.indir, tumor)
    shuffled = list(records)
    rng.shuffle(shuffled)
    unsorted_header = dict(header)
    unsorted_header["HD"] = dict(header["HD"], SO="unsorted")
    write_bam(os.path.join(case_dir, "normal.bam"), unsorted_header, shuffled, index=False)

    # Normal BAM is sorted but its index is missing
    case_dir = new_case(args.outdir, "missing_index_normal", args.indir, tumor)
    write_bam(os.path.join(case_dir, "normal.bam"), header, records, index=False)

    # Normal BAM is truncated in the middle of the file, the index belongs to the complete file
    case_dir = new_case(args.outdir, "truncated_normal", args.indir, tumor)
    path = os.path.join(case_dir, "normal.bam")
    write_bam(path, header, records)
    with open(path, "rb") as stream:
        content = stream.read()
    with open(path, "wb") as stream:
        stream.write(content[:int(len(content) * 0.6)])

    # Normal file is not a BAM at all
    case_dir = new_case(args.outdir, "not_a_bam_normal", args.indir, tumor)
    with open(os.path.join(case_dir, "normal.bam"), "w") as stream:
        stream.write("This is not an alignment file\n" * 20)

    # Normal file is empty
    case_dir = new_case(args.outdir, "empty_file_normal", args.indir, tumor)
    open(os.path.join(case_dir, "normal.bam"), "w").close()

    # Normal file does not exist
    new_case(args.outdir, "nonexistent_normal", args.indir, tumor)

    # Reads with IUPAC ambiguity codes other than N
    case_dir = new_case(args.outdir, "invalid_bases_normal", args.indir, tumor)
    modified = [pysam.AlignedSegment.from_dict(record.to_dict(), record.header) for record in records]
    for i in rng.sample(primary_proper_records(modified), 5):
        record = modified[i]
        quals = record.query_qualities
        bases = list(record.query_sequence)
        bases[40] = "R"
        bases[90] = "Y"
        record.query_sequence = "".join(bases)
        record.query_qualities = quals
    write_bam(os.path.join(case_dir, "normal.bam"), header, modified)

    # Reads whose CIGAR is longer than their sequence. samtools refuses to index such a file, so both this BAM and
    # a correct twin are written without compression, which makes the index of the twin valid for both
    case_dir = new_case(args.outdir, "cigar_longer_than_sequence_normal", args.indir, tumor)
    path = os.path.join(case_dir, "normal.bam")
    write_bam(path, header, records, mode="wb0")
    modified = [pysam.AlignedSegment.from_dict(record.to_dict(), record.header) for record in records]
    for i in rng.sample(primary_proper_records(modified), 5):
        record = modified[i]
        record.cigartuples = [(0, record.query_length + 30)]
    write_bam(path, header, modified, index=False, mode="wb0")

    # Reads without base qualities
    case_dir = new_case(args.outdir, "missing_qualities_normal", args.indir, tumor)
    modified = [pysam.AlignedSegment.from_dict(record.to_dict(), record.header) for record in records]
    for record in modified:
        record.query_qualities = None
    write_bam(os.path.join(case_dir, "normal.bam"), header, modified)

    # Reads aligned beyond the end of their contig
    case_dir = new_case(args.outdir, "read_beyond_contig_end_normal", args.indir, tumor)
    modified = [pysam.AlignedSegment.from_dict(record.to_dict(), record.header) for record in records]
    contig_length = [sq["LN"] for sq in header["SQ"] if sq["SN"] == SUBSET_CONTIG][0]
    last = primary_proper_records(modified)[-1]
    modified[last].reference_start = contig_length - 50
    modified.sort(key=lambda record: record.reference_start)
    write_bam(os.path.join(case_dir, "normal.bam"), header, modified)

    # Normal BAM aligned to a reference with other contig names
    case_dir = new_case(args.outdir, "contig_names_mismatch_normal", args.indir, tumor)
    renamed_header = dict(header)
    renamed_header["SQ"] = [dict(sq, SN="chr" + sq["SN"]) for sq in header["SQ"]]
    write_bam(os.path.join(case_dir, "normal.bam"), renamed_header, records)

    # Normal BAM aligned to a reference with other contig lengths (same names)
    case_dir = new_case(args.outdir, "contig_lengths_mismatch_normal", args.indir, tumor)
    resized_header = dict(header)
    resized_header["SQ"] = [dict(sq, LN=sq["LN"] + 1000) for sq in header["SQ"]]
    write_bam(os.path.join(case_dir, "normal.bam"), resized_header, records)

    # Single-end reads: pairing information removed
    case_dir = new_case(args.outdir, "single_end_normal", args.indir, tumor)
    modified = [pysam.AlignedSegment.from_dict(record.to_dict(), record.header) for record in records]
    for record in modified:
        record.flag = record.flag & ~(0x1 | 0x2 | 0x8 | 0x20 | 0x40 | 0x80)
        record.next_reference_id = -1
        record.next_reference_start = -1
        record.template_length = 0
    write_bam(os.path.join(case_dir, "normal.bam"), header, modified)

    # Normal BAM without read groups in the header nor in the reads
    case_dir = new_case(args.outdir, "no_read_groups_normal", args.indir, tumor)
    modified = [pysam.AlignedSegment.from_dict(record.to_dict(), record.header) for record in records]
    for record in modified:
        record.set_tag("RG", None)
    headerless = {key: value for key, value in header.items() if key != "RG"}
    write_bam(os.path.join(case_dir, "normal.bam"), headerless, modified)

    # Tumor BAM without index, normal correct
    case_dir = new_case(args.outdir, "missing_index_tumor", args.indir, tumor)
    write_bam(os.path.join(case_dir, "normal.bam"), header, records)
    os.remove(os.path.join(case_dir, "tumor.bam.bai"))

    # Reference without its .fai index
    case_dir = new_case(args.outdir, "missing_reference_index", args.indir, tumor, reference_index=False)
    write_bam(os.path.join(case_dir, "normal.bam"), header, records)

    # Reference of another genome: contigs of the alignments are missing
    case_dir = new_case(args.outdir, "reference_of_other_genome", args.indir, tumor, reference_index=False)
    write_bam(os.path.join(case_dir, "normal.bam"), header, records)
    other = pysam.FastaFile(os.path.join(args.indir, "reference.fa"))
    with open(os.path.join(case_dir, "reference.fa"), "w") as out:
        for name in ("sim1", "sim2"):
            sequence = other.fetch(name)
            out.write(">other_" + name + "\n")
            for i in range(0, len(sequence), 60):
                out.write(sequence[i:i + 60] + "\n")
    other.close()
    pysam.faidx(os.path.join(case_dir, "reference.fa"))

    print("Malformed cases: " + ", ".join(sorted(os.listdir(args.outdir))))


if __name__ == "__main__":
    main()
