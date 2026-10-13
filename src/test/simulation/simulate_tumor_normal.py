#!/usr/bin/env python3
"""
Simulates a small tumor-normal pair of Illumina-like paired-end WGS read sets with spiked-in variation.

The script writes a reference genome, two pairs of FASTQ files (normal and tumor), the truth table of every
spiked-in variant and a BED file describing the special regions of the dataset. Alignment is done afterwards by
build_datasets.sh (bwa mem + samtools), so MAPQ, soft clips, supplementary alignments and discordant pairs are the
ones a real aligner produces from the simulated haplotypes.

Model:
  - Every sample is a mixture of diploid genomes. The normal sample is 100% germline. The tumor sample is made of
    normal cells plus three nested clones (FOUNDER > SUB1 > SUB2), which gives somatic variants different VAFs.
  - A haplotype is a list of segments of the reference (or literal inserted sequence), so SNVs, indels, deletions,
    insertions, inversions, tandem duplications and reciprocal translocations are all applied on reference
    coordinates and reads are sampled from the rearranged sequence.
  - Coverage is modulated per sample with multipliers on reference coordinates, to produce regions covered by only
    one of the samples, by none, or at higher depth.
  - Base qualities follow the binned scheme of current Illumina instruments (37, 23, 12, 2) with quality decaying
    along the read, more in read 2. Substitution errors are drawn from the quality of each base.

Everything is driven by a single seed, so the dataset is reproducible.
"""
import argparse
import bisect
import os
import random

# ---------------------------------------------------------------------------------------------------------------
# Dataset layout. All coordinates in this block are 0-based, half-open, over the reference.
# ---------------------------------------------------------------------------------------------------------------
CONTIGS = [("sim1", 20000), ("sim2", 15000), ("sim3", 10000), ("sim4", 8000)]

READ_LENGTH = 150
INSERT_MEAN = 420
INSERT_SD = 70
INSERT_MIN = 180
BASE_COVERAGE = 10.0
DUPLICATE_FRACTION = 0.04

# Tumor composition: fraction of cells carrying the variants of each clone level
TUMOR_PURITY = 0.8
CLONE_CELL_FRACTION = {"FOUNDER": 1.0, "SUB1": 0.5, "SUB2": 0.2}
LEVELS = ["FOUNDER", "SUB1", "SUB2"]

# Repetitive and ambiguous regions of the reference
SEGDUP_IDENTICAL = (("sim1", 6000, 7200), ("sim3", 3000))
SEGDUP_DIVERGED = (("sim2", 9000, 9800), ("sim4", 5000), 0.02)
STR_CA = ("sim2", 2000, "CA", 40)
HOMOPOLYMER = ("sim1", 15000, "A", 25)
TANDEM_REPEAT = ("sim3", 7000, 35, 12)
N_BLOCK = ("sim4", 7000, 7100)

# Coverage multipliers: (name, contig, start, end, normal multiplier, tumor multiplier)
COVERAGE_REGIONS = [
    ("tumor_only", "sim1", 10000, 11000, 0.0, 1.0),
    ("normal_only", "sim1", 12500, 13300, 1.0, 0.0),
    ("normal_low_tumor_none_1", "sim2", 800, 1500, 0.3, 0.0),
    ("deep", "sim2", 4000, 7000, 1.8, 1.8),
    ("uncovered_both", "sim3", 8000, 8300, 0.0, 0.0),
    ("normal_low_tumor_none_2", "sim4", 1500, 2200, 0.3, 0.0),
    ("tumor_low_normal_none", "sim4", 4000, 4500, 0.0, 0.3),
]

# Structural variants: dict with id, type, origin, level, haps and coordinates
STRUCTURAL_VARIANTS = [
    dict(id="G_SVDEL_1", type="SV_DEL", origin="germline", level="GERMLINE", haps=(0,), contig="sim1", start=3000, end=3400),
    dict(id="G_SVDEL_HOM", type="SV_DEL", origin="germline", level="GERMLINE", haps=(0, 1), contig="sim4", start=6200, end=6500),
    dict(id="S_SVDEL_1", type="SV_DEL", origin="somatic", level="FOUNDER", haps=(1,), contig="sim1", start=17000, end=17600),
    dict(id="G_SVINV_1", type="SV_INV", origin="germline", level="GERMLINE", haps=(0,), contig="sim2", start=11000, end=12200),
    dict(id="S_SVINV_1", type="SV_INV", origin="somatic", level="FOUNDER", haps=(0,), contig="sim3", start=5000, end=6000),
    dict(id="S_SVDUP_1", type="SV_DUP", origin="somatic", level="FOUNDER", haps=(0,), contig="sim2", start=13000, end=13500),
    dict(id="G_SVINS_1", type="SV_INS", origin="germline", level="GERMLINE", haps=(1,), contig="sim3", start=1500, end=1500, length=130),
    dict(id="S_SVINS_1", type="SV_INS", origin="somatic", level="FOUNDER", haps=(1,), contig="sim4", start=3000, end=3000, length=120),
    dict(id="G_SVTRA_1", type="SV_TRA", origin="germline", level="GERMLINE", haps=(1,), contig="sim2", start=14200, end=14200,
         mate_contig="sim4", mate_pos=600),
    dict(id="S_SVTRA_1", type="SV_TRA", origin="somatic", level="FOUNDER", haps=(0,), contig="sim1", start=18800, end=18800,
         mate_contig="sim3", mate_pos=9300),
]

# Small variants placed by hand in the special regions: (id, type, origin, level, haps, contig, start, extra)
# extra: for SNV None, for DEL the deleted length, for INS the inserted length (or the literal sequence)
EXPLICIT_SMALL_VARIANTS = [
    # Region covered only by the tumor: nothing can be told apart from germline there
    ("G_SNV_TONLY_HET", "SNV", "germline", "GERMLINE", (0,), "sim1", 10200, None),
    ("G_SNV_TONLY_HOM", "SNV", "germline", "GERMLINE", (0, 1), "sim1", 10400, None),
    ("G_DEL_TONLY_HOM", "DEL", "germline", "GERMLINE", (0, 1), "sim1", 10600, 5),
    ("S_SNV_TONLY_HET", "SNV", "somatic", "FOUNDER", (1,), "sim1", 10300, None),
    ("S_SNV_TONLY_HOM", "SNV", "somatic", "FOUNDER", (0, 1), "sim1", 10800, None),
    # Region covered only by the normal
    ("G_SNV_NONLY_HOM", "SNV", "germline", "GERMLINE", (0, 1), "sim1", 12700, None),
    ("G_SNV_NONLY_HET", "SNV", "germline", "GERMLINE", (1,), "sim1", 12900, None),
    ("G_INS_NONLY_HOM", "INS", "germline", "GERMLINE", (0, 1), "sim1", 13100, 3),
    # Regions with a few normal reads and no tumor reads
    ("G_SNV_NLOW1_HOM_A", "SNV", "germline", "GERMLINE", (0, 1), "sim2", 1000, None),
    ("G_SNV_NLOW1_HOM_B", "SNV", "germline", "GERMLINE", (0, 1), "sim2", 1200, None),
    ("G_SNV_NLOW1_HET", "SNV", "germline", "GERMLINE", (0,), "sim2", 1300, None),
    ("G_SNV_NLOW2_HOM", "SNV", "germline", "GERMLINE", (0, 1), "sim4", 1800, None),
    ("G_DEL_NLOW2_HOM", "DEL", "germline", "GERMLINE", (0, 1), "sim4", 2000, 2),
    # Region with a few tumor reads and no normal reads
    ("G_SNV_TLOW_HOM", "SNV", "germline", "GERMLINE", (0, 1), "sim4", 4200, None),
    ("S_SNV_TLOW_HOM", "SNV", "somatic", "FOUNDER", (0, 1), "sim4", 4300, None),
    # Repeats
    ("G_SNV_SEGDUP_HOM", "SNV", "germline", "GERMLINE", (0, 1), "sim1", 6500, None),
    ("S_SNV_SEGDUP", "SNV", "somatic", "FOUNDER", (0,), "sim1", 6800, None),
    ("G_SNV_DIVDUP_HET", "SNV", "germline", "GERMLINE", (1,), "sim2", 9400, None),
    ("S_SNV_DIVDUP", "SNV", "somatic", "FOUNDER", (1,), "sim4", 5400, None),
    ("G_DEL_STR_HET", "DEL", "germline", "GERMLINE", (0,), "sim2", 2000, 4),
    ("S_INS_HOMOPOLYMER", "INS", "somatic", "FOUNDER", (1,), "sim1", 15000, "AA"),
    ("G_DEL_TANDEM_HOM", "DEL", "germline", "GERMLINE", (0, 1), "sim3", 7035, 35),
    # Around the block of N in the reference
    ("G_SNV_NBLOCK_HOM", "SNV", "germline", "GERMLINE", (0, 1), "sim4", 6985, None),
    ("S_SNV_NBLOCK", "SNV", "somatic", "FOUNDER", (0,), "sim4", 7110, None),
    # Clustered variants
    ("G_SNV_MNV_A", "SNV", "germline", "GERMLINE", (0,), "sim1", 1000, None),
    ("G_SNV_MNV_B", "SNV", "germline", "GERMLINE", (0,), "sim1", 1001, None),
    ("G_SNV_NEAR_SOM_CIS", "SNV", "germline", "GERMLINE", (0,), "sim1", 1500, None),
    ("S_SNV_NEAR_GERM_CIS", "SNV", "somatic", "FOUNDER", (0,), "sim1", 1503, None),
    ("G_SNV_NEAR_SOM_TRANS", "SNV", "germline", "GERMLINE", (0,), "sim1", 2000, None),
    ("S_SNV_NEAR_GERM_TRANS", "SNV", "somatic", "FOUNDER", (1,), "sim1", 2004, None),
    ("G_DEL_NEAR_SOM_HOM", "DEL", "germline", "GERMLINE", (0, 1), "sim3", 500, 3),
    ("S_SNV_NEAR_GERM_DEL", "SNV", "somatic", "FOUNDER", (0,), "sim3", 506, None),
    # Contig ends
    ("G_SNV_START_HOM", "SNV", "germline", "GERMLINE", (0, 1), "sim3", 40, None),
    ("G_SNV_END_HET", "SNV", "germline", "GERMLINE", (1,), "sim1", 19960, None),
    ("S_SNV_END_HOM", "SNV", "somatic", "FOUNDER", (0, 1), "sim4", 7955, None),
    # Inside structural variants
    ("G_SNV_IN_S_INV_HOM", "SNV", "germline", "GERMLINE", (0, 1), "sim3", 5500, None),
    ("S_SNV_IN_G_INV", "SNV", "somatic", "FOUNDER", (1,), "sim2", 11600, None),
    ("G_SNV_IN_S_DUP_HET", "SNV", "germline", "GERMLINE", (1,), "sim2", 13250, None),
]

# Randomly placed small variants: (id prefix, type, origin, level, zygosity, count, only inside the deep region)
RANDOM_SMALL_VARIANTS = [
    ("G_SNV_HET", "SNV", "germline", "GERMLINE", "het", 30, False),
    ("G_SNV_HOM", "SNV", "germline", "GERMLINE", "hom", 15, False),
    ("G_INDEL_HET", "INDEL", "germline", "GERMLINE", "het", 8, False),
    ("G_INDEL_HOM", "INDEL", "germline", "GERMLINE", "hom", 6, False),
    ("S_SNV_FOUNDER", "SNV", "somatic", "FOUNDER", "het", 14, False),
    ("S_SNV_SUB1", "SNV", "somatic", "SUB1", "het", 14, False),
    ("S_SNV_SUB2", "SNV", "somatic", "SUB2", "het", 10, False),
    ("S_SNV_FOUNDER_HOM", "SNV", "somatic", "FOUNDER", "hom", 3, False),
    ("S_INDEL_FOUNDER", "INDEL", "somatic", "FOUNDER", "het", 4, False),
    ("S_INDEL_SUB1", "INDEL", "somatic", "SUB1", "het", 3, False),
    ("S_INDEL_SUB2", "INDEL", "somatic", "SUB2", "het", 2, False),
    ("G_SNV_DEEP_HET", "SNV", "germline", "GERMLINE", "het", 5, True),
    ("S_SNV_DEEP_FOUNDER", "SNV", "somatic", "FOUNDER", "het", 5, True),
    ("S_SNV_DEEP_SUB1", "SNV", "somatic", "SUB1", "het", 6, True),
    ("S_SNV_DEEP_SUB2", "SNV", "somatic", "SUB2", "het", 6, True),
    ("S_SNV_DEEP_FOUNDER_HOM", "SNV", "somatic", "FOUNDER", "hom", 2, True),
    ("S_INDEL_DEEP_SUB1", "INDEL", "somatic", "SUB1", "het", 2, True),
]
INDEL_LENGTHS = [1, 1, 1, 2, 2, 3, 4, 5, 8, 12]
RANDOM_VARIANT_MIN_SPACING = 30
RANDOM_VARIANT_SV_MARGIN = 250
RANDOM_VARIANT_CONTIG_MARGIN = 200

# Sequencing qualities (binned Illumina scheme)
QUALITY_BINS = [37, 23, 12, 2]
LOW_QUALITY_SPLIT = [0.62, 0.37, 0.01]

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def random_sequence(rng, length, gc=0.41):
    weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choices(BASES, weights=weights, k=length))


def other_base(rng, base):
    return rng.choice([b for b in BASES if b != base.upper()])


# ---------------------------------------------------------------------------------------------------------------
# Reference
# ---------------------------------------------------------------------------------------------------------------
def build_reference(rng):
    """
    Returns two dicts of contig -> sequence: the donor reference, where every base is known, and the published
    reference, that is the same but with the N block masked.
    """
    donor = {name: list(random_sequence(rng, length)) for name, length in CONTIGS}
    # Identical segmental duplication
    (src, s, e), (dst, d) = SEGDUP_IDENTICAL
    donor[dst][d:d + (e - s)] = donor[src][s:e]
    # Diverged segmental duplication
    (src, s, e), (dst, d), divergence = SEGDUP_DIVERGED
    copy = donor[src][s:e]
    for i in range(len(copy)):
        if rng.random() < divergence:
            copy[i] = other_base(rng, copy[i])
    donor[dst][d:d + (e - s)] = copy
    # Short tandem repeat, homopolymer and minisatellite-like tandem repeat
    contig, pos, unit, copies = STR_CA
    donor[contig][pos:pos + len(unit) * copies] = list(unit * copies)
    contig, pos, unit, copies = HOMOPOLYMER
    donor[contig][pos:pos + copies] = list(unit * copies)
    contig, pos, unit_length, copies = TANDEM_REPEAT
    unit = random_sequence(rng, unit_length)
    donor[contig][pos:pos + unit_length * copies] = list(unit * copies)
    donor = {name: "".join(seq) for name, seq in donor.items()}
    published = dict(donor)
    contig, s, e = N_BLOCK
    published[contig] = donor[contig][:s] + "N" * (e - s) + donor[contig][e:]
    return donor, published


# ---------------------------------------------------------------------------------------------------------------
# Variants
# ---------------------------------------------------------------------------------------------------------------
def expected_vafs(origin, level, haps):
    allele_fraction = len(haps) / 2.0
    if origin == "germline":
        return allele_fraction, allele_fraction
    return 0.0, allele_fraction * TUMOR_PURITY * CLONE_CELL_FRACTION[level]


def make_small_variant(rng, donor, var_id, var_type, origin, level, haps, contig, start, extra):
    """Builds the dict describing a SNV, DEL or INS. start is the first affected base (0-based)."""
    seq = donor[contig]
    variant = dict(id=var_id, type=var_type, origin=origin, level=level, haps=tuple(haps), contig=contig, start=start)
    if var_type == "SNV":
        variant["end"] = start + 1
        variant["ref"] = seq[start]
        variant["alt"] = other_base(rng, seq[start])
        variant["vcf_pos"] = start + 1
    elif var_type == "DEL":
        length = extra
        variant["end"] = start + length
        variant["ref"] = seq[start - 1:start + length]
        variant["alt"] = seq[start - 1]
        variant["vcf_pos"] = start
    else:
        inserted = extra if isinstance(extra, str) else random_sequence(rng, extra)
        variant["end"] = start
        variant["inserted"] = inserted
        variant["ref"] = seq[start - 1]
        variant["alt"] = seq[start - 1] + inserted
        variant["vcf_pos"] = start
    return variant


def forbidden_intervals():
    """Intervals of the reference where no random variant is placed."""
    forbidden = {name: [] for name, _ in CONTIGS}
    for name, length in CONTIGS:
        forbidden[name].append((0, RANDOM_VARIANT_CONTIG_MARGIN))
        forbidden[name].append((length - RANDOM_VARIANT_CONTIG_MARGIN, length))
    for sv in STRUCTURAL_VARIANTS:
        forbidden[sv["contig"]].append((sv["start"] - RANDOM_VARIANT_SV_MARGIN, sv["end"] + RANDOM_VARIANT_SV_MARGIN))
        if "mate_contig" in sv:
            forbidden[sv["mate_contig"]].append((sv["mate_pos"] - RANDOM_VARIANT_SV_MARGIN,
                                                 sv["mate_pos"] + RANDOM_VARIANT_SV_MARGIN))
    for name, contig, start, end, normal_mult, tumor_mult in COVERAGE_REGIONS:
        if name != "deep":
            forbidden[contig].append((start - 200, end + 200))
    (src, s, e), (dst, d) = SEGDUP_IDENTICAL
    forbidden[src].append((s - 50, e + 50))
    forbidden[dst].append((d - 50, d + (e - s) + 50))
    (src, s, e), (dst, d), divergence = SEGDUP_DIVERGED
    forbidden[src].append((s - 50, e + 50))
    forbidden[dst].append((d - 50, d + (e - s) + 50))
    contig, pos, unit, copies = STR_CA
    forbidden[contig].append((pos - 50, pos + len(unit) * copies + 50))
    contig, pos, unit, copies = HOMOPOLYMER
    forbidden[contig].append((pos - 50, pos + copies + 50))
    contig, pos, unit_length, copies = TANDEM_REPEAT
    forbidden[contig].append((pos - 50, pos + unit_length * copies + 50))
    contig, s, e = N_BLOCK
    forbidden[contig].append((s - 150, e + 150))
    return forbidden


def build_small_variants(rng, donor):
    variants = []
    for var_id, var_type, origin, level, haps, contig, start, extra in EXPLICIT_SMALL_VARIANTS:
        variants.append(make_small_variant(rng, donor, var_id, var_type, origin, level, haps, contig, start, extra))
    forbidden = forbidden_intervals()
    occupied = {name: [] for name, _ in CONTIGS}
    for variant in variants:
        occupied[variant["contig"]].append((variant["start"], variant["end"]))
    deep = [region for region in COVERAGE_REGIONS if region[0] == "deep"][0]
    contig_names = [name for name, _ in CONTIGS]
    contig_lengths = [length for _, length in CONTIGS]

    def is_free(contig, start, end):
        for s, e in forbidden[contig]:
            if start < e and end > s:
                return False
        for s, e in occupied[contig]:
            if start < e + RANDOM_VARIANT_MIN_SPACING and end > s - RANDOM_VARIANT_MIN_SPACING:
                return False
        return True

    for prefix, var_type, origin, level, zygosity, count, in_deep in RANDOM_SMALL_VARIANTS:
        for n in range(1, count + 1):
            while True:
                if in_deep:
                    contig = deep[1]
                    start = rng.randrange(deep[2] + 100, deep[3] - 100)
                else:
                    contig = rng.choices(contig_names, weights=contig_lengths)[0]
                    start = rng.randrange(0, len(donor[contig]))
                    if deep[1] == contig and deep[2] - 200 < start < deep[3] + 200:
                        continue
                concrete_type = var_type
                extra = None
                if var_type == "INDEL":
                    concrete_type = rng.choice(["INS", "DEL"])
                    extra = rng.choice(INDEL_LENGTHS)
                end = start + extra if concrete_type == "DEL" else start + 1
                if is_free(contig, start, end):
                    break
            haps = (0, 1) if zygosity == "hom" else (rng.randrange(2),)
            var_id = prefix.replace("INDEL", concrete_type) + "_" + str(n)
            variant = make_small_variant(rng, donor, var_id, concrete_type, origin, level, haps, contig, start, extra)
            variants.append(variant)
            occupied[contig].append((variant["start"], max(variant["end"], variant["start"] + 1)))
    return variants


def build_structural_variants(rng):
    variants = []
    for sv in STRUCTURAL_VARIANTS:
        variant = dict(sv)
        if sv["type"] == "SV_INS":
            variant["inserted"] = random_sequence(rng, sv["length"])
        variants.append(variant)
    return variants


# ---------------------------------------------------------------------------------------------------------------
# Haplotypes
# ---------------------------------------------------------------------------------------------------------------
class Segment:
    """A piece of a haplotype: either an interval of a reference contig in a given strand, or literal sequence."""
    def __init__(self, contig=None, start=0, end=0, strand=1, literal=None):
        self.contig = contig
        self.start = start
        self.end = end
        self.strand = strand
        self.literal = literal

    def is_literal(self):
        return self.literal is not None

    def length(self):
        return len(self.literal) if self.is_literal() else self.end - self.start

    def reversed(self):
        if self.is_literal():
            return Segment(literal=revcomp(self.literal))
        return Segment(self.contig, self.start, self.end, -self.strand)


class Haplotype:
    def __init__(self, donor):
        self.donor = donor
        self.snvs = {}
        self.chromosomes = {name: [Segment(name, 0, len(seq), 1)] for name, seq in donor.items()}
        self.sequences = None
        self.offsets = None

    def split(self, contig, pos):
        for segments in self.chromosomes.values():
            for i, seg in enumerate(segments):
                if not seg.is_literal() and seg.contig == contig and seg.start < pos < seg.end:
                    if seg.strand != 1:
                        raise ValueError("Cannot split a reversed segment at " + contig + ":" + str(pos))
                    segments[i:i + 1] = [Segment(contig, seg.start, pos, 1), Segment(contig, pos, seg.end, 1)]
                    return

    def find_run(self, contig, start, end):
        """Returns the chromosome key and the indexes [i, j) of the segments covering contig:start-end."""
        self.split(contig, start)
        self.split(contig, end)
        for key, segments in self.chromosomes.items():
            idxs = [i for i, seg in enumerate(segments)
                    if not seg.is_literal() and seg.contig == contig and seg.start >= start and seg.end <= end]
            if idxs:
                return key, idxs[0], idxs[-1] + 1
        raise ValueError("Region not found in haplotype: " + contig + ":" + str(start) + "-" + str(end))

    def find_boundary(self, contig, pos):
        """Returns the chromosome key and the index of the first segment starting at contig:pos."""
        self.split(contig, pos)
        for key, segments in self.chromosomes.items():
            for i, seg in enumerate(segments):
                if not seg.is_literal() and seg.contig == contig and seg.start == pos and seg.strand == 1:
                    return key, i
        raise ValueError("Boundary not found in haplotype: " + contig + ":" + str(pos))

    def apply(self, variant):
        var_type = variant["type"]
        contig = variant["contig"]
        if var_type == "SNV":
            self.snvs[(contig, variant["start"])] = variant["alt"]
        elif var_type in ("DEL", "SV_DEL"):
            key, i, j = self.find_run(contig, variant["start"], variant["end"])
            del self.chromosomes[key][i:j]
        elif var_type in ("INS", "SV_INS"):
            key, i = self.find_boundary(contig, variant["start"])
            self.chromosomes[key].insert(i, Segment(literal=variant["inserted"]))
        elif var_type == "SV_INV":
            key, i, j = self.find_run(contig, variant["start"], variant["end"])
            segments = self.chromosomes[key]
            segments[i:j] = [seg.reversed() for seg in reversed(segments[i:j])]
        elif var_type == "SV_DUP":
            key, i, j = self.find_run(contig, variant["start"], variant["end"])
            segments = self.chromosomes[key]
            segments[j:j] = [Segment(seg.contig, seg.start, seg.end, seg.strand, seg.literal) for seg in segments[i:j]]
        elif var_type == "SV_TRA":
            key1, i = self.find_boundary(contig, variant["start"])
            key2, j = self.find_boundary(variant["mate_contig"], variant["mate_pos"])
            first = self.chromosomes[key1]
            second = self.chromosomes[key2]
            self.chromosomes[key1] = first[:i] + second[j:]
            self.chromosomes[key2] = second[:j] + first[i:]
        else:
            raise ValueError("Unknown variant type " + var_type)

    def materialize(self):
        mutated = {}
        for name, seq in self.donor.items():
            mutated[name] = list(seq)
        for (contig, pos), alt in self.snvs.items():
            mutated[contig][pos] = alt
        mutated = {name: "".join(seq) for name, seq in mutated.items()}
        self.sequences = {}
        self.offsets = {}
        for key, segments in self.chromosomes.items():
            pieces = []
            offsets = []
            total = 0
            for seg in segments:
                offsets.append(total)
                if seg.is_literal():
                    piece = seg.literal
                else:
                    piece = mutated[seg.contig][seg.start:seg.end]
                    if seg.strand == -1:
                        piece = revcomp(piece)
                pieces.append(piece)
                total += len(piece)
            self.sequences[key] = "".join(pieces)
            self.offsets[key] = offsets

    def to_reference(self, key, pos):
        """Maps a position of this haplotype to (contig, position) on the reference, or None for inserted bases."""
        idx = bisect.bisect_right(self.offsets[key], pos) - 1
        seg = self.chromosomes[key][idx]
        if seg.is_literal():
            return None
        offset = pos - self.offsets[key][idx]
        if seg.strand == 1:
            return seg.contig, seg.start + offset
        return seg.contig, seg.end - 1 - offset


def build_genome(donor, small_variants, structural_variants, levels):
    """Builds the two haplotypes of a genome that carries the germline variants plus the given somatic levels."""
    haplotypes = [Haplotype(donor), Haplotype(donor)]
    for variant in small_variants + structural_variants:
        if variant["level"] != "GERMLINE" and variant["level"] not in levels:
            continue
        for hap in variant["haps"]:
            haplotypes[hap].apply(variant)
    for haplotype in haplotypes:
        haplotype.materialize()
    return haplotypes


# ---------------------------------------------------------------------------------------------------------------
# Reads
# ---------------------------------------------------------------------------------------------------------------
def build_coverage_tracks(sample_idx):
    tracks = {name: [1.0] * length for name, length in CONTIGS}
    for name, contig, start, end, normal_mult, tumor_mult in COVERAGE_REGIONS:
        mult = normal_mult if sample_idx == 0 else tumor_mult
        for i in range(start, end):
            tracks[contig][i] = mult
    return tracks


def sequence_read(rng, template, is_second):
    """Adds Illumina-like qualities and substitution errors to the bases of a read."""
    n = len(template)
    bases = []
    quals = []
    errors = 0
    for i in range(n):
        cycle = i / float(n)
        high_quality_prob = 0.93 - 0.15 * cycle * cycle if is_second else 0.955 - 0.10 * cycle * cycle
        if rng.random() < high_quality_prob:
            qual = QUALITY_BINS[0]
        else:
            qual = rng.choices(QUALITY_BINS[1:], weights=LOW_QUALITY_SPLIT)[0]
        base = template[i]
        if qual == 2:
            base = "N"
        elif rng.random() < 10 ** (-qual / 10.0):
            base = other_base(rng, base)
            errors += 1
        bases.append(base)
        quals.append(chr(qual + 33))
    return "".join(bases), "".join(quals), errors


def simulate_sample(rng, genomes, weights, sample_idx, out_prefix, flowcell):
    """Samples read pairs from a mixture of genomes and writes them to a pair of FASTQ files."""
    tracks = build_coverage_tracks(sample_idx)
    max_mult = max(max(track) for track in tracks.values())
    genome_size = sum(length for _, length in CONTIGS)
    candidates = int(BASE_COVERAGE * genome_size * max_mult / (2.0 * READ_LENGTH))
    stats = dict(pairs=0, duplicates=0, bases=0, errors=0)
    serial = 0
    with open(out_prefix + "_R1.fastq", "w") as out1, open(out_prefix + "_R2.fastq", "w") as out2:
        for c in range(candidates):
            genome = rng.choices(genomes, weights=weights)[0]
            haplotype = genome[rng.randrange(2)]
            keys = list(haplotype.sequences.keys())
            key = rng.choices(keys, weights=[len(haplotype.sequences[k]) for k in keys])[0]
            sequence = haplotype.sequences[key]
            fragment_length = max(INSERT_MIN, int(round(rng.gauss(INSERT_MEAN, INSERT_SD))))
            if fragment_length > len(sequence):
                continue
            start = rng.randrange(0, len(sequence) - fragment_length + 1)
            # The coverage multiplier of a fragment is the lowest one among its ends and its midpoint
            mult = max_mult
            for pos in (start, start + fragment_length // 2, start + fragment_length - 1):
                ref_coordinate = haplotype.to_reference(key, pos)
                if ref_coordinate is not None:
                    mult = min(mult, tracks[ref_coordinate[0]][ref_coordinate[1]])
            if rng.random() >= mult / max_mult:
                continue
            fragment = sequence[start:start + fragment_length]
            if rng.random() < 0.5:
                fragment = revcomp(fragment)
            copies = 2 if rng.random() < DUPLICATE_FRACTION else 1
            for copy in range(copies):
                serial += 1
                read1, qual1, errors1 = sequence_read(rng, fragment[:READ_LENGTH], False)
                read2, qual2, errors2 = sequence_read(rng, revcomp(fragment[-READ_LENGTH:]), True)
                tile = 1101 + (serial % 8)
                name = "A00741:" + str(213 + sample_idx) + ":" + flowcell + ":" + str(1 + sample_idx) + ":" + str(tile) \
                       + ":" + str(rng.randrange(1000, 32000)) + ":" + str(rng.randrange(1000, 36000)) + ":" + str(serial)
                out1.write("@" + name + "/1\n" + read1 + "\n+\n" + qual1 + "\n")
                out2.write("@" + name + "/2\n" + read2 + "\n+\n" + qual2 + "\n")
                stats["pairs"] += 1
                stats["bases"] += len(read1) + len(read2)
                stats["errors"] += errors1 + errors2
                if copy > 0:
                    stats["duplicates"] += 1
    return stats


# ---------------------------------------------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------------------------------------------
def write_fasta(path, sequences):
    with open(path, "w") as out:
        for name, _ in CONTIGS:
            out.write(">" + name + "\n")
            seq = sequences[name]
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")


def region_tags(contig, start, end):
    """Names of the special regions overlapping an interval, used to annotate the truth table."""
    tags = []
    for name, region_contig, s, e, normal_mult, tumor_mult in COVERAGE_REGIONS:
        if region_contig == contig and start < e and end > s:
            tags.append(name)
    for name, region_contig, s, e in special_regions():
        if region_contig == contig and start < e and end > s:
            tags.append(name)
    return ",".join(tags) if tags else "."


def special_regions():
    regions = []
    (src, s, e), (dst, d) = SEGDUP_IDENTICAL
    regions.append(("segdup_identical", src, s, e))
    regions.append(("segdup_identical", dst, d, d + (e - s)))
    (src, s, e), (dst, d), divergence = SEGDUP_DIVERGED
    regions.append(("segdup_diverged", src, s, e))
    regions.append(("segdup_diverged", dst, d, d + (e - s)))
    contig, pos, unit, copies = STR_CA
    regions.append(("str_ca", contig, pos, pos + len(unit) * copies))
    contig, pos, unit, copies = HOMOPOLYMER
    regions.append(("homopolymer", contig, pos, pos + copies))
    contig, pos, unit_length, copies = TANDEM_REPEAT
    regions.append(("tandem_repeat", contig, pos, pos + unit_length * copies))
    contig, s, e = N_BLOCK
    regions.append(("n_block", contig, s, e))
    return regions


def write_truth(path, small_variants, structural_variants):
    header = ["id", "type", "origin", "level", "zygosity", "contig", "pos", "end", "ref", "alt", "length",
              "mate_contig", "mate_pos", "expected_normal_vaf", "expected_tumor_vaf", "regions"]
    with open(path, "w") as out:
        out.write("#" + "\t".join(header) + "\n")
        for variant in small_variants + structural_variants:
            normal_vaf, tumor_vaf = expected_vafs(variant["origin"], variant["level"], variant["haps"])
            zygosity = "hom" if len(variant["haps"]) == 2 else "het"
            var_type = variant["type"]
            # pos and end are 1-based inclusive. For indels pos is the anchor base before the event, as in VCF
            if var_type == "SNV":
                pos, end, ref, alt, length = variant["vcf_pos"], variant["vcf_pos"], variant["ref"], variant["alt"], 1
            elif var_type == "DEL":
                length = variant["end"] - variant["start"]
                pos, end, ref, alt = variant["vcf_pos"], variant["vcf_pos"] + length, variant["ref"], variant["alt"]
            elif var_type == "INS":
                length = len(variant["inserted"])
                pos, end, ref, alt = variant["vcf_pos"], variant["vcf_pos"], variant["ref"], variant["alt"]
            elif var_type == "SV_INS":
                length = len(variant["inserted"])
                pos, end, ref, alt = variant["start"], variant["start"], ".", variant["inserted"]
            else:
                length = variant["end"] - variant["start"]
                pos, end, ref, alt = variant["start"], max(variant["end"], variant["start"]), ".", "."
            mate_contig = variant.get("mate_contig", ".")
            mate_pos = variant.get("mate_pos", ".")
            fields = [variant["id"], var_type, variant["origin"], variant["level"], zygosity, variant["contig"], pos, end,
                      ref, alt, length, mate_contig, mate_pos, "%.3f" % normal_vaf, "%.3f" % tumor_vaf,
                      region_tags(variant["contig"], variant["start"], max(variant["end"], variant["start"] + 1))]
            out.write("\t".join(str(f) for f in fields) + "\n")


def write_regions(path):
    with open(path, "w") as out:
        for name, contig, start, end, normal_mult, tumor_mult in COVERAGE_REGIONS:
            out.write("\t".join([contig, str(start), str(end), name, "coverage",
                                 "normal=" + str(normal_mult) + ";tumor=" + str(tumor_mult)]) + "\n")
        for name, contig, start, end in special_regions():
            out.write("\t".join([contig, str(start), str(end), name, "reference", "."]) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument("-s", "--seed", type=int, default=20261009, help="Seed of the simulation")
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.seed)

    donor, published = build_reference(rng)
    write_fasta(os.path.join(args.outdir, "reference.fa"), published)
    small_variants = build_small_variants(rng, donor)
    structural_variants = build_structural_variants(rng)
    write_truth(os.path.join(args.outdir, "truth_variants.tsv"), small_variants, structural_variants)
    write_regions(os.path.join(args.outdir, "regions.bed"))

    # Genomes of the mixture: germline cells and the three nested tumor clones
    germline = build_genome(donor, small_variants, structural_variants, [])
    clones = [build_genome(donor, small_variants, structural_variants, LEVELS[:i + 1]) for i in range(len(LEVELS))]
    sub1 = TUMOR_PURITY * CLONE_CELL_FRACTION["SUB1"]
    sub2 = TUMOR_PURITY * CLONE_CELL_FRACTION["SUB2"]
    tumor_weights = [1.0 - TUMOR_PURITY, TUMOR_PURITY - sub1, sub1 - sub2, sub2]

    normal_stats = simulate_sample(rng, [germline], [1.0], 0, os.path.join(args.outdir, "normal"), "HSIMNDSX2")
    tumor_stats = simulate_sample(rng, [germline] + clones, tumor_weights, 1, os.path.join(args.outdir, "tumor"), "HSIMTDSX3")
    for label, stats in (("normal", normal_stats), ("tumor", tumor_stats)):
        print(label + ": pairs=" + str(stats["pairs"]) + " duplicates=" + str(stats["duplicates"])
              + " substitution_error_rate=%.4f" % (stats["errors"] / float(stats["bases"])))
    print("variants: small=" + str(len(small_variants)) + " structural=" + str(len(structural_variants)))


if __name__ == "__main__":
    main()
