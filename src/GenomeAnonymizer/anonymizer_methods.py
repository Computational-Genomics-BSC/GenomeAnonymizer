# @author: Nicolas Gaitan
import logging
from typing import Protocol, Dict, List, Tuple, Generator, Set
import pysam
from variant_extractor.variants import VariantType
from src.GenomeAnonymizer.variants import CalledGenomicVariant, SomaticVariationType
from src.GenomeAnonymizer.variation_classifier import classify_variation_in_pileup_column, DATASET_IDX_TUMORAL, \
    DATASET_IDX_NORMAL, PAIR_1_IDX, PAIR_2_IDX


class AnonymizedRead:
    def __init__(self, read_alignment: pysam.AlignedSegment, dataset_idx: int):
        self.anonymized_sequence_list = []
        self.anonymized_qualities_list = []
        self.read_id = read_alignment.query_name
        self.is_read1 = read_alignment.is_read1
        self.is_read2 = read_alignment.is_read2
        self.is_reverse = read_alignment.is_reverse
        self.read_alignment = read_alignment
        self.original_sequence = read_alignment.query_sequence
        self.original_qualities = read_alignment.query_qualities
        self.set_original_sequence()
        self.set_original_qualities()
        self.dataset_idx = dataset_idx
        # If AnonymizedReads are initialized correctly, this only happens once if a supp. begins before the primary
        self.has_only_supplementary = read_alignment.is_supplementary
        # Left over variants that were found in the suppl. but the primary mapping was not found yet
        self.left_over_variants_to_mask: List[Tuple[int, str]] = []
        self.has_left_overs_to_mask = False

    def update_from_primary_mapping(self, aln: pysam.AlignedSegment):
        if aln.is_supplementary:
            raise ValueError("Trying to update AnonymizedRead using a supplementary alignment: "
                             "The update should always be called only if the primary mapping appears")
        self.read_alignment = aln
        self.original_sequence = aln.query_sequence
        self.original_qualities = aln.query_qualities
        self.set_original_sequence()
        self.set_original_qualities()
        self.has_only_supplementary = False

    def set_original_sequence(self):
        self.anonymized_sequence_list = list(self.original_sequence)

    def set_original_qualities(self):
        self.anonymized_qualities_list = list(self.original_qualities)

    def mask_or_modify_base_pair(self, pos_in_read: int, new_base: str, modify_qualities: bool = False,
                                 new_quality: int = 0):
        self.anonymized_sequence_list[pos_in_read] = new_base
        if modify_qualities:
            self.anonymized_qualities_list[pos_in_read] = new_quality

    def reverse_complement(self):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        self.anonymized_sequence_list = [complement[base] for base in reversed(self.anonymized_sequence_list)]
        self.anonymized_qualities_list = [base for base in reversed(self.anonymized_qualities_list)]

    def get_anonymized_fastq_record(self) -> pysam.FastxRecord:
        if self.is_reverse:
            self.reverse_complement()
        anonymized_read_seq = ''.join(self.anonymized_sequence_list)
        anonimized_read_qual = ''.join(map(lambda x: chr(x + 33), self.anonymized_qualities_list))
        if len(anonymized_read_seq) != len(self.original_sequence):
            raise ValueError("Anonymized read length does not match original read length")
        if len(anonimized_read_qual) != len(self.original_qualities):
            raise ValueError("Anonymized qualities length does not match original qualities length")
        anonymized_read: pysam.FastxRecord = pysam.FastxRecord(name=self.read_id, sequence=anonymized_read_seq,
                                                               quality=anonimized_read_qual)
        return anonymized_read

    def add_left_over_variant(self, var_pos_in_read, new_base):
        if not self.has_only_supplementary:
            raise ValueError(
                f'Trying to add left over variant to AnonymizedRead {self.read_id} containing a primary mapping'
                f'\n all variants can be masked already')
        self.left_over_variants_to_mask.append((var_pos_in_read, new_base))
        self.has_left_overs_to_mask = True

    def mask_or_anonymize_left_over_variants(self):
        if self.has_only_supplementary:
            raise ValueError(
                f'Trying to mask left over variants in AnonymizedRead {self.read_id} without a primary mapping')
        if not self.has_left_overs_to_mask or len(self.left_over_variants_to_mask) == 0:
            raise ValueError(
                f'Trying to mask left over variants in AnonymizedRead {self.read_id} with no left over variants to mask')
        for var_pos_in_read, new_base in self.left_over_variants_to_mask:
            self.mask_or_modify_base_pair(var_pos_in_read, new_base)
        self.has_left_overs_to_mask = False


class Anonymizer(Protocol):
    """
    Protocol that dictates the methods an anonymizer class must implement.
    """

    def anonymize(self, variant_record, tumor_reads_pileup, normal_reads_pileup, ref_genome) -> None:
        pass

    def mask_germline_snvs(self, pileup_column, called_genomic_variants,
                           variant_record) -> None:
        pass

    def get_anonymized_reads(self):
        pass

    def yield_anonymized_reads(self):
        pass

    def reset(self):
        pass


def decode_specific_read_pair_name(specific_read_pair_name) -> Tuple[str, int]:
    """Decodes the read pair and the read name from the coding applied in the variation classifier"""
    split_name = specific_read_pair_name.split(';')
    read_name = split_name[0]
    pair_number = int(split_name[1])
    return read_name, pair_number


def add_anonymized_read_pair_to_collection(anonymized_reads: Dict[str, List[AnonymizedRead]], aln: pysam.AlignedSegment,
                                           dataset_idx: int):
    if aln.query_name not in anonymized_reads:
        anonymized_reads[aln.query_name] = [None, None]
        paired_anonymized_read_list = anonymized_reads[aln.query_name]
        new_anonymized_read = AnonymizedRead(aln, dataset_idx)
        if aln.is_read1:
            paired_anonymized_read_list[PAIR_1_IDX] = new_anonymized_read
        if aln.is_read2:
            paired_anonymized_read_list[PAIR_2_IDX] = new_anonymized_read
    else:
        paired_anonymized_read_list = anonymized_reads[aln.query_name]
        new_anonymized_read = AnonymizedRead(aln, dataset_idx)
        if aln.is_read1:
            if paired_anonymized_read_list[PAIR_1_IDX] is None:
                paired_anonymized_read_list[PAIR_1_IDX] = new_anonymized_read
            new_anonymized_read = paired_anonymized_read_list[PAIR_1_IDX]
        if aln.is_read2:
            if paired_anonymized_read_list[PAIR_2_IDX] is None:
                paired_anonymized_read_list[PAIR_2_IDX] = new_anonymized_read
            new_anonymized_read = paired_anonymized_read_list[PAIR_2_IDX]
        if not aln.is_supplementary and new_anonymized_read.has_only_supplementary:
            new_anonymized_read.update_from_primary_mapping(aln)


class CompleteGermlineAnonymizer:
    def __init__(self):
        self.anonymized_reads: Dict[str, List[AnonymizedRead]] = dict()
        self.has_anonymized_reads = False

    def reset(self):
        self.anonymized_reads = dict()
        self.has_anonymized_reads = False

    def anonymize(self, variant_record, tumor_reads_pileup, normal_reads_pileup, ref_genome):
        called_genomic_variants = {}
        left_overs_to_mask: Set[AnonymizedRead] = set()
        for dataset_idx, current_pileup in enumerate((tumor_reads_pileup, normal_reads_pileup)):
            seen_read_alns = set()
            for pileup_column in current_pileup:
                classify_variation_in_pileup_column(pileup_column, dataset_idx, seen_read_alns, ref_genome,
                                                    called_genomic_variants)
                # TODO: mask indels
                for pileup_read in pileup_column.pileups:
                    aln = pileup_read.alignment
                    add_anonymized_read_pair_to_collection(self.anonymized_reads, aln, dataset_idx)
                if dataset_idx == DATASET_IDX_NORMAL:
                    pos = pileup_column.reference_pos
                    variants_in_column: List[CalledGenomicVariant] = called_genomic_variants.get(pos, None)
                    if variants_in_column is None:
                        continue
                    self.mask_germline_snvs(variants_in_column, variant_record, left_overs_to_mask)
        # Mask leftovers, referring to reads that were updated from their primary alignment,
        # after initializing with their supplementary
        for anonymized_read in left_overs_to_mask:
            if not anonymized_read.has_only_supplementary:  # and anonymized_read.has_left_overs_to_mask:
                anonymized_read.mask_or_anonymize_left_over_variants()
            else:
                logging.warning(
                    f'Leftover found with only supplementary alignment: {anonymized_read.read_id}'
                    f' Primary mapping may or not be inside of a region to include'
                )
        self.has_anonymized_reads = True

    def mask_germline_snvs(self, variants_in_column, variant_record,
                           left_overs_to_mask: Set[AnonymizedRead]):
        variant_to_keep = CalledGenomicVariant.from_variant_record(variant_record)
        for called_variant in variants_in_column:
            if (called_variant.somatic_variation_type == SomaticVariationType.TUMORAL_NORMAL_VARIANT
                    and called_variant != variant_to_keep):
                # TODO: mask indels
                if called_variant.variant_type == VariantType.SNV:
                    for specific_read_id, var_read_pos in called_variant.supporting_reads.items():
                        read_id, pair = decode_specific_read_pair_name(specific_read_id)
                        anonymized_read = self.anonymized_reads.get(read_id)[pair]
                        if anonymized_read.has_only_supplementary:
                            anonymized_read.add_left_over_variant(var_read_pos, called_variant.ref_allele)
                            left_overs_to_mask.add(anonymized_read)
                            continue
                        anonymized_read.mask_or_modify_base_pair(var_read_pos, called_variant.ref_allele)

    def get_anonymized_reads(self) -> Tuple[List[AnonymizedRead], List[AnonymizedRead]]:
        anonymized_reads_by_dataset = ([], [])
        if not self.has_anonymized_reads:
            raise ValueError("No reads have been anonymized, call anonymize() first")
        for read_id, anonymized_read_obj_pair in self.anonymized_reads.items():
            dataset_records = anonymized_reads_by_dataset[anonymized_read_obj_pair[0].dataset_idx]
            dataset_records.append(anonymized_read_obj_pair)
        return anonymized_reads_by_dataset

    def yield_anonymized_reads(self) -> Generator[Tuple[AnonymizedRead, AnonymizedRead], None, None]:
        """Generator that returns the anonymized reads in pairs (Tuple), indexed by their corresponding pair number:
        0 if pair1, 1 if pair2"""
        if not self.has_anonymized_reads:
            raise ValueError("No reads have been anonymized, call anonymize() first")
        for read_id, anonymized_read_obj_pair in self.anonymized_reads.items():
            if len(anonymized_read_obj_pair) != 2:
                raise ValueError(
                    f'Unexpected number of anonymized reads {len(anonymized_read_obj_pair)} in read id: {read_id}')
            yield anonymized_read_obj_pair
