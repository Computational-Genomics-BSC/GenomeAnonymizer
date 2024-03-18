# @author: Nicolas Gaitan

from typing import Protocol

import pysam
from variants import CalledGenomicVariant, SomaticVariationType
from variant_extractor.variants import VariantRecord
from variant_extractor.variants import VariantType
from variation_classifier import classify_variation_in_pileup_column
from variation_classifier import DATASET_IDX_NORMAL, DATASET_IDX_TUMORAL


class Anonymizer(Protocol):
    """
    Protocol that dictates the methods an anonymizer class must implement.

    Args:
        variant_record: The variant not to be anonymized.
        normal_reads_pileup: The pileup of normal reads.
        tumor_reads_pileup: The pileup of tumor reads.

    Returns:
        None
    """

    def anonymize(self, variant_record, tumor_reads_pileup, normal_reads_pileup) -> None:
        pass

    def mask_germline_snvs(self, pileup_column, called_genomic_variants, variant_record) -> None:
        pass


class CompleteGermlineAnonymizer:
    def __init__(self):
        self.anonymized_reads = dict()

    def anonymize(self, variant_record, tumor_reads_pileup, normal_reads_pileup):
        called_genomic_variants = {}
        seen_read_alns = set()
        for dataset_idx, current_pileup in enumerate((tumor_reads_pileup, normal_reads_pileup)):
            for pileup_column in current_pileup:
                classify_variation_in_pileup_column(pileup_column, dataset_idx, seen_read_alns, called_genomic_variants)
                self.mask_germline_snvs(pileup_column, dataset_idx, called_genomic_variants, variant_record)

    def mask_germline_snvs(self, pileup_column, dataset_idx, called_genomic_variants, variant_record):
        # TODO: mask indels
        pos = pileup_column.reference_pos
        for pileup_read in pileup_column.pileups:
            aln = pileup_read.alignment
            if aln.query_name not in self.anonymized_reads:
                self.anonymized_reads[aln.query_name] = AnonymizedRead(aln, dataset_idx)
        variants_in_column = called_genomic_variants.get(pos, None)
        if variants_in_column is None:
            # TODO: Still collect reads without variation
            return
        variant_to_keep = CalledGenomicVariant.from_variant_record(variant_record)
        for called_variant in variants_in_column:
            if (called_variant.variant_type == VariantType.SNV and
                    called_variant.somatic_variation_type == SomaticVariationType.TUMORAL_NORMAL_VARIANT
                    and called_variant != variant_to_keep):
                for read_id in called_variant.supporting_read_ids:
                    # TODO: modify reads according to the variation
                    pass
                pass


class AnonymizedRead:
    def __init__(self, read_alignment: pysam.AlignedSegment, dataset_idx: int):
        self.read_id = read_alignment.query_name
        self.read_alignment = read_alignment
        self.original_sequence = read_alignment.query_sequence
        self.original_qualities = read_alignment.query_qualities
        self.anonymized_sequence_list = list(self.original_sequence)
        self.anonymized_qualities_list = list(self.original_qualities)
        self.dataset_idx = dataset_idx

    def modify_base_pair(self, pos_in_read: int, new_base: str, modify_qualities: bool = False, new_quality: int = 0):
        self.anonymized_sequence_list[pos_in_read] = new_base
        if modify_qualities:
            self.anonymized_qualities_list[pos_in_read] = new_quality

    def get_anonymized_read(self) -> pysam.FastxRecord:
        anonymized_read_seq = ''.join(self.anonymized_sequence_list)
        anonimized_read_qual = ''.join(self.anonymized_qualities_list)
        if len(anonymized_read_seq) != len(self.original_sequence):
            raise ValueError("Anonymized read length does not match original read length")
        if len(anonimized_read_qual) != len(self.original_qualities):
            raise ValueError("Anonymized qualities length does not match original qualities length")
        anonymized_read: pysam.FastxRecord = pysam.FastxRecord(name=self.read_id, sequence=anonymized_read_seq,
                                                               quality=anonimized_read_qual)
        return anonymized_read
