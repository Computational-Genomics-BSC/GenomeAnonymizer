# @author: Nicolas Gaitan

from typing import Protocol
from variants import CalledGenomicVariant, SomaticVariationType
from variation_classifier import classify_variation_in_pileup_column


class Anonymizer(Protocol):
    """
    Interface that dictates the methods an anonymizer class must follow.

    Args:
        variant_record: The variant not to be anonymized.
        normal_reads_pileup: The pileup of normal reads.
        tumor_reads_pileup: The pileup of tumor reads.

    Returns:
        None
    """

    def anonymize(self, variant_record, tumor_reads_pileup, normal_reads_pileup) -> None:
        pass


class CompleteGermlineAnonymizer:
    def __init__(self):
        self.anonymized_reads = []

    def anonymize(self, variant_record, tumor_reads_pileup, normal_reads_pileup):
        called_genomic_variants = {}
        seen_read_alns = set()
        for dataset_idx, current_pileup in enumerate((tumor_reads_pileup, normal_reads_pileup)):
            for pileup_column in current_pileup:
                classify_variation_in_pileup_column(pileup_column, dataset_idx, seen_read_alns, called_genomic_variants)
                self.mask_germline_snvs(pileup_column, called_genomic_variants, variant_record)

    def mask_germline_snvs(self, pileup_column, called_genomic_variants, variant_record):
        pos = pileup_column.reference_pos
        variants_in_column = called_genomic_variants.get(pos, None)
        if variants_in_column is None:
            return
        for called_variant in variants_in_column:
            # TODO: mask indels
            if (called_variant.var_type == CalledGenomicVariant.TYPE_SNV and
                    called_variant.somatic_variation_type == SomaticVariationType.TUMORAL_NORMAL_VARIANT):
                pass

