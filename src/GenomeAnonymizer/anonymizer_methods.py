# @author: Nicolas Gaitan
from typing import Protocol
from variation_classifier import classify_variation


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

    def anonymize(self, variant_record, normal_reads_pileup, tumor_reads_pileup) -> None:
        pass


class CompleteGermlineAnonymizer:
    def anonymize(self, variant_record, normal_reads_pileup, tumor_reads_pileup):
        classify_variation(variant_record, normal_reads_pileup, tumor_reads_pileup)
        pass
