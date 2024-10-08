# @author: Nicolas Gaitan
import math
from enum import Enum

from variant_extractor.variants import VariantRecord
from variant_extractor.variants import VariantType


def compare(seq_idx1: int, first1: int, last1: int, seq_idx2: int, first2: int, last2: int) -> int:
    overlap = first2 <= last1 and last2 >= first1
    if seq_idx1 < seq_idx2:
        return -3
    if seq_idx1 > seq_idx2:
        return 3
    # For these cases seq_idx1 == seq_idx2
    if last1 < last2:
        return -1 if overlap else -2
    if last2 < last1:
        return 1 if overlap else 2
    # For these cases last1 == last2
    if first1 < first2:
        return -1
    if first2 < first1:
        return 1
    return 0


def estimate_euclidean_distance(x1: int, y1: int, z1: int, x2: int, y2: int, z2: int) -> float:
    square_sum = (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
    return math.sqrt(square_sum)


class SomaticVariationType(Enum):
    UNCLASSIFIED = 0
    NORMAL_SINGLE_READ_VARIANT = 1
    TUMORAL_SINGLE_READ_VARIANT = 2
    NORMAL_ONLY_VARIANT = 3
    TUMORAL_ONLY_VARIANT = 4
    TUMORAL_NORMAL_VARIANT = 5


class CalledGenomicVariant:

    def __init__(self, seq_name, pos, end, var_type, length, allele, ref_allele):
        self.seq_name: str = seq_name
        self.pos: int = pos
        self.end: int = end
        self.variant_type: VariantType = var_type
        self.length = length
        self.allele: str = allele
        self.ref_allele = ref_allele
        self.somatic_variation_type = SomaticVariationType.UNCLASSIFIED
        # self.has_diffused = False
        self.is_linked_to_another_germline = False
        # Dictionary with the supporting read ids as key and the position of the variant in the read as value
        self.supporting_reads = dict()

    @classmethod
    def from_variant_record(cls, variant_record: VariantRecord):
        # VariantRecord coordinates are 1-based, while CalledGenomicVariant coordinates are 0-based
        return cls(variant_record.contig, variant_record.pos - 1, variant_record.end - 1,
                   variant_record.variant_type, variant_record.length, variant_record.alt, variant_record.ref)

    def add_supporting_read(self, read_id, var_read_pos):
        self.supporting_reads[read_id] = var_read_pos

    # def set_diffused_status(self):
    #     self.has_diffused = True

    def set_link_to_another_germline(self):
        self.is_linked_to_another_germline = True

    def is_candidate_for_diffusion(self):
        # if self.somatic_variation_type == SomaticVariationType.TUMORAL_NORMAL_VARIANT:
        #     return False
        if self.is_linked_to_another_germline:
            return False
        return True

    def calculate_distance_to_another(self, variant2):
        return estimate_euclidean_distance(self.pos, self.end, self.length, variant2.pos, variant2.end, variant2.length)

    def __eq__(self, var2):
        if self.seq_name != var2.seq_name:
            return False
        if self.variant_type != var2.variant_type:
            return False
        if self.pos != var2.pos:
            return False
        if self.end != var2.end:
            return False
        if self.length != var2.length:
            return False
        if self.allele != var2.allele:
            return False
        return True

    def __str__(self):
        return (f'seq_name: {self.seq_name} pos: {self.pos} end: {self.end} var_type: {self.variant_type} '
                f'length: {self.length} alt_allele: {self.allele} ref_allele: {self.ref_allele} '
                f'somatic_variation_type: {self.somatic_variation_type}')
