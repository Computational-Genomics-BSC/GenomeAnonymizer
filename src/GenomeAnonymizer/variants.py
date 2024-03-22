# @author: Nicolas Gaitan
from enum import Enum
from variant_extractor.variants import VariantRecord
from variant_extractor.variants import VariantType


class SomaticVariationType(Enum):
    UNCLASSIFIED = 0
    NORMAL_SINGLE_READ_VARIANT = 1
    TUMORAL_SINGLE_READ_VARIANT = 2
    NORMAL_ONLY_VARIANT = 3
    TUMORAL_ONLY_VARIANT = 4
    TUMORAL_NORMAL_VARIANT = 5


class CalledGenomicVariant:

    def __init__(self, seq_name, pos, end, var_type, length, allele):
        self.seq_name: str = seq_name
        self.pos: int = pos
        self.end: int = end
        self.variant_type: VariantType = var_type
        self.length = length
        self.allele: str = allele
        self.somatic_variation_type = SomaticVariationType.UNCLASSIFIED
        # Dictionary with the supporting read ids as key and the position of the variant in the read as value
        self.supporting_reads = dict()

    @classmethod
    def from_variant_record(cls, variant_record: VariantRecord):
        return cls(variant_record.contig, variant_record.pos, variant_record.end,
                   variant_record.variant_type, variant_record.length, variant_record.alt)

    def add_supporting_read(self, read_id, var_read_pos):
        self.supporting_reads[read_id] = var_read_pos

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
