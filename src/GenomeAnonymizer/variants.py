# @author: Nicolas Gaitan
from enum import Enum


class SomaticVariationType(Enum):
    UNCLASSIFIED = 0
    NORMAL_SINGLE_READ_VARIANT = 1
    TUMORAL_SINGLE_READ_VARIANT = 2
    NORMAL_ONLY_VARIANT = 3
    TUMORAL_ONLY_VARIANT = 4
    TUMORAL_NORMAL_VARIANT = 5


class CalledGenomicVariant:
    TYPE_SNV = 'SNV'
    TYPE_INDEL = 'INDEL'

    def __init__(self, seq_name, pos, end, var_type, length, allele):
        self.seq_name: str = seq_name
        self.pos: int = pos
        self.end: int = end
        self.var_type: str = var_type
        self.length = length
        self.allele: str = allele
        self.somatic_variation_type = SomaticVariationType.UNCLASSIFIED

    def __eq__(self, var2):
        if self.seq_name != var2.seq_name:
            return False
        if self.var_type != var2.var_type:
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
