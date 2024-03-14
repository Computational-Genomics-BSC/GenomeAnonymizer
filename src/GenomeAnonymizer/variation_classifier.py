import re
from typing import Dict, List
import pysam
from variant_extractor import VariantExtractor
from variants import CalledGenomicVariant, SomaticVariationType


def process_indels(pileups, dataset_idx, seen_read_alns, called_genomic_variants):
    regexp = r"(?<=[a-zA-Z=])(?=[0-9])|(?<=[0-9])(?=[a-zA-Z=])"  # regex to split cigar string
    cigar_indels = {"I", "D"}  # cigar operations to report
    ref_consuming = {'M', 'D', 'N', '=', 'X'}  # stores reference consuming cigar operations
    for pileup_read in pileups:
        aln = pileup_read.alignment
        if aln.query_name not in seen_read_alns:
            return
        seen_read_alns.add(aln.query_name)
        cigar_list = re.split(regexp, aln.cigarstring)
        start_ref_pos = aln.reference_start
        current_cigar_len = 0
        seq_name = aln.reference_name
        for cigar_list_idx, symbol in enumerate(cigar_list):
            if symbol.isdigit():
                cigar_op = cigar_list[cigar_list_idx + 1]
                if cigar_op in cigar_indels:
                    pos = start_ref_pos + current_cigar_len
                    length = int(symbol)
                    var_type = CalledGenomicVariant.TYPE_INDEL
                    end = pos + 1 if symbol == 'I' else pos + length
                    called_indel = CalledGenomicVariant(seq_name, pos, end, var_type, length, "")
                    if called_indel.pos not in called_genomic_variants:
                        called_genomic_variants[called_indel.pos] = []
                    indel_pos_list = called_genomic_variants[called_indel.pos]
                    indel_exists = False
                    for var_indel in indel_pos_list:
                        if called_indel == var_indel:
                            called_indel = var_indel
                            indel_exists = True
                            break
                    if not indel_exists:
                        if dataset_idx == 0:
                            called_indel.somatic_variation_type = SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT
                        if dataset_idx == 1:
                            called_indel.somatic_variation_type = SomaticVariationType.NORMAL_SINGLE_READ_VARIANT
                        called_genomic_variants[called_indel.pos].append(called_indel)
                    else:  # if the variation is already annotated in the dict
                        var_code = called_indel.somatic_variation_type
                        if dataset_idx == 0:
                            if var_code == SomaticVariationType.NORMAL_SINGLE_READ_VARIANT or var_code == SomaticVariationType.NORMAL_ONLY_VARIANT:
                                called_indel.somatic_variation_type = SomaticVariationType.TUMORAL_NORMAL_VARIANT
                            if var_code == SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT:
                                called_indel.somatic_variation_type = SomaticVariationType.TUMORAL_ONLY_VARIANT
                        if dataset_idx == 1:
                            if var_code == SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT or var_code == SomaticVariationType.TUMORAL_ONLY_VARIANT:
                                called_indel.somatic_variation_type = SomaticVariationType.TUMORAL_NORMAL_VARIANT
                            if var_code == SomaticVariationType.NORMAL_SINGLE_READ_VARIANT:
                                called_indel.somatic_variation_type = SomaticVariationType.NORMAL_ONLY_VARIANT
                if cigar_op in ref_consuming:
                    current_cigar_len += int(symbol)


def classify_variation(variant_record, tumor_reads_pileup, normal_reads_pileup) -> Dict[str, List[CalledGenomicVariant]]:
    """
    Classify the read variation returning a dictionary with every INDEL and SNV CalledGenomicVariant by coordinate.
    """
    called_genomic_variants = {}
    seen_read_alns = set()
    for dataset_idx, current_pileup in enumerate((tumor_reads_pileup, normal_reads_pileup)):
        for pileup_column in current_pileup:
            process_indels(pileup_column.pileups, dataset_idx, seen_read_alns, called_genomic_variants)
    return called_genomic_variants
