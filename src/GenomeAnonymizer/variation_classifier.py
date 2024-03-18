# @author: Nicolas Gaitan

import re
from typing import Dict, List
import pysam
from variant_extractor import VariantExtractor
from variant_extractor.variants import VariantType
from variants import CalledGenomicVariant, SomaticVariationType


# constants
DATASET_IDX_TUMORAL = 0
DATASET_IDX_NORMAL = 1


def process_indels(pileup_column, dataset_idx, seen_read_alns, called_genomic_variants):
    regexp = r"(?<=[a-zA-Z=])(?=[0-9])|(?<=[0-9])(?=[a-zA-Z=])"  # regex to split cigar string
    cigar_indels = {"I", "D"}  # cigar operations to report
    ref_consuming = {'M', 'D', 'N', '=', 'X'}  # stores reference consuming cigar operations
    pileups = pileup_column.pileups
    for pileup_read in pileups:
        aln = pileup_read.alignment
        if aln.query_name in seen_read_alns:
            continue
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
                    var_type = VariantType.INS if cigar_op == 'I' else VariantType.DEL
                    end = pos + 1 if var_type == VariantType.INS else pos + length
                    called_indel = CalledGenomicVariant(seq_name, pos, end, var_type, length, "")
                    called_indel.add_supporting_read_id(aln.query_name)
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
                        if dataset_idx == DATASET_IDX_TUMORAL:
                            called_indel.somatic_variation_type = SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT
                        if dataset_idx == DATASET_IDX_NORMAL:
                            called_indel.somatic_variation_type = SomaticVariationType.NORMAL_SINGLE_READ_VARIANT
                        called_genomic_variants[called_indel.pos].append(called_indel)
                    else:  # if the variation is already annotated in the dict
                        var_code = called_indel.somatic_variation_type
                        if dataset_idx == DATASET_IDX_TUMORAL:
                            if (var_code == SomaticVariationType.NORMAL_SINGLE_READ_VARIANT or
                                    var_code == SomaticVariationType.NORMAL_ONLY_VARIANT):
                                called_indel.somatic_variation_type = SomaticVariationType.TUMORAL_NORMAL_VARIANT
                            if var_code == SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT:
                                called_indel.somatic_variation_type = SomaticVariationType.TUMORAL_ONLY_VARIANT
                        if dataset_idx == DATASET_IDX_NORMAL:
                            if (var_code == SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT or
                                    var_code == SomaticVariationType.TUMORAL_ONLY_VARIANT):
                                called_indel.somatic_variation_type = SomaticVariationType.TUMORAL_NORMAL_VARIANT
                            if var_code == SomaticVariationType.NORMAL_SINGLE_READ_VARIANT:
                                called_indel.somatic_variation_type = SomaticVariationType.NORMAL_ONLY_VARIANT
                if cigar_op in ref_consuming:
                    current_cigar_len += int(symbol)


def process_snvs(pileup_column, dataset_idx, called_genomic_variants):
    column_bases = pileup_column.get_query_sequences(mark_matches=True)
    seq_name = pileup_column.reference_name
    column_pos = pileup_column.reference_pos
    ignore_bases = {'.', ',', 'N'}
    read_names = pileup_column.get_query_names
    for i, base in enumerate(column_bases):
        base = base.upper()
        if base in ignore_bases:
            continue
        called_snv = CalledGenomicVariant(seq_name, column_pos, column_pos, VariantType.SNV, 1,
                                          base)
        called_snv.add_supporting_read_id(read_names[i])
        if called_snv.pos not in called_genomic_variants:
            called_genomic_variants[called_snv.pos] = []
        snv_pos_list = called_genomic_variants[called_snv.pos]
        snv_exists = False
        for var_snv in snv_pos_list:
            if called_snv == var_snv:
                called_snv = var_snv
                snv_exists = True
                break
        if not snv_exists:
            if dataset_idx == DATASET_IDX_TUMORAL:
                called_snv.somatic_variation_type = SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT
            if dataset_idx == DATASET_IDX_NORMAL:
                called_snv.somatic_variation_type = SomaticVariationType.NORMAL_SINGLE_READ_VARIANT
            snv_pos_list.append(called_snv)
        else:
            var_code = called_snv.somatic_variation_type
            if dataset_idx == DATASET_IDX_TUMORAL:
                if (var_code == SomaticVariationType.NORMAL_SINGLE_READ_VARIANT or
                        var_code == SomaticVariationType.NORMAL_ONLY_VARIANT):
                    called_snv.somatic_variation_type = SomaticVariationType.TUMORAL_NORMAL_VARIANT
                if var_code == SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT:
                    called_snv.somatic_variation_type = SomaticVariationType.TUMORAL_ONLY_VARIANT
            if dataset_idx == DATASET_IDX_NORMAL:
                if (var_code == SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT or
                        var_code == SomaticVariationType.TUMORAL_ONLY_VARIANT):
                    called_snv.somatic_variation_type = SomaticVariationType.TUMORAL_NORMAL_VARIANT
                if var_code == SomaticVariationType.NORMAL_SINGLE_READ_VARIANT:
                    called_snv.somatic_variation_type = SomaticVariationType.NORMAL_ONLY_VARIANT


def classify_variation_in_pileup_column(pileup_column, dataset_idx, seen_read_alns, called_genomic_variants):
    """
    Classify the read variation returning a dictionary with every INDEL and SNV CalledGenomicVariant by coordinate.
    """
    process_indels(pileup_column, dataset_idx, seen_read_alns, called_genomic_variants)
    process_snvs(pileup_column, dataset_idx, called_genomic_variants)
