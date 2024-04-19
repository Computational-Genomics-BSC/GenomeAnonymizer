# @author: Nicolas Gaitan
import logging
import re
from typing import List
from timeit import default_timer as timer
from pysam import PileupColumn, FastaFile, AlignedSegment, PileupRead
from variant_extractor.variants import VariantType
from src.GenomeAnonymizer.variants import CalledGenomicVariant, SomaticVariationType

# constants
DATASET_IDX_TUMORAL = 0
DATASET_IDX_NORMAL = 1

PAIR_1_IDX = 0
PAIR_2_IDX = 1

DEBUG_TOTAL_TIMES = {'anonymize_windows': 0, 'anonymize_call': 0, 'write_pairs': 0,
                     'unpaired_searches': 0, 'process_indels': 0, 'process_snvs': 0,
                     'mask_germlines': 0, 'mask_germline_snvs': 0, 'mask_germlines_left_overs_in_window': 0, 'classify_variants': 0}


def generate_pair_name(aln):
    return f'{aln.query_name};{PAIR_1_IDX}' if aln.is_read1 else f'{aln.query_name};{PAIR_2_IDX}'


def process_indels(aln: AlignedSegment, specific_pair_query_name, dataset_idx, ref_genome, called_genomic_variants):
    regexp = r"(?<=[a-zA-Z=])(?=[0-9])|(?<=[0-9])(?=[a-zA-Z=])"  # regex to split cigar string
    cigar_indels = {"I", "D"}  # cigar operations to report
    ref_consuming = {'M', 'D', 'N', '=', 'X'}  # stores reference consuming cigar operations
    read_consuming_only = ['S', 'H', 'I']  # stores read consuming cigar operations
    cigar_list = re.split(regexp, aln.cigarstring)
    start_ref_pos = aln.reference_start
    current_cigar_len = 0
    read_consumed_bases = 0
    seq_name = aln.reference_name
    for cigar_list_idx, symbol in enumerate(cigar_list):
        if symbol.isdigit():
            cigar_op = cigar_list[cigar_list_idx + 1]
            if cigar_op in cigar_indels:
                pos = start_ref_pos + current_cigar_len
                in_read_pos = current_cigar_len + read_consumed_bases
                length = int(symbol)
                var_type = VariantType.INS if cigar_op == 'I' else VariantType.DEL
                end = pos + 1 if var_type == VariantType.INS else pos + length
                # TODO: Fix alleles for indel masking
                called_indel = CalledGenomicVariant(seq_name, pos, end, var_type, length, "", "")
                if called_indel.pos not in called_genomic_variants:
                    called_genomic_variants[called_indel.pos] = []
                indel_pos_list = called_genomic_variants[called_indel.pos]
                indel_exists = False
                for var_indel in indel_pos_list:
                    if called_indel == var_indel:
                        called_indel = var_indel
                        indel_exists = True
                        break
                called_indel.add_supporting_read(specific_pair_query_name, in_read_pos)
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
            if cigar_op in read_consuming_only:
                read_consumed_bases += int(symbol)
            if cigar_op == 'D':
                read_consumed_bases -= int(symbol)


def process_snvs(aln: AlignedSegment, specific_pair_query_name, reference_pos, ref_base, in_read_position, dataset_idx,
                 called_genomic_variants):
    seq_name = aln.reference_name
    ignore_bases = {'N'}
    base = aln.query_sequence[in_read_position].upper()
    if base in ignore_bases or base == ref_base:
        return
    called_snv = CalledGenomicVariant(seq_name, reference_pos, reference_pos, VariantType.SNV, 1,
                                      allele=base, ref_allele=ref_base)
    if called_snv.pos not in called_genomic_variants:
        called_genomic_variants[called_snv.pos] = []
    snv_pos_list = called_genomic_variants[called_snv.pos]
    snv_exists = False
    for var_snv in snv_pos_list:
        if called_snv == var_snv:
            called_snv = var_snv
            snv_exists = True
            break
    called_snv.add_supporting_read(specific_pair_query_name, in_read_position)
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


def classify_variation_in_pileup_column(pileup_column: PileupColumn, dataset_idx, seen_read_alns, ref_genome: FastaFile,
                                        called_genomic_variants):
    """
    Classify the read variation returning a dictionary with every INDEL and SNV CalledGenomicVariant by coordinate.
    """
    pileups: List[PileupRead] = pileup_column.pileups
    reference_pos = pileup_column.reference_pos
    ref_base = \
    ref_genome.fetch(pileup_column.reference_name, pileup_column.reference_pos, pileup_column.reference_pos + 1)[0]
    ref_base = ref_base.upper()
    for pileup_read in pileups:
        aln: AlignedSegment = pileup_read.alignment
        specific_pair_query_name = generate_pair_name(aln)
        if specific_pair_query_name not in seen_read_alns:
            start1 = timer()
            process_indels(aln, specific_pair_query_name, dataset_idx, ref_genome, called_genomic_variants)
            end1 = timer()
            logging.debug("Time to process indels: " + str(end1 - start1))
            DEBUG_TOTAL_TIMES['process_indels'] += end1 - start1
            seen_read_alns.add(specific_pair_query_name)
        in_read_position = pileup_read.query_position
        if in_read_position is None:
            continue
        start2 = timer()
        process_snvs(aln, specific_pair_query_name, reference_pos, ref_base, in_read_position, dataset_idx,
                     called_genomic_variants)
        end2 = timer()
        logging.debug("Time to process SNVs: " + str(end2 - start2))
        DEBUG_TOTAL_TIMES['process_snvs'] += end2 - start2
