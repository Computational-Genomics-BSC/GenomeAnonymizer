# @author: Nicolas Gaitan
import logging
import sys
from typing import Protocol, Dict, List, Tuple, Generator, Set, Union
import itertools as it
import numpy as np
import pysam
from timeit import default_timer as timer

from pysam import PileupColumn
from variant_extractor.variants import VariantType
from src.GenomeAnonymizer.variants import CalledGenomicVariant, SomaticVariationType
from src.GenomeAnonymizer.variation_classifier import classify_variation_in_pileup_column, DATASET_IDX_TUMORAL, \
    DATASET_IDX_NORMAL, PAIR_1_IDX, PAIR_2_IDX, DEBUG_TOTAL_TIMES

#  = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
reverses = {ord('A'): ord('T'), ord('C'): ord('G'), ord('G'): ord('C'), ord('T'): ord('A'), ord('N'): ord('N')}
INT_TO_BASE = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'N'}


# def encode_base(base: str) -> int:
#    return BASE_TO_INT[base]

def encode_base(base: str):
    return bytearray(base.encode())[0]


def decode_base(int_base: int) -> str:
    return chr(int_base)


def get_pileup_pair_in_order(pileup1, pileup2) -> Generator[Tuple[PileupColumn, PileupColumn], None, None]:
    p1: PileupColumn = next(pileup1, None)
    p2: PileupColumn = next(pileup2, None)
    while p1 is not None and p2 is not None:
        if p1.reference_pos < p2.reference_pos or p2 is None:
            yield p1, None
            p1 = next(pileup1, None)
        elif p1.reference_pos > p2.reference_pos or p1 is None:
            yield None, p2
            p2 = next(pileup2, None)
        else:
            yield p1, p2
            p1 = next(pileup1, None)
            p2 = next(pileup2, None)


class AnonymizedRead:
    def __init__(self, read_alignment: pysam.AlignedSegment, dataset_idx: int):
        self.anonymized_sequence_array: np.ndarray = []
        self.anonymized_qualities_array: np.ndarray = []
        self.query_name = read_alignment.query_name
        self.is_read1 = read_alignment.is_read1
        self.is_read2 = read_alignment.is_read2
        self.is_reverse = read_alignment.is_reverse
        # self.read_alignment = read_alignment
        self.set_original_sequence(read_alignment.query_sequence)
        self.set_original_qualities(read_alignment.query_qualities)
        self.dataset_idx = dataset_idx
        # If AnonymizedReads are initialized correctly, this only happens once if a supp. begins before the primary
        self.is_supplementary = read_alignment.is_supplementary
        # Left over variants that were found in the suppl. but the primary mapping was not found yet
        self.left_over_variants_to_mask: List[Tuple[int, int]] = []
        self.has_left_overs_to_mask = False

    def get_pair_idx(self):
        if self.is_read1:
            return PAIR_1_IDX
        if self.is_read2:
            return PAIR_2_IDX

    def update_from_primary_mapping(self, aln: pysam.AlignedSegment):
        if aln.is_supplementary:
            raise ValueError("Trying to update AnonymizedRead using a supplementary alignment: "
                             "The update should always be called only if the primary mapping appears")
        # self.read_alignment = aln
        self.set_original_sequence(aln.query_sequence)
        self.set_original_qualities(aln.query_qualities)
        self.is_supplementary = False

    def set_original_sequence(self, original_sequence: str):
        # self.anonymized_sequence_array: List[int] = [encode_base(base) for base in original_sequence]
        self.anonymized_sequence_array: np.ndarray = np.frombuffer(bytearray(original_sequence.encode()),
                                                                   dtype=np.uint8)

    def set_original_qualities(self, original_qualities):
        # self.anonymized_qualities_array = list(original_qualities)
        self.anonymized_qualities_array: np.ndarray = np.frombuffer(original_qualities,
                                                                    dtype=np.uint8)

    def mask_or_modify_base_pair(self, pos_in_read: int, new_base: str, modify_qualities: bool = False,
                                 new_quality: int = 0):
        # self.anonymized_sequence_array[pos_in_read] = encode_base(new_base)
        np.put(self.anonymized_sequence_array, pos_in_read, encode_base(new_base), mode='raise')
        if modify_qualities:
            np.put(self.anonymized_qualities_array, pos_in_read, new_quality, mode='raise')

    def reverse_complement(self):
        # complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        # complement = {0: 3, 1: 2, 2: 1, 3: 0, 4: 4}
        self.anonymized_sequence_array = np.flip(np.vectorize(reverses.get)(self.anonymized_sequence_array))
        # self.anonymized_sequence_array = [encode_base(complement[decode_base(base)]) for base in
        #                                  reversed(self.anonymized_sequence_array)]
        self.anonymized_qualities_array = np.flip(self.anonymized_qualities_array)

    def get_anonymized_fastq_record(self) -> pysam.FastxRecord:
        if self.is_reverse:
            self.reverse_complement()
        # anonymized_read_seq = ''.join(map(lambda x: decode_base(x), self.anonymized_sequence_array))
        anonymized_read_seq: str = np.array2string(self.anonymized_sequence_array,
                                                   formatter={'int': lambda x: decode_base(x)},
                                                   separator='',
                                                   max_line_width=sys.maxsize).strip('[]')
        # anonimized_read_qual: str = ''.join(map(lambda x: chr(x + 33), self.anonymized_qualities_array))
        # anonimized_read_qual: str = np.array2string(self.anonymized_qualities_array,
        #                                            formatter={'int': lambda x: chr(x + 33)},
        #                                            separator='',
        #                                            max_line_width=sys.maxsize).strip('[]')
        anonimized_read_qual: str = ''.join([chr(x + 33) for x in self.anonymized_qualities_array])
        # if len(anonymized_read_seq) != len(self.original_sequence):
        #    raise ValueError("Anonymized read length does not match original read length")
        # if len(anonimized_read_qual) != len(self.original_qualities):
        #    raise ValueError("Anonymized qualities length does not match original qualities length")
        read_pair_name = f'{self.query_name}/{PAIR_1_IDX + 1}' if self.is_read1 else f'{self.query_name}/{PAIR_2_IDX + 1}'
        anonymized_read: pysam.FastxRecord = pysam.FastxRecord(name=read_pair_name, sequence=anonymized_read_seq,
                                                               quality=anonimized_read_qual)
        return anonymized_read

    def add_left_over_variant(self, var_pos_in_read, new_base):
        if not self.is_supplementary:
            raise ValueError(
                f'Trying to add left over variant to AnonymizedRead {self.query_name} containing a primary mapping'
                f'\n all variants can be masked already')
        self.left_over_variants_to_mask.append((var_pos_in_read, ord(new_base)))
        self.has_left_overs_to_mask = True

    def mask_or_anonymize_left_over_variants(self):
        if self.is_supplementary:
            raise ValueError(
                f'Trying to mask left over variants in AnonymizedRead {self.query_name} without a primary mapping')
        if not self.has_left_overs_to_mask or len(self.left_over_variants_to_mask) == 0:
            logging.warning(
                f'Trying to mask left over variants in AnonymizedRead {self.query_name} with no left over variants to mask')
        for var_pos_in_read, new_base in self.left_over_variants_to_mask:
            self.mask_or_modify_base_pair(var_pos_in_read, decode_base(new_base))
        self.has_left_overs_to_mask = False

    def __hash__(self):
        return hash((self.query_name, self.get_pair_idx()))


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


def add_anonymized_read_pair_to_collection_from_alignment(anonymized_reads: Dict[str, List[AnonymizedRead]],
                                                          aln: pysam.AlignedSegment,
                                                          dataset_idx: int):
    if aln.query_name not in anonymized_reads:
        anonymized_reads[aln.query_name] = [None, None]
        paired_anonymized_read_list = anonymized_reads[aln.query_name]
        new_anonymized_read = AnonymizedRead(aln, dataset_idx)
        pair_idx = new_anonymized_read.get_pair_idx()
        paired_anonymized_read_list[pair_idx] = new_anonymized_read
    else:
        paired_anonymized_read_list = anonymized_reads[aln.query_name]
        new_anonymized_read = AnonymizedRead(aln, dataset_idx)
        pair_idx = new_anonymized_read.get_pair_idx()
        if paired_anonymized_read_list[pair_idx] is None:
            paired_anonymized_read_list[pair_idx] = new_anonymized_read
        new_anonymized_read = paired_anonymized_read_list[pair_idx]
        # This only happens if the aln is an AlignedSegment, is not a supplementary mapping,
        # and the new_anonymized_read does not come from a primary mapping
        if not aln.is_supplementary and new_anonymized_read.is_supplementary:
            new_anonymized_read.update_from_primary_mapping(aln)


def add_or_update_anonymized_read_from_other(anonymized_reads: Dict[str, List[AnonymizedRead]],
                                             anonymized_read: AnonymizedRead):
    if anonymized_read.query_name not in anonymized_reads:
        anonymized_reads[anonymized_read.query_name] = [None, None]
        paired_anonymized_read_list = anonymized_reads[anonymized_read.query_name]
        pair_idx = anonymized_read.get_pair_idx()
        paired_anonymized_read_list[pair_idx] = anonymized_read
    else:
        paired_anonymized_read_list = anonymized_reads[anonymized_read.query_name]
        pair_idx = anonymized_read.get_pair_idx()
        if paired_anonymized_read_list[pair_idx] is None:
            paired_anonymized_read_list[pair_idx] = anonymized_read
            return
        saved_anonymized_read = paired_anonymized_read_list[pair_idx]
        if saved_anonymized_read.is_supplementary and not anonymized_read.is_supplementary:
            anonymized_read.left_over_variants_to_mask.extend(saved_anonymized_read.left_over_variants_to_mask)
            if len(anonymized_read.left_over_variants_to_mask) > 0:
                anonymized_read.has_left_overs_to_mask = True
            paired_anonymized_read_list[pair_idx] = anonymized_read
            return
        if anonymized_read.has_left_overs_to_mask:
            saved_anonymized_read.left_over_variants_to_mask.extend(anonymized_read.left_over_variants_to_mask)
            if len(saved_anonymized_read.left_over_variants_to_mask) > 0:
                saved_anonymized_read.has_left_overs_to_mask = True


def anonymized_read_pair_is_writeable(anonymized_read_pair1: AnonymizedRead,
                                      anonymized_read_pair2: AnonymizedRead) -> bool:
    if anonymized_read_pair1 is None or anonymized_read_pair2 is None:
        return False
    if anonymized_read_pair1.is_supplementary or anonymized_read_pair2.is_supplementary:
        return False
    return True


def is_candidate_to_yield(aln, pile_up_column_ref_pos):
    # if aln.has_tag('SA'):
    # if aln.get_tag('SA') != '':
    #     if aln.query_name == 'HWI-ST1324:58:D1D2VACXX:6:2212:3144:42486':
    #         print('SA tag found')
    #     return False
    if aln.reference_end != pile_up_column_ref_pos:
        return False
    return True


def mask_left_over_variants_in_pair(anonymized_read_pair1: AnonymizedRead, anonymized_read_pair2: AnonymizedRead):
    if anonymized_read_pair1 is not None:
        if not anonymized_read_pair1.is_supplementary and anonymized_read_pair1.has_left_overs_to_mask:
            anonymized_read_pair1.mask_or_anonymize_left_over_variants()
    if anonymized_read_pair2 is not None:
        if not anonymized_read_pair2.is_supplementary and anonymized_read_pair2.has_left_overs_to_mask:
            anonymized_read_pair2.mask_or_anonymize_left_over_variants()


class CompleteGermlineAnonymizer:
    def __init__(self):
        self.anonymized_reads: Dict[str, List[AnonymizedRead]] = dict()
        # self.has_anonymized_reads = False

    def reset(self):
        self.anonymized_reads = dict()
        # self.has_anonymized_reads = False

    def anonymize(self, variant_record, tumor_reads_pileup, normal_reads_pileup, ref_genome) -> \
            Generator[Tuple[AnonymizedRead, AnonymizedRead], None, None]:
        called_genomic_variants = {}
        start0 = timer()
        to_yield_anonymized_reads: Set[str] = set()
        pileup_pair_iter = get_pileup_pair_in_order(tumor_reads_pileup, normal_reads_pileup)
        seen_read_alns = set()
        for pileup_pair in pileup_pair_iter:
            # for dataset_idx, current_pileup in enumerate((tumor_reads_pileup, normal_reads_pileup)):
            for dataset_idx, pileup_column in enumerate(pileup_pair):
                if pileup_column is None:
                    continue
                # is_in_normal = dataset_idx == DATASET_IDX_NORMAL
                #     for pileup_column in current_pileup:
                is_in_normal = dataset_idx == DATASET_IDX_NORMAL
                start1 = timer()
                # TODO: Multithread this part
                classify_variation_in_pileup_column(pileup_column, dataset_idx, seen_read_alns, ref_genome,
                                                    called_genomic_variants)
                end1 = timer()
                logging.debug(f"Classify variation time: {end1 - start1}")
                DEBUG_TOTAL_TIMES['classify_variants'] += end1 - start1
                # TODO: mask indels
                start4 = timer()
                for pileup_read in pileup_column.pileups:
                    aln = pileup_read.alignment
                    add_anonymized_read_pair_to_collection_from_alignment(self.anonymized_reads, aln, dataset_idx)
                    if is_candidate_to_yield(aln, pileup_column.reference_pos):
                        to_yield_anonymized_reads.add(aln.query_name)
                if is_in_normal:
                    pos = pileup_column.reference_pos
                    variants_in_column: List[CalledGenomicVariant] = called_genomic_variants.get(pos)
                    if variants_in_column is not None:
                        start2 = timer()
                        self.mask_germline_snvs(variants_in_column, variant_record)
                        end2 = timer()
                        logging.debug(f"Mask germline snvs time: {end2 - start2}")
                        DEBUG_TOTAL_TIMES['mask_germline_snvs'] += end2 - start2
                    yielded_reads = set()
                    for read_id in to_yield_anonymized_reads:
                        candidate_pair = self.anonymized_reads.get(read_id)
                        if anonymized_read_pair_is_writeable(candidate_pair[PAIR_1_IDX], candidate_pair[PAIR_2_IDX]):
                            mask_left_over_variants_in_pair(candidate_pair[PAIR_1_IDX], candidate_pair[PAIR_2_IDX])
                            yield candidate_pair
                            self.anonymized_reads.pop(read_id)
                            yielded_reads.add(read_id)
                    to_yield_anonymized_reads -= yielded_reads
                end4 = timer()
                DEBUG_TOTAL_TIMES['mask_germlines'] += end4 - start4
        end0 = timer()
        DEBUG_TOTAL_TIMES['anonymize_with_pileup'] += end0 - start0
        # Mask leftovers, referring to reads that were updated from their primary alignment,
        # after initializing with their supplementary
        start3 = timer()
        end3 = timer()
        logging.debug(f"Mask left over time: {end3 - start3}")
        DEBUG_TOTAL_TIMES['mask_germlines_left_overs_in_window'] += end3 - start3
        for read_id, anonymized_read_pair in self.anonymized_reads.items():
            mask_left_over_variants_in_pair(anonymized_read_pair[PAIR_1_IDX], anonymized_read_pair[PAIR_2_IDX])
            yield anonymized_read_pair
        self.reset()
        # self.has_anonymized_reads = True

    def mask_germline_snvs(self, variants_in_column, variant_record):
        variant_to_keep = CalledGenomicVariant.from_variant_record(variant_record)
        for called_variant in variants_in_column:
            if (called_variant.somatic_variation_type == SomaticVariationType.TUMORAL_NORMAL_VARIANT
                    and called_variant != variant_to_keep):
                # TODO: mask indels
                if called_variant.variant_type == VariantType.SNV:
                    for specific_read_id, var_read_pos in called_variant.supporting_reads.items():
                        read_id, pair = decode_specific_read_pair_name(specific_read_id)
                        anonymized_read = self.anonymized_reads.get(read_id)[pair]
                        if anonymized_read.is_supplementary:
                            anonymized_read.add_left_over_variant(var_read_pos, called_variant.ref_allele)
                            continue
                        anonymized_read.mask_or_modify_base_pair(var_read_pos, called_variant.ref_allele)
