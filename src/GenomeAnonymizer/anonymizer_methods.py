# @author: Nicolas Gaitan
import logging
import sys
from typing import Protocol, Dict, List, Tuple, Generator, Set, Union, NamedTuple, Any
import itertools as it
from numpy import strings
import numpy as np
import psutil
import pysam
from timeit import default_timer as timer
from pysam import PileupColumn
from variant_extractor.variants import VariantType
from src.GenomeAnonymizer.variants import CalledGenomicVariant, SomaticVariationType
from src.GenomeAnonymizer.variation_classifier import classify_variation_in_pileup_column, DATASET_IDX_TUMORAL, \
    DATASET_IDX_NORMAL, PAIR_1_IDX, PAIR_2_IDX, DEBUG_TOTAL_TIMES

# DEBUG/
process = psutil.Process()
# \DEBUG


reverses = {ord('A'): ord('T'), ord('C'): ord('G'), ord('G'): ord('C'), ord('T'): ord('A'), ord('N'): ord('N')}
# BASE_TO_INT = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
# INT_TO_BASE = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'N'}
complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


def encode_base(base: str):
    return bytearray(base.encode())[0]


def encode_sequence_as_np_array(sequence: str) -> np.ndarray[Any, np.dtype[np.bytes_]]:
    return np.frombuffer(bytearray(sequence.encode()), dtype=np.uint8)


# def encode_base(base: str):
#     return base.encode('ascii')


# def encode_base(base: str) -> int:
#     return BASE_TO_INT[base]


def decode_base(int_base: int) -> str:
    return chr(int_base)


# def decode_base(int_base: int) -> str:
#     return INT_TO_BASE[int_base]


# def get_numeric_complement(fw_base: int) -> int:
#     rs_compl: str = complements[decode_base(fw_base)]
#     return encode_base(rs_compl)


def generate_anonymized_read(name, sequence, quality):
    return f'@{name}\n{sequence}\n+\n{quality}'


def get_supplementary_hash_from_aln(aln: pysam.AlignedSegment):
    return f'{aln.reference_name};{aln.reference_start};{aln.cigarstring};{aln.query_sequence}'


# reverses = {encode_base(base): encode_base(compl_base) for base, compl_base in complements.items()}


"""def get_pileup_pair_in_order(pileup1, pileup2) -> Generator[Tuple[PileupColumn, PileupColumn], None, None]:
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
            p2 = next(pileup2, None)"""


class AnonymizedRead:
    def __init__(self, read_alignment: pysam.AlignedSegment, dataset_idx: int):
        self.anonymized_sequence_array: np.ndarray = []
        # self.anonymized_sequence_array: np.ndarray[Any, np.dtype[np.bytes_]] = np.array([], dtype=np.uint8)
        self.anonymized_qualities_array: list[int] = []
        self.query_name = read_alignment.query_name
        self.is_read1 = read_alignment.is_read1
        self.is_read2 = read_alignment.is_read2
        self.is_reverse = read_alignment.is_reverse
        # self.read_alignment = read_alignment
        self.set_original_sequence(read_alignment.query_sequence)
        self.set_original_qualities(read_alignment.get_forward_qualities())
        # print(read_alignment.get_forward_qualities())
        self.dataset_idx = dataset_idx
        # An AnonymizedRead is_supplementary, if the read has a supp. alignment and the primary mapping has not been processed
        self.is_supplementary = read_alignment.is_supplementary
        self.has_supplementary = read_alignment.has_tag('SA')
        self.supplementary_hashes = set()
        self.n_supplementaries = 0
        if self.has_supplementary:
            supplementaries = read_alignment.get_tag('SA').rstrip(';').split(';')
            self.n_supplementaries = len(supplementaries) - 1
            if self.is_supplementary:
                self.record_supplementary_aln(get_supplementary_hash_from_aln(read_alignment))
        # DEBUG/
        # if self.query_name == 'HWI-ST1133:217:D1D4WACXX:2:1204:12012:86477':
        #     logging.info(f'# n_supplementaries={self.n_supplementaries}'
        #                  f' supplementaries_recorded={len(self.supplementary_hashes)}')
        # \DEBUG
        # self.supplementary_recorded = self.is_supplementary and self.has_supplementary
        # Left over variants that were found in the suppl. but the primary mapping was not found yet
        self.left_over_variants_to_mask: List[Tuple[int, CalledGenomicVariant]] = []
        self.has_left_overs_to_mask = False

    def get_pair_idx(self):
        if self.is_read1:
            return PAIR_1_IDX
        if self.is_read2:
            return PAIR_2_IDX

    def anonymized_read_is_complete(self) -> bool:
        if self.is_supplementary:
            return False
        if self.has_supplementary:
            # DEBUG/
            # assert (self.supplementaries_recorded <= self.n_supplementaries,
            #         f'Recorded supplementary alignments: {self.supplementaries_recorded}'
            #         f' cannot be more than the total: {self.n_supplementaries}'
            #         f' supplementary alignments.')
            # \DEBUG
            if len(self.supplementary_hashes) < self.n_supplementaries:
                return False
        return True

    def record_supplementary_aln(self, supplementary_hash: str):
        self.supplementary_hashes.add(supplementary_hash)

    def update_from_primary_mapping(self, aln: pysam.AlignedSegment):
        if aln.is_supplementary:
            raise ValueError("Trying to update AnonymizedRead using a supplementary alignment: "
                             "The update should always be called only if the primary mapping appears")
        # self.read_alignment = aln
        self.set_original_sequence(aln.query_sequence)
        self.set_original_qualities(aln.get_forward_qualities())
        self.is_supplementary = False

    def set_original_sequence(self, original_sequence: str):
        # self.anonymized_sequence_array: np.ndarray[Any, np.dtype[np.bytes_]] = encode_sequence_as_np_array(
        #     original_sequence
        # )
        # self.anonymized_sequence_array: List[int] = [encode_base(base) for base in original_sequence]
        # self.anonymized_sequence_array: np.array = np.array([encode_base(base) for base in original_sequence],
        #                                                     dtype=np.uint8)
        # self.anonymized_sequence_array: np.ndarray[Any, np.dtype[np.bytes_]] = np.strings.encode(
        #     [bp for bp in original_sequence],
        #     encoding='ascii')
        # self.anonymized_sequence_array: np.ndarray = np.frombuffer(bytearray(original_sequence.encode()),
        #                                                            dtype=np.uint8)
        self.anonymized_sequence_array: np.ndarray = encode_sequence_as_np_array(original_sequence)

    def set_original_qualities(self, original_qualities):
        self.anonymized_qualities_array = original_qualities
        # self.anonymized_qualities_array: np.ndarray = np.frombuffer(original_qualities,
        #                                                             dtype=np.uint8)

    def mask_or_modify_base_pair(self, pos_in_read: int, new_base: str, modify_qualities: bool = False,
                                 new_quality: int = 0):
        # self.anonymized_sequence_array[pos_in_read] = encode_base(new_base)
        # self.anonymized_sequence_array[pos_in_read] = new_base
        np.put(self.anonymized_sequence_array, pos_in_read, encode_base(new_base), mode='raise')
        if modify_qualities:
            self.anonymized_qualities_array[pos_in_read] = new_quality

    def mask_or_modify_indel(self, var_pos_in_read, variant):
        # TODO: Modify to adapt to numpy array instead of regular list
        original_sequence = self.anonymized_sequence_array
        if variant.variant_type == VariantType.INS:
            # Deletes the insertion event
            new_sequence: list[int] = original_sequence[0:var_pos_in_read] + original_sequence[
                                                                             var_pos_in_read + variant.length:]
        elif variant.variant_type == VariantType.DEL:
            # Adds the deleted reference bases
            new_sequence: list[int] = original_sequence[0:var_pos_in_read] + [base for base in
                                                                              variant.ref_allele] + original_sequence[
                                                                                                    variant.end:]
        else:
            # Placeholder for SVs
            new_sequence = original_sequence
        self.anonymized_sequence_array = new_sequence

    def reverse_complement(self):
        # complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        # complement = {0: 3, 1: 2, 2: 1, 3: 0, 4: 4}
        # self.anonymized_sequence_array = np.flip(np.vectorize(
        #     lambda x: complements[x.decode()])(self.anonymized_sequence_array))
        # self.anonymized_sequence_array = [get_numeric_complement(base) for base in
        #                                   reversed(self.anonymized_sequence_array)]
        self.anonymized_sequence_array = np.flip(np.vectorize(reverses.get)(self.anonymized_sequence_array))
        self.anonymized_qualities_array = reversed(self.anonymized_qualities_array)

    def get_anonymized_fastq_record(self) -> str:
        if self.is_reverse:
            self.reverse_complement()
        # else:
        #     self.anonymized_sequence_array = strings.decode(self.anonymized_sequence_array, encoding='ascii')
        read_pair_name = f'{self.query_name}/{PAIR_1_IDX + 1}' if self.is_read1 else f'{self.query_name}/{PAIR_2_IDX + 1}'
        # anonymized_read_seq: str = ''.join(self.anonymized_sequence_array)
        # anonymized_read_seq: str = np.strings.decode(self.anonymized_sequence_array, encoding='ascii')
        anonymized_read_seq = ''.join(map(lambda x: decode_base(x), self.anonymized_sequence_array))
        # anonymized_read_seq: str = np.array2string(self.anonymized_sequence_array,
        #                                            formatter={'int': lambda x: decode_base(x)},
        #                                            separator='',
        #                                            max_line_width=sys.maxsize).strip('[]')
        # anonimized_read_qual: str = ''.join(map(lambda x: chr(x + 33), self.anonymized_qualities_array))
        # anonimized_read_qual: str = np.array2string(self.anonymized_qualities_array,
        #                                            formatter={'int': lambda x: chr(x + 33)},
        #                                            separator='',
        #                                            max_line_width=sys.maxsize).strip('[]')
        anonimized_read_qual: str = ''.join([chr(x + 33) for x in self.anonymized_qualities_array])
        # anonimized_read_qual: str = ''.join([chr(x) for x in self.anonymized_qualities_array])
        # if len(anonymized_read_seq) != len(self.original_sequence):
        #    raise ValueError("Anonymized read length does not match original read length")
        # if len(anonimized_read_qual) != len(self.original_qualities):
        #    raise ValueError("Anonymized qualities length does not match original qualities length")
        # anonymized_read: pysam.FastxRecord = pysam.FastxRecord(name=read_pair_name, sequence=anonymized_read_seq,
        #                                                        quality=anonimized_read_qual)
        anonymized_read = generate_anonymized_read(name=read_pair_name, sequence=anonymized_read_seq,
                                                   quality=anonimized_read_qual)
        return anonymized_read

    def add_left_over_variant(self, var_pos_in_read: int, variant: CalledGenomicVariant):
        if not self.is_supplementary:
            if variant.variant_type == VariantType.SNV:
                raise ValueError(
                    f'Trying to add left over SNV variant to AnonymizedRead {self.query_name} containing a primary mapping'
                    f'\n all SNVs can be masked already')
        self.left_over_variants_to_mask.append((var_pos_in_read, variant))
        self.has_left_overs_to_mask = True

    def mask_or_anonymize_left_over_variants(self):
        if self.is_supplementary:
            raise ValueError(
                f'Trying to mask left over variants in AnonymizedRead {self.query_name} without a primary mapping')
        if not self.has_left_overs_to_mask or len(self.left_over_variants_to_mask) == 0:
            logging.warning(
                f'Trying to mask left over variants in AnonymizedRead {self.query_name} '
                f'with no left over variants to mask')
        for var_pos_in_read, variant in self.left_over_variants_to_mask:
            if variant.variant_type == VariantType.SNV:
                self.mask_or_modify_base_pair(var_pos_in_read, variant.allele)
            if variant.variant_type == VariantType.DEL or variant.variant_type == VariantType.INS:
                # TODO: Implement masking to reference for indels
                # self.mask_or_modify_indel(var_pos_in_read, variant)
                pass
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
            # new_anonymized_read.supplementary_recorded = True
        if aln.is_supplementary:
            new_anonymized_read.record_supplementary_aln(get_supplementary_hash_from_aln(aln))
        # DEBUG/
        # if aln.query_name == 'HWI-ST1133:217:D1D4WACXX:2:1204:12012:86477':
        #     logging.info(f'# from aln - n_supplementaries={new_anonymized_read.n_supplementaries}'
        #                  f' supplementaries_recorded={len(new_anonymized_read.supplementary_hashes)}'
        #                  f' cigar={aln.cigarstring}')
        # \DEBUG


def add_or_update_anonymized_read_from_other(anonymized_reads: Dict[str, List[AnonymizedRead]],
                                             anonymized_read: AnonymizedRead):
    if anonymized_read.query_name not in anonymized_reads:
        anonymized_reads[anonymized_read.query_name] = [None, None]
        paired_anonymized_read_list = anonymized_reads.get(anonymized_read.query_name)
        pair_idx = anonymized_read.get_pair_idx()
        paired_anonymized_read_list[pair_idx] = anonymized_read
    else:
        paired_anonymized_read_list = anonymized_reads.get(anonymized_read.query_name)
        pair_idx = anonymized_read.get_pair_idx()
        if paired_anonymized_read_list[pair_idx] is None:
            paired_anonymized_read_list[pair_idx] = anonymized_read
            return
        saved_anonymized_read = paired_anonymized_read_list[pair_idx]
        # Saved anonymized read is supplementary and the new one is a primary mapping
        if saved_anonymized_read.is_supplementary and not anonymized_read.is_supplementary:
            anonymized_read.left_over_variants_to_mask.extend(saved_anonymized_read.left_over_variants_to_mask)
            if len(anonymized_read.left_over_variants_to_mask) > 0:
                anonymized_read.has_left_overs_to_mask = True
            paired_anonymized_read_list[pair_idx] = anonymized_read
            return
        # In any case, update leftover bases to mask from the new anonymized read
        if anonymized_read.has_left_overs_to_mask:
            saved_anonymized_read.left_over_variants_to_mask.extend(anonymized_read.left_over_variants_to_mask)
            if len(saved_anonymized_read.left_over_variants_to_mask) > 0:
                saved_anonymized_read.has_left_overs_to_mask = True
        # If the new anonymized read is supplementary update the supp. record
        if anonymized_read.is_supplementary:
            for suppl_hash in anonymized_read.supplementary_hashes:
                saved_anonymized_read.record_supplementary_aln(suppl_hash)
        # DEBUG/
        # if anonymized_read.query_name == 'HWI-ST1133:217:D1D4WACXX:2:1204:12012:86477':
        #     logging.info(f'# from anon_read - n_supplementaries={anonymized_read.n_supplementaries}'
        #                  f' supplementaries_recorded={len(anonymized_read.supplementary_hashes)}')
        # \DEBUG


def anonymized_read_pair_is_writeable(anonymized_read_pair1: AnonymizedRead,
                                      anonymized_read_pair2: AnonymizedRead) -> bool:
    # DEBUG/
    # available_pair = anonymized_read_pair1 if anonymized_read_pair1 is not None else \
    #     anonymized_read_pair2
    # if available_pair.query_name == 'HWI-ST0738:314:D1E5KACXX:5:1204:15477:168431':
    #     logging.info(f'# Tracked read assessing if writeable\n'
    #                  f'pair1 is None: {anonymized_read_pair1 is None}\n'
    #                  f'pair2 is None: {anonymized_read_pair2 is None}\n'
    #                  f'pair1 is complete: {anonymized_read_pair1.anonymized_read_is_complete()}\n'
    #                  f'pair2 is complete: {anonymized_read_pair2.anonymized_read_is_complete()}')
    # \DEBUG
    if anonymized_read_pair1 is None or anonymized_read_pair2 is None:
        return False
    # if anonymized_read_pair1.is_supplementary or anonymized_read_pair2.is_supplementary:
    #     return False
    if not anonymized_read_pair1.anonymized_read_is_complete() or not anonymized_read_pair2.anonymized_read_is_complete():
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

    def anonymize(self, validated_source_variant: CalledGenomicVariant, tumor_normal_pileup, ref_genome,
                  stats_recorder=None) -> Generator[Tuple[AnonymizedRead, AnonymizedRead], None, None]:
        called_genomic_variants = {}
        start0 = timer()
        to_yield_anonymized_reads: Dict[str, int] = dict()
        # pileup_pair_iter = get_pileup_pair_in_order(tumor_reads_pileup, normal_reads_pileup)
        seen_read_alns = set()
        for pileup_pair in tumor_normal_pileup:
            # for dataset_idx, current_pileup in enumerate((tumor_reads_pileup, normal_reads_pileup)):
            # print(f'{pileup_pair} - len: {len(pileup_pair)}')
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
                # logging.debug(f"Classify variation time: {end1 - start1}")
                DEBUG_TOTAL_TIMES['classify_variants'] += end1 - start1
                # TODO: mask indels
                start4 = timer()
                for pileup_read in pileup_column.pileups:
                    aln = pileup_read.alignment
                    add_anonymized_read_pair_to_collection_from_alignment(self.anonymized_reads, aln, dataset_idx)
                    if aln.query_name not in to_yield_anonymized_reads:
                        to_yield_anonymized_reads[aln.query_name] = aln.reference_end
                    else:
                        to_yield_anonymized_reads[aln.query_name] = max(to_yield_anonymized_reads[aln.query_name],
                                                                        aln.reference_end)
                if is_in_normal:
                    pos = pileup_column.reference_pos
                    variants_in_column: List[CalledGenomicVariant] = called_genomic_variants.get(pos)
                    if variants_in_column is not None:
                        start2 = timer()
                        self.mask_germline_snvs(variants_in_column, validated_source_variant,
                                                stats_recorder=stats_recorder)
                        end2 = timer()
                        # logging.debug(f"Mask germline snvs time: {end2 - start2}")
                        DEBUG_TOTAL_TIMES['mask_germline_snvs'] += end2 - start2
                    new_yielded_reads = set()
                    for read_id, right_most_end in to_yield_anonymized_reads.items():
                        candidate_pair = self.anonymized_reads.get(read_id)
                        is_candidate_to_yield = right_most_end < pos
                        if is_candidate_to_yield and anonymized_read_pair_is_writeable(candidate_pair[PAIR_1_IDX],
                                                                                       candidate_pair[PAIR_2_IDX]):
                            mask_left_over_variants_in_pair(candidate_pair[PAIR_1_IDX], candidate_pair[PAIR_2_IDX])
                            # DEBUG/
                            # available_pair = candidate_pair[PAIR_1_IDX] if candidate_pair[
                            #                                                          PAIR_1_IDX] is not None else \
                            #     candidate_pair[PAIR_2_IDX]
                            # if available_pair.query_name == 'HWI-ST0738:314:D1E5KACXX:5:1204:15477:168431':
                            #     logging.info(f'# Tracked read is about to be yielded in stream')
                            # \DEBUG
                            yield candidate_pair
                            self.anonymized_reads.pop(read_id)
                            new_yielded_reads.add(read_id)
                        # DEBUG/
                        # else:
                        #     logging.debug(f'# mem_usage={process.memory_info().rss / (1024 * 1024)} MB\t'
                        #                   f'n_left_overs={len(self.anonymized_reads)}')
                        # \DEBUG
                    for read_id in new_yielded_reads:
                        to_yield_anonymized_reads.pop(read_id)
                    # to_yield_anonymized_reads -= new_yielded_reads
                end4 = timer()
                DEBUG_TOTAL_TIMES['mask_germlines'] += end4 - start4
        end0 = timer()
        DEBUG_TOTAL_TIMES['anonymize_with_pileup'] += end0 - start0
        # Mask leftovers, referring to reads that were updated from their primary alignment,
        # after initializing with their supplementary
        start3 = timer()
        for read_id, anonymized_read_pair in self.anonymized_reads.items():
            start3 = timer()
            mask_left_over_variants_in_pair(anonymized_read_pair[PAIR_1_IDX], anonymized_read_pair[PAIR_2_IDX])
            end3 = timer()
            DEBUG_TOTAL_TIMES['mask_germlines_left_overs_in_window'] += end3 - start3
            # DEBUG/
            # available_pair = anonymized_read_pair[PAIR_1_IDX] if anonymized_read_pair[PAIR_1_IDX] is not None else \
            #     anonymized_read_pair[PAIR_2_IDX]
            # if available_pair.query_name == 'HWI-ST0738:314:D1E5KACXX:5:1204:15477:168431':
            #     logging.info(f'# Tracked read is about to be yielded from leftovers')
            # \DEBUG
            yield anonymized_read_pair
        # logging.debug(f"Mask left over time: {end3 - start3}")
        self.reset()
        # self.has_anonymized_reads = True

    def mask_germline_snvs(self, variants_in_column: List[CalledGenomicVariant], variant_to_keep, stats_recorder=None):
        """Mask all germline SNVs and save the information from incomplete pairs and other variants to mask after"""
        for called_variant in variants_in_column:
            if (called_variant.somatic_variation_type == SomaticVariationType.TUMORAL_NORMAL_VARIANT
                    and called_variant != variant_to_keep):
                for specific_read_id, var_read_pos in called_variant.supporting_reads.items():
                    read_id, pair = decode_specific_read_pair_name(specific_read_id)
                    anonymized_read = self.anonymized_reads.get(read_id)[pair]
                    if anonymized_read.is_supplementary or called_variant.variant_type != VariantType.SNV:
                        anonymized_read.add_left_over_variant(var_read_pos, called_variant)
                        continue
                    anonymized_read.mask_or_modify_base_pair(var_read_pos, called_variant.ref_allele)
                # 11/06/24: For now, only SNVs are anonymized, but indels are also counted
                stats_recorder.count_variant(called_variant)
