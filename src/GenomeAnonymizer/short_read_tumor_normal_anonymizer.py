# @author: Nicolas Gaitan
import itertools
import logging
import re
import shutil
from dataclasses import dataclass
import numpy
import numpy as np
import pileup_io
import psutil
import pysam
from concurrent.futures import ProcessPoolExecutor, as_completed
from timeit import default_timer as timer
from typing import List, Dict, Tuple
from pysam import FastaFile
from variant_extractor import VariantExtractor
from variant_extractor.variants import VariantType, VariantRecord
from src.GenomeAnonymizer.anonymizer_methods import AnonymizedRead, \
    add_anonymized_read_pair_to_collection_from_alignment, \
    add_or_update_anonymized_read_from_other, anonymized_read_pair_is_writeable, Anonymizer
from src.GenomeAnonymizer.variants import CalledGenomicVariant, compare
from src.GenomeAnonymizer.variation_classifier import DATASET_IDX_TUMORAL, DATASET_IDX_NORMAL, PAIR_1_IDX, PAIR_2_IDX, \
    DEBUG_TOTAL_TIMES

"""Module for anonymizing short read tumor-normal pair genomes, using any object that implements the Anonymizer 
protocol"""

# DEBUG
process = psutil.Process()


# Define windows as namedtuples
# Window = namedtuple('Window', 'sequence first last variant')

@dataclass
class Window:
    sequence: str
    first: int
    last: int
    variant: CalledGenomicVariant = None

    def set_last_pos(self, updated_last):
        self.last = updated_last

    def is_variant_window(self):
        return self.variant is not None

    def __str__(self):
        if self.variant is None:
            return ','.join(map(str, (self.sequence, self.first, self.last)))
        else:
            return ','.join(map(str, (self.sequence, self.first, self.last, self.variant)))


def name_output(sample):
    output_suffix = '.anonymized'
    sample_name = re.sub('.bam|.sam|.cram', output_suffix, sample)
    return sample_name


def get_ref_idxs(ref_genome: FastaFile) -> Dict[str, int]:
    ref_sequences = ref_genome.references
    n_sequences = len(ref_sequences)
    return {k: v for (k, v) in zip(ref_sequences, range(n_sequences))}


def sort_window_list(windows: list[Window], ref_sequences_dict):
    windows.sort(key=lambda x: (ref_sequences_dict.get(x.sequence), x.first, x.last))


def get_windows(variants, ref_sequences_dict, window_size=2000) -> List[Window]:
    # The VCF must be sorted by ref genome qualified sequence and numeric coordinates
    half_window = int(window_size / 2)
    windows = list()
    # seq_idxs_dict = dict()
    for variant_record in variants:
        # if variant_record.contig not in seq_idxs_dict:
        #    seq_idxs_dict[variant_record.contig] = seq_idx
        #    seq_idx += 1
        called_variant = CalledGenomicVariant.from_variant_record(variant_record)
        end = variant_record.end
        if variant_record.alt_sv_breakend:
            end_chrom = variant_record.alt_sv_breakend.contig
            if variant_record.contig != end_chrom:
                end = variant_record.alt_sv_breakend.pos
        else:
            end_chrom = variant_record.contig
        if variant_record.variant_type == VariantType.INV:
            if variant_record.pos + half_window > variant_record.end - half_window:
                window = Window(sequence=variant_record.contig, first=variant_record.pos - half_window,
                                last=variant_record.end + half_window + 1, variant=called_variant)
                windows.append(window)
                # windows.append((variant_record.contig, variant_record.pos - half_window, variant_record.end + half_window, variant_record))
            else:
                window1 = Window(sequence=variant_record.contig, first=variant_record.pos - half_window,
                                 last=variant_record.pos + half_window + 1, variant=called_variant)
                window2 = Window(sequence=variant_record.contig, first=variant_record.end - half_window,
                                 last=variant_record.end + half_window + 1, variant=called_variant)
                windows.append(window1)
                windows.append(window2)
                # windows.append((variant_record.contig, variant_record.pos - half_window, variant_record.pos + half_window, variant_record))
                # windows.append((variant_record.contig, variant_record.end - half_window, variant_record.end + half_window, variant_record))
        elif variant_record.variant_type == VariantType.TRA:
            window1 = Window(sequence=variant_record.contig, first=variant_record.pos - half_window,
                             last=variant_record.pos + half_window + 1, variant=called_variant)
            window2 = Window(sequence=end_chrom, first=end - half_window,
                             last=end + half_window + 1, variant=called_variant)
            windows.append(window1)
            windows.append(window2)
            # windows.append((variant_record.contig, variant_record.pos - half_window, variant_record.pos + half_window, variant_record))
            # windows.append((variant_record.contig, variant_record.end - half_window, variant_record.end + half_window, variant_record))
        elif variant_record.variant_type == VariantType.SNV:
            window = Window(sequence=variant_record.contig, first=variant_record.pos - half_window,
                            last=variant_record.pos + half_window + 1, variant=called_variant)
            windows.append(window)
            # windows.append((variant_record.contig, variant_record.pos - half_window, variant_record.pos + half_window, variant_record))
        else:
            if variant_record.length < 100_000:
                window = Window(sequence=variant_record.contig, first=variant_record.pos - half_window,
                                last=variant_record.end + half_window + 1, variant=called_variant)
                windows.append(window)
            else:
                window1 = Window(sequence=variant_record.contig, first=variant_record.pos - half_window,
                                 last=variant_record.pos + half_window + 1, variant=called_variant)
                window2 = Window(sequence=end_chrom, first=end - half_window,
                                 last=end + half_window + 1, variant=called_variant)
                windows.append(window1)
                windows.append(window2)
    # windows.sort(key=lambda x: (ref_sequences_dict.get(x.sequence), x.first, x.last))
    sort_window_list(windows, ref_sequences_dict)
    return windows


def write_pair(indexed_writer_streams, anonymized_read_pair1, anonymized_read_pair2, written_read_ids=None):
    if written_read_ids is not None:
        # If this method is called, the AnonymizedRead pair must be writable (pairs are not None, and are complete)
        read_id = anonymized_read_pair1.query_name
        # DEBUG/
        # if read_id == 'HWI-ST1133:223:C1940ACXX:4:2314:12959:99156':
        #     logging.info(f'# assessing writing read')
        # DEBUG/
        if read_id in written_read_ids:
            # DEBUG/
            # if read_id == 'HWI-ST1133:223:C1940ACXX:4:2314:12959:99156':
            #     logging.info(f'# read in set')
            # DEBUG/
            return
        else:
            # DEBUG/
            # if read_id == 'HWI-ST1133:223:C1940ACXX:4:2314:12959:99156':
            #     logging.info(f'# read not in set')
            # DEBUG/
            written_read_ids.add(read_id)
    try:
        fastq_record_pair1 = str(anonymized_read_pair1.get_anonymized_fastq_record())
        fastq_record_pair2 = str(anonymized_read_pair2.get_anonymized_fastq_record())
        dataset_idx = anonymized_read_pair1.dataset_idx
        indexed_writer_streams[dataset_idx][PAIR_1_IDX].write(f'{fastq_record_pair1}\n')
        indexed_writer_streams[dataset_idx][PAIR_2_IDX].write(f'{fastq_record_pair2}\n')
    except TypeError as e:
        logging.error(f'Exception {str(e)} from trying to write anonymized read:\n'
                      f'Read: {anonymized_read_pair1.query_name}'
                      f'Pair 1 - SEQ = {anonymized_read_pair1.anonymized_sequence_array} - QUALs = {anonymized_read_pair1.anonymized_qualities_array}\n'
                      f'Pair 2 - SEQ = {anonymized_read_pair2.anonymized_sequence_array} - QUALs = {anonymized_read_pair2.anonymized_qualities_array}')
        raise


def close_paired_streams(indexed_pair_writer_streams):
    for pair_writer in indexed_pair_writer_streams:
        pair_writer[PAIR_1_IDX].close()
        pair_writer[PAIR_2_IDX].close()
    indexed_pair_writer_streams = []


class AnonymizedVariantsStatistics:

    outside_windows_str: str = 'outside_windows,-,-,-'

    def __init__(self, file_output: str):
        self.file_output = file_output
        # self.window_names = set()
        # self.total_counts_by_variant_type = []*len(VariantType)
        self.window_var_counts = dict()
        self.window_var_counts[self.outside_windows_str] = [0] * len(VariantType)
        self.current_window = ''
        # self.n_windows = 0

    def add_window(self, window: Window):
        # window_str = ','.join(map(str, (window.sequence, window.first, window.last, window.variant)))
        window_str = str(window)
        # window_str = ','.join(map(str, []))
        # self.window_names.add(window_str)
        self.window_var_counts[window_str] = [0] * len(VariantType)
        # Set the new window to actively compute variant statistics
        self.set_current_window(window_str)
        # self.n_windows += 1

    def count_variant(self, called_variant: CalledGenomicVariant):
        variant_type_enum = called_variant.variant_type
        var_type_idx = variant_type_enum.value - 1
        # logging.debug(f'{variant_type_enum}: {variant_type_enum.name}: {var_type_idx}')
        window_counts_by_type = self.window_var_counts.get(self.current_window)
        window_counts_by_type[var_type_idx] += 1
        # self.total_counts_by_variant_type[var_type_idx] += 1

    def set_current_window(self, window_str):
        self.current_window = window_str

    def set_outside_windows_as_current_window(self):
        self.current_window = self.outside_windows_str

    def write_statistics(self):
        var_counts_by_type = [[] for _ in range(len(VariantType))]
        stats = ['total_counts', 'average_counts', 'median_counts', 'max_counts', 'min_counts']
        with open(self.file_output, 'w') as statistics_file:
            # Variant types follow the order from their enum creation on VariantExtractor,
            # following then idxs from 0 to n types
            statistics_file.write('\t'.join(['#SEQ', '#FIRST', '#LAST', '#SNV', '#DEL', '#INS', '#DUP',
                                             '#INV', '#CNV', '#TRA', '#SGL']) + '\n')
            for window_info_key, window_counts_by_type in self.window_var_counts.items():
                window_fields = window_info_key.split(',')[:-1]
                statistics_file.write('\t'.join(map(str, itertools.chain(window_fields, window_counts_by_type))) + '\n')
                for var_type_idx, window_var_type_count in enumerate(window_counts_by_type):
                    var_counts_by_type[var_type_idx].append(window_var_type_count)
            var_types_header = '\t'.join(['#SNV', '#DEL', '#INS', '#DUP', '#INV', '#CNV', '#TRA', '#SGL']) + '\n'
            statistics_file.write('### Overall statistics:\n')
            statistics_file.write(var_types_header)
            arrays_by_type = [numpy.array(var_counts, dtype=numpy.int64) for var_counts in var_counts_by_type]
            for stat in stats:
                statistics_file.write(f'#{stat}\t')
                if stat == 'total_counts':
                    stat_by_type = [np.sum(arr) for arr in arrays_by_type]
                if stat == 'average_counts':
                    stat_by_type = [arr.mean() for arr in arrays_by_type]
                if stat == 'median_counts':
                    stat_by_type = [np.median(arr) for arr in arrays_by_type]
                if stat == 'max_counts':
                    stat_by_type = [arr.max() for arr in arrays_by_type]
                if stat == 'min_counts':
                    stat_by_type = [arr.min() for arr in arrays_by_type]
                statistics_file.write('\t'.join(map(str, stat_by_type)) + '\n')
                # f'counts: {self.total_counts_by_variant_type[variant_type.value]}')


def get_genome_sections(windows_in_sample: list[Window], ref_genome: FastaFile) -> list[Window]:
    sections: list[Window] = []
    sequences = ref_genome.references
    lengths = ref_genome.lengths
    ref_idxs = get_ref_idxs(ref_genome)
    seq_lengths = {seq: length for seq, length in zip(sequences, lengths)}
    window_dict = {k: [] for k in sequences}
    for window in windows_in_sample:
        window_dict[window.sequence].append(window)
    for seq in sequences:
        inter_window_first = 1
        seq_windows = window_dict[seq]
        if not seq_windows:
            # In this case, the whole chromosome is an inter_window
            inter_window = Window(sequence=seq, first=0, last=0)
            sections.append(inter_window)
            continue
        for window in seq_windows:
            inter_window_last = window.first - 1
            inter_window = Window(sequence=seq, first=inter_window_first, last=inter_window_last)
            inter_window_first = window.last + 1
            sections.append(inter_window)
            sections.append(window)
        inter_window_last = seq_lengths[seq] - 1
        last_inter_window = Window(sequence=seq, first=inter_window_first, last=inter_window_last)
        sections.append(last_inter_window)
    sort_window_list(sections, ref_idxs)
    # DEBUG/
    # for window in sections:
    #     logging.debug(f'{str(window)}')
    # \DEBUG
    return sections


def anonymize_window(window: Window, tumor_bam, normal_bam, tumor_output_fastq, normal_output_fastq, ref_genome,
                     anonymizer, to_pair_anonymized_reads, written_read_ids, mem_debug_writer, stats_recorder=None):
    # DEBUG/
    # logging.info(f'Anonymizing window {seq_name}:{window[0]}-{window[1]}')
    # \DEBUG
    # tumor_reads_pileup = tumor_bam.pileup(contig=seq_name, start=window[0], stop=window[1],
    #                                       fastafile=ref_genome, min_base_quality=0, max_depth=100000)
    # normal_reads_pileup = normal_bam.pileup(contig=seq_name, start=window[0], stop=window[1],
    #                                         fastafile=ref_genome, min_base_quality=0, max_depth=100000)
    start2 = timer()
    tumor_normal_pileup = pileup_io.iter_pileups(tumor_bam, normal_bam, ref_genome, seq_name=window.sequence,
                                                 start=window.first,
                                                 stop=window.last)
    anonymized_reads_generator = anonymizer.anonymize(window.variant, tumor_normal_pileup, ref_genome,
                                                      stats_recorder=stats_recorder)
    end2 = timer()
    DEBUG_TOTAL_TIMES['anonymize_call'] += end2 - start2
    # Anonymized reads generated per window
    indexed_pair_writer_streams = [
        [open(tumor_output_fastq + '.1.fastq', 'a'), open(tumor_output_fastq + '.2.fastq', 'a')],
        [open(normal_output_fastq + '.1.fastq', 'a'), open(normal_output_fastq + '.2.fastq', 'a')]]
    start3 = timer()
    # with open(tumor_output_fastq + '.1.fastq', 'a') as t1, open(tumor_output_fastq + '.2.fastq', 'a') as t2, \
    #        open(normal_output_fastq + '.1.fastq', 'a') as n1, open(normal_output_fastq + '.2.fastq', 'a') as n2:
    # indexed_pair_writer_streams = [[t1, t2], [n1, n2]]
    for anonymized_read_pair in anonymized_reads_generator:
        anonymized_read_pair1 = anonymized_read_pair[PAIR_1_IDX]
        anonymized_read_pair2 = anonymized_read_pair[PAIR_2_IDX]
        # gc.collect(2)
        # TODO: Manage anonymized pairs that have a supplementary in another window
        # TODO: Fix missing mappings with pair in another chr
        if anonymized_read_pair_is_writeable(anonymized_read_pair1, anonymized_read_pair2):
            write_pair(indexed_pair_writer_streams, anonymized_read_pair1, anonymized_read_pair2,
                       written_read_ids=written_read_ids)
        else:
            # DEBUG/
            # available_pair = anonymized_read_pair[PAIR_1_IDX] if anonymized_read_pair[PAIR_1_IDX] is not None else \
            #     anonymized_read_pair[PAIR_2_IDX]
            # if available_pair.query_name == 'HWI-ST1133:223:C1940ACXX:4:2314:12959:99156':
            #     logging.info(f'# Tracked read found in pileup')
            # \DEBUG
            if anonymized_read_pair1 is not None:
                dataset_idx = anonymized_read_pair1.dataset_idx
                add_or_update_anonymized_read_from_other(to_pair_anonymized_reads,
                                                         anonymized_read_pair1)
                read_id = anonymized_read_pair1.query_name
                # DEBUG/
                # if available_pair.query_name == 'HWI-ST1133:217:D1D4WACXX:1:2103:3953:77371':
                #     logging.info(f'# Tracked read is pair1')
                # \DEBUG
            if anonymized_read_pair2 is not None:
                dataset_idx = anonymized_read_pair2.dataset_idx
                add_or_update_anonymized_read_from_other(to_pair_anonymized_reads,
                                                         anonymized_read_pair2)
                read_id = anonymized_read_pair2.query_name
            # A read_aln pair is present in the collection, but may be the same pair, if it has supplementary alignments in other windows
            updated_anonymized_read_pair = to_pair_anonymized_reads.get(read_id)
            updated_anonymized_read_pair1 = updated_anonymized_read_pair[PAIR_1_IDX]
            updated_anonymized_read_pair2 = updated_anonymized_read_pair[PAIR_2_IDX]
            # DEBUG/
            # if available_pair.query_name == 'HWI-ST1133:223:C1940ACXX:4:2314:12959:99156':
            #     if anonymized_read_pair1 is not None:
            #         logging.info(f'# Updated Tracked read pair1 in pileup is: {updated_anonymized_read_pair1.get_status()}')
            #     if anonymized_read_pair2 is not None:
            #         logging.info(f'# Updated Tracked read pair2 in pileup is: {updated_anonymized_read_pair2.get_status()}')
            # \DEBUG
            # DEBUG/
            # if available_pair.query_name == 'HWI-ST1133:223:C1940ACXX:4:2314:12959:99156':
            #     logging.info(f'# Tracked read is being tested')
            #     if anonymized_read_pair_is_writeable(updated_anonymized_read_pair1, updated_anonymized_read_pair2):
            #         logging.info(f'# Tracked read is writable')
            #     else:
            #         logging.info(f'# Tracked read is not writable')
            # \DEBUG
            if anonymized_read_pair_is_writeable(updated_anonymized_read_pair1, updated_anonymized_read_pair2):
                if updated_anonymized_read_pair1.has_left_overs_to_mask:
                    updated_anonymized_read_pair1.mask_or_anonymize_left_over_variants()
                if updated_anonymized_read_pair2.has_left_overs_to_mask:
                    updated_anonymized_read_pair2.mask_or_anonymize_left_over_variants()
                write_pair(indexed_pair_writer_streams, updated_anonymized_read_pair1,
                           updated_anonymized_read_pair2, written_read_ids=written_read_ids)
                to_pair_anonymized_reads.pop(read_id)

    end3 = timer()
    # logging.debug(
    #    f'Time to write anonymized reads in window {seq_name} {window[0]}-{window[1]} type {window[2].variant_type.name}: {end3 - start3}')
    DEBUG_TOTAL_TIMES['mask_germlines_left_overs_in_window'] += end3 - start3
    close_paired_streams(indexed_pair_writer_streams)
    memory_usage = process.memory_info().rss / (1024 * 1024)  # in MB
    mem_debug_writer.write(
        f'Memory usage after window {window.sequence}-{window.first}-{window.last}'
        f': {memory_usage} MB\n')
    # Reset the anonymizer for the next window
    # anonymizer.reset()


def pair_unmapped_or_non_pileup_pairs_and_write(to_pair_anonymized_reads, read_aln, dataset_idx,
                                                indexed_pair_writer_streams,
                                                written_read_ids):
    add_anonymized_read_pair_to_collection_from_alignment(to_pair_anonymized_reads, read_aln,
                                                          dataset_idx)
    updated_anonymized_read_pair = to_pair_anonymized_reads.get(read_aln.query_name)
    updated_anonymized_read_pair1 = updated_anonymized_read_pair[PAIR_1_IDX]
    updated_anonymized_read_pair2 = updated_anonymized_read_pair[PAIR_2_IDX]
    # DEBUG/
    # if read_aln.query_name == 'HWI-ST1133:223:C1940ACXX:4:2314:12959:99156':
    #     logging.info(f'# Tracked read found')
    #     if read_aln.is_read1:
    #         logging.info(f'# Tracked read pair 1 found to add')
    #     if read_aln.is_read2:
    #         logging.info(f'# Tracked read pair 2 found to add')
    #     if updated_anonymized_read_pair1 is None:
    #         logging.info(f'# Tracked read pair1 is None')
    #     else:
    #         logging.info(f'# Tracked read pair1 is: {updated_anonymized_read_pair1.get_status()}')
    #     if updated_anonymized_read_pair2 is None:
    #         logging.info(f'# Tracked read pair2 is None')
    #     else:
    #         logging.info(f'# Tracked read pair2 is: {updated_anonymized_read_pair1.get_status()}')
    # \DEBUG
    if anonymized_read_pair_is_writeable(updated_anonymized_read_pair1,
                                         updated_anonymized_read_pair2):
        if updated_anonymized_read_pair1.has_left_overs_to_mask:
            updated_anonymized_read_pair1.mask_or_anonymize_left_over_variants()
        if updated_anonymized_read_pair2.has_left_overs_to_mask:
            updated_anonymized_read_pair2.mask_or_anonymize_left_over_variants()
        write_pair(indexed_pair_writer_streams, updated_anonymized_read_pair1,
                   updated_anonymized_read_pair2, written_read_ids)


def pair_unpaired_or_supplementaries(to_pair_anonymized_reads, tumor_bam_file, normal_bam_file,
                                     tumor_output_fastq, normal_output_fastq, ref_genome_file,
                                     anonymizer, written_read_ids, mem_debug_writer, threads_per_file):
    with pysam.AlignmentFile(tumor_bam_file, reference_filename=ref_genome_file) as tumor_reads_file_fetch, \
            pysam.AlignmentFile(normal_bam_file, reference_filename=ref_genome_file) as normal_reads_file_fetch, \
            pysam.AlignmentFile(tumor_bam_file, reference_filename=ref_genome_file) as tumor_reads_file_pileup, \
            pysam.AlignmentFile(normal_bam_file, reference_filename=ref_genome_file) as normal_reads_file_pileup:
        # tumor_reads_stream = tumor_reads_file.fetch(until_eof=True)
        # normal_reads_stream = normal_reads_file.fetch(until_eof=True)
        tumor_normal_fetcher = pileup_io.iter_fetch_pair(tumor_reads_file_fetch, normal_reads_file_fetch,
                                                         {*to_pair_anonymized_reads}, ref_genome_file)
        ref_genome = pysam.FastaFile(ref_genome_file)
        # ref_seq_idxs = get_ref_idxs(ref_genome)
        indexed_pair_writer_streams = [
            [open(tumor_output_fastq + '.1.fastq', 'a'), open(tumor_output_fastq + '.2.fastq', 'a')],
            [open(normal_output_fastq + '.1.fastq', 'a'), open(normal_output_fastq + '.2.fastq', 'a')]]
        # window = Window(sequence=ref_genome.references[0], first=0, last=0, )
        for fetched in tumor_normal_fetcher:
            if fetched[DATASET_IDX_NORMAL] is not None and fetched[DATASET_IDX_TUMORAL] is not None:
                seq, left, right = fetched[2]
                # logging.debug(f'seq {seq} of type {type(seq)}\nleft {left} of type {type(left)} \nright {right} of type {type(right)}')
                window = Window(sequence=seq, first=left, last=right, )
                # DEBUG/
                # if window.first == 66023812 and window.last == 66024052:
                #     logging.info(f'# Window {str(window)} is processed for tracked read')
                # \DEBUG
                anonymize_window(window, tumor_reads_file_pileup, normal_reads_file_pileup, tumor_output_fastq,
                                 normal_output_fastq,
                                 ref_genome, anonymizer, to_pair_anonymized_reads, written_read_ids,
                                 mem_debug_writer)
            elif fetched[DATASET_IDX_NORMAL] is None and fetched[DATASET_IDX_TUMORAL] is None:
                for dataset_idx in (DATASET_IDX_TUMORAL, DATASET_IDX_NORMAL):
                    unmapped_reads_tuple = fetched[2]
                    for read_aln in unmapped_reads_tuple[dataset_idx]:
                        # DEBUG/
                        # if read_aln.query_name == 'HWI-ST1133:217:D1D4WACXX:1:2103:3953:77371' and read_aln.is_read2:
                        #     logging.info(f'# Tracked read pair 2 found in unmapped')
                        # \DEBUG
                        pair_unmapped_or_non_pileup_pairs_and_write(to_pair_anonymized_reads, read_aln, dataset_idx,
                                                                    indexed_pair_writer_streams, written_read_ids)
            else:
                if fetched[DATASET_IDX_TUMORAL] is not None:
                    dataset_idx = DATASET_IDX_TUMORAL
                elif fetched[DATASET_IDX_NORMAL] is not None:
                    dataset_idx = DATASET_IDX_NORMAL
                for read_aln in fetched[dataset_idx]:
                    # DEBUG/
                    # if read_aln.query_name == 'HWI-ST1133:217:D1D4WACXX:1:2103:3953:77371':
                    #     logging.info(f'# Tracked read pair1 {read_aln.is_read1} or pair2 {read_aln.is_read2} found in non-pileup')
                    # \DEBUG
                    pair_unmapped_or_non_pileup_pairs_and_write(to_pair_anonymized_reads, read_aln, dataset_idx,
                                                                indexed_pair_writer_streams, written_read_ids)
        """for dataset_idx, current_reads_stream in enumerate((tumor_reads_stream, normal_reads_stream)):
            for read_aln in current_reads_stream:
                read_id = read_aln.query_name
                if read_id in to_pair_anonymized_reads:
                    if read_aln.is_unmapped:
                        pair_unmapped_or_non_pileup_pairs_and_write(to_pair_anonymized_reads, read_aln, dataset_idx, indexed_pair_writer_streams, written_read_ids)
                        continue
                    intersects, righmost_bp = compare_read_aln_window(read_aln, window, ref_seq_idxs)
                    if intersects:
                        window.set_last_pos(righmost_bp)
                    else:
                        anonymize_window(window, tumor_reads_file_pileup, normal_reads_file_pileup, tumor_output_fastq,
                                         normal_output_fastq,
                                         ref_genome, anonymizer, to_pair_anonymized_reads, written_read_ids,
                                         mem_debug_writer)
                        window = Window(sequence=read_aln.reference_name, first=read_aln.reference_start,
                                        last=read_aln.reference_end, )
                    # add_anonymized_read_pair_to_collection_from_alignment(to_pair_anonymized_reads, read_aln,
                    #                                                       dataset_idx)
                    # updated_anonymized_read_pair = to_pair_anonymized_reads.get(read_id)
                    # updated_anonymized_read_pair1 = updated_anonymized_read_pair[PAIR_1_IDX]
                    # updated_anonymized_read_pair2 = updated_anonymized_read_pair[PAIR_2_IDX]
                    # if anonymized_read_pair_is_writeable(updated_anonymized_read_pair1, updated_anonymized_read_pair2):
                    #     if updated_anonymized_read_pair1.has_left_overs_to_mask:
                    #         updated_anonymized_read_pair1.mask_or_anonymize_left_over_variants()
                    #     if updated_anonymized_read_pair2.has_left_overs_to_mask:
                    #         updated_anonymized_read_pair2.mask_or_anonymize_left_over_variants()
                    #     start4 = timer()
                    #     write_pair(indexed_pair_writer_streams, updated_anonymized_read_pair1,
                    #                updated_anonymized_read_pair2, written_read_ids)
                    #     end4 = timer()
                    #     DEBUG_TOTAL_TIMES['write_pairs'] += end4 - start4
                    #     to_pair_anonymized_reads.pop(read_id)"""
        close_paired_streams(indexed_pair_writer_streams)
        ref_genome.close()


def anonymize_inter_window_region(window: Window, to_pair_anonymized_reads, tumor_bam_pileup, normal_bam_pileup,
                                  tumor_bam_fetch, normal_bam_fetch, tumor_output_fastq, normal_output_fastq,
                                  ref_genome_file, anonymizer, written_read_ids, mem_debug_writer, stats_recorder=None):
    # tumor_reads_stream = tumor_reads_file.fetch(until_eof=True)
    # normal_reads_stream = normal_reads_file.fetch(until_eof=True)
    sequence = window.sequence
    first = window.first
    last = window.last
    if first + last == 0:
        first = None
        last = None
    # tumor_normal_fetcher = pileup_io.iter_fetch_pair(tumor_bam_fetch, normal_bam_fetch, {*to_pair_anonymized_reads},
    #                                                  ref_genome_file,
    #                                                  seq=sequence, first=first, last=last)
    tumor_normal_fetcher = pileup_io.iter_fetch_pair(tumor_bam_fetch, normal_bam_fetch, ref_genome_file,
                                                     seq=sequence, first=first, last=last)
    ref_genome = pysam.FastaFile(ref_genome_file)
    # ref_seq_idxs = get_ref_idxs(ref_genome)
    indexed_pair_writer_streams = [
        [open(tumor_output_fastq + '.1.fastq', 'a'), open(tumor_output_fastq + '.2.fastq', 'a')],
        [open(normal_output_fastq + '.1.fastq', 'a'), open(normal_output_fastq + '.2.fastq', 'a')]]
    # window = Window(sequence=ref_genome.references[0], first=0, last=0, )
    for fetched in tumor_normal_fetcher:
        if fetched is None:
            break
        if fetched[DATASET_IDX_NORMAL] is not None and fetched[DATASET_IDX_TUMORAL] is not None:
            seq, left, right = fetched[2]
            # logging.debug(f'seq {seq} of type {type(seq)}\nleft {left} of type {type(left)} \nright {right} of type {type(right)}')
            window = Window(sequence=seq, first=left, last=right, )
            # DEBUG/
            # if window.first == 66023812 and window.last == 66024052:
            #     logging.info(f'# Window {str(window)} is processed for tracked read')
            # \DEBUG
            anonymize_window(window, tumor_bam_pileup, normal_bam_pileup, tumor_output_fastq,
                             normal_output_fastq,
                             ref_genome, anonymizer, to_pair_anonymized_reads, written_read_ids,
                             mem_debug_writer, stats_recorder=stats_recorder)
        elif fetched[DATASET_IDX_NORMAL] is None and fetched[DATASET_IDX_TUMORAL] is None:
            for dataset_idx in (DATASET_IDX_TUMORAL, DATASET_IDX_NORMAL):
                unmapped_reads_tuple = fetched[2]
                for read_aln in unmapped_reads_tuple[dataset_idx]:
                    # DEBUG/
                    # if read_aln.query_name == 'HWI-ST1133:223:C1940ACXX:4:2314:12959:99156' and read_aln.is_read2:
                    #     logging.info(f'# Tracked read pair 2 found in unmapped')
                    # \DEBUG
                    pair_unmapped_or_non_pileup_pairs_and_write(to_pair_anonymized_reads, read_aln, dataset_idx,
                                                                indexed_pair_writer_streams, written_read_ids)
        else:
            if fetched[DATASET_IDX_TUMORAL] is not None:
                dataset_idx = DATASET_IDX_TUMORAL
            elif fetched[DATASET_IDX_NORMAL] is not None:
                dataset_idx = DATASET_IDX_NORMAL
            for read_aln in fetched[dataset_idx]:
                # DEBUG/
                # if read_aln.query_name == 'HWI-ST1133:223:C1940ACXX:4:2314:12959:99156':
                #     logging.info(
                #         f'# Tracked read pair1 {read_aln.is_read1} or pair2 {read_aln.is_read2} found in non-pileup')
                # \DEBUG
                pair_unmapped_or_non_pileup_pairs_and_write(to_pair_anonymized_reads, read_aln, dataset_idx,
                                                            indexed_pair_writer_streams, written_read_ids)
    close_paired_streams(indexed_pair_writer_streams)


def pair_unmapped_mates(windows: list[Window], to_pair_anonymized_reads, tumor_bam_file, normal_bam_file,
                        tumor_output_fastq,
                        normal_output_fastq, ref_genome_file, written_read_ids):
    indexed_pair_writer_streams = [
        [open(tumor_output_fastq + '.1.fastq', 'a'), open(tumor_output_fastq + '.2.fastq', 'a')],
        [open(normal_output_fastq + '.1.fastq', 'a'), open(normal_output_fastq + '.2.fastq', 'a')]]
    with pysam.AlignmentFile(tumor_bam_file, reference_filename=ref_genome_file) as tumor_reads_file_fetch, \
            pysam.AlignmentFile(normal_bam_file, reference_filename=ref_genome_file) as normal_reads_file_fetch:
        for window in windows:
            tumor_fetcher = tumor_reads_file_fetch.fetch(reference=window.sequence, start=window.first - 1,
                                                         stop=window.last)
            normal_fetcher = normal_reads_file_fetch.fetch(reference=window.sequence, start=window.first - 1,
                                                           stop=window.last)
            for dataset_idx, fetcher in enumerate((tumor_fetcher, normal_fetcher)):
                for read_aln in fetcher:
                    # DEBUG/
                    # if read_aln.query_name == 'HWI-ST1133:223:C1940ACXX:4:2314:12959:99156':
                    #     logging.debug(f'Tracked read found to pair: is unmapped={read_aln.is_unmapped}'
                    #                   f' is in dict={read_aln.query_name in to_pair_anonymized_reads} '
                    #                   f' saved pair=\n {to_pair_anonymized_reads[read_aln.query_name][PAIR_1_IDX]}')
                    # \DEBUG
                    if read_aln.is_unmapped and read_aln.query_name in to_pair_anonymized_reads:
                        # logging.debug(f'read_aln unmapped to pair: {read_aln}')
                        pair_unmapped_or_non_pileup_pairs_and_write(to_pair_anonymized_reads, read_aln, dataset_idx,
                                                                    indexed_pair_writer_streams, written_read_ids)
        """tumor_fetcher = tumor_reads_file_fetch.fetch(until_eof=True)
        normal_fetcher = normal_reads_file_fetch.fetch(until_eof=True)
        for dataset_idx, fetcher in enumerate((tumor_fetcher, normal_fetcher)):
            for read_aln in fetcher:
                # if read_aln.query_name in to_pair_anonymized_reads:  #
                if read_aln.query_name == 'HWI-ST1133:217:D1D4WACXX:2:2311:14232:47191':
                    logging.debug(f'Tracked read found to pair: is unmapped={read_aln.is_unmapped}'
                                  f' is in dict={read_aln.query_name in to_pair_anonymized_reads} '
                                  f' saved pair=\n {to_pair_anonymized_reads[read_aln.query_name][PAIR_1_IDX]}')
                if read_aln.is_unmapped and read_aln.query_name in to_pair_anonymized_reads:
                    # logging.debug(f'read_aln unmapped to pair: {read_aln}')
                    pair_unmapped_or_non_pileup_pairs_and_write(to_pair_anonymized_reads, read_aln, dataset_idx,
                                                                indexed_pair_writer_streams, written_read_ids)"""

    close_paired_streams(indexed_pair_writer_streams)


def write_single_end_reads(to_pair_anonymized_reads, tumor_output_fastq, normal_output_fastq):
    with (open(tumor_output_fastq + '.single_end.fastq', 'w') as tumor_fastq_writer_single_end,
          open(normal_output_fastq + '.single_end.fastq', 'w') as normal_fastq_writer_single_end):
        for read_id, anonymized_read_pair in to_pair_anonymized_reads.items():
            if anonymized_read_pair[PAIR_1_IDX] is not None:
                single_anonymized_read = anonymized_read_pair[PAIR_1_IDX]
                logging.warning('Single pair2 read not found for read id: ' + read_id)
            elif anonymized_read_pair[PAIR_2_IDX] is not None:
                single_anonymized_read = anonymized_read_pair[PAIR_2_IDX]
                logging.warning('Single pair1 read not found for read id: ' + read_id)
            if single_anonymized_read.is_supplementary:
                continue
            if single_anonymized_read.has_left_overs_to_mask:
                single_anonymized_read.mask_or_anonymize_left_over_variants()
            dataset_idx = single_anonymized_read.dataset_idx
            fastq_record = str(single_anonymized_read.get_anonymized_fastq_record())
            if dataset_idx == DATASET_IDX_TUMORAL:
                tumor_fastq_writer_single_end.write(f'{fastq_record}\n')
            elif dataset_idx == DATASET_IDX_NORMAL:
                normal_fastq_writer_single_end.write(f'{fastq_record}\n')


def anonymize_genome(windows_in_sample: List, tumor_bam_file: str, normal_bam_file: str,
                     ref_genome_file: str, anonymizer: Anonymizer, tumor_output_fastq: str,
                     normal_output_fastq: str, tumor_bam_to_pair: str, normal_bam_to_pair: str,
                     record_statistics: bool, available_threads):
    """
    Anonymizes genomic data using the provided VCF variants, normal and tumor BAM files, reference genome,
    classifier, and anonymizer object.
    """
    mem_debug_writer = open(f'{tumor_output_fastq.split("/")[-1]}_{normal_output_fastq.split("/")[-1]}.mem_debug', 'w')
    # vcf_variants = VariantExtractor(vcf_variant_file)
    threads_per_file = available_threads  # max(available_threads  - 1 // 2, 1)
    # remaining_thread = available_threads - threads_per_file
    # tumor_bam = pysam.AlignmentFile(tumor_bam_file, threads=threads_per_file+remaining_threads)
    # normal_bam = pysam.AlignmentFile(normal_bam_file, threads=threads_per_file)
    recorder = None
    if record_statistics:
        recorder = AnonymizedVariantsStatistics(f'{normal_bam_file}.statistics.txt')
    ref_genome = pysam.FastaFile(ref_genome_file)
    start1 = timer()
    # windows = get_windows(vcf_variants, ref_genome)
    end1 = timer()
    # vcf_variants.close()
    logging.debug(f'Time to get windows: {end1 - start1} s')
    # This collection contains anonymized reads that are unpaired, because their pairs are not present in the same
    # window
    to_pair_anonymized_reads: Dict[str, List[AnonymizedRead]] = dict()
    # fetch_limits_per_chr = {k: [-1,  sys.maxsize] for k in windows.keys()}
    open(tumor_output_fastq + '.1.fastq', 'w').close()
    open(tumor_output_fastq + '.2.fastq', 'w').close()
    open(normal_output_fastq + '.1.fastq', 'w').close()
    open(normal_output_fastq + '.2.fastq', 'w').close()
    memory_usage = process.memory_info().rss / (1024 * 1024)  # in MB
    mem_debug_writer.write(f'Memory usage before windows: {memory_usage} MB\n')
    written_read_ids = set()
    # current_sequence = ''
    genome_sections = get_genome_sections(windows_in_sample, ref_genome)
    with pysam.AlignmentFile(tumor_bam_file, reference_filename=ref_genome_file) as tumor_bam_pileup_windows, \
            pysam.AlignmentFile(normal_bam_file, reference_filename=ref_genome_file) as normal_bam_pileup_windows, \
            pysam.AlignmentFile(tumor_bam_file, reference_filename=ref_genome_file) as tumor_bam_fetch, \
            pysam.AlignmentFile(normal_bam_file, reference_filename=ref_genome_file) as normal_bam_fetch:  #, \
        # pysam.AlignmentFile(tumor_bam_file, reference_filename=ref_genome_file) as tumor_bam_pileup_non_windows, \
        # pysam.AlignmentFile(normal_bam_file, reference_filename=ref_genome_file) as normal_bam_pileup_non_windows:
        for window in genome_sections:
            # if current_sequence != window.sequence:
            #     written_read_ids.clear()
            # current_sequence = window.sequence
            # Window is a tuple (start, end, VariantRecord)
            start2 = timer()
            if window.is_variant_window():
                if record_statistics:
                    recorder.add_window(window)
                logging.debug(f'Anonymizing window: {str(window)}')
                anonymize_window(window, tumor_bam_pileup_windows, normal_bam_pileup_windows, tumor_output_fastq,
                                 normal_output_fastq,
                                 ref_genome, anonymizer, to_pair_anonymized_reads, written_read_ids, mem_debug_writer,
                                 stats_recorder=recorder)
            else:
                logging.debug(f'Anonymizing inter-window region: {str(window)}')
                if record_statistics:
                    recorder.set_outside_windows_as_current_window()
                anonymize_inter_window_region(window, to_pair_anonymized_reads,
                                              tumor_bam_pileup_windows, normal_bam_pileup_windows,
                                              # tumor_bam_pileup_non_windows, normal_bam_pileup_non_windows,
                                              tumor_bam_fetch, normal_bam_fetch,
                                              tumor_output_fastq, normal_output_fastq, ref_genome_file, anonymizer,
                                              written_read_ids, mem_debug_writer, stats_recorder=recorder)
            # else:
            #     anonymize_window(window, tumor_bam, normal_bam, tumor_output_fastq, normal_output_fastq,
            #                      ref_genome, anonymizer, to_pair_anonymized_reads, written_read_ids, mem_debug_writer)
            end2 = timer()
            # logging.debug(
            #    f'Time to anonymize reads in window {seq_name} {window[0]}-{window[1]} type {window[2].variant_type.name}: {end2 - start2}')
            DEBUG_TOTAL_TIMES['anonymize_windows'] += end2 - start2
    # Deal with any remaining anonymized reads, which have pairs or primary mappings in non-window regions
    memory_usage = process.memory_info().rss / (1024 * 1024)  # in MB
    mem_debug_writer.write(
        f'Memory usage after windows: {memory_usage} MB\n')
    ref_genome.close()
    """if to_pair_anonymized_reads:
        logging.info('Searching for remaining unpaired anonymized reads')
        start3 = timer()
        pair_unpaired_or_supplementaries(to_pair_anonymized_reads, tumor_bam_to_pair, normal_bam_to_pair,
                                         tumor_output_fastq, normal_output_fastq, ref_genome_file, anonymizer,
                                         written_read_ids, mem_debug_writer, threads_per_file)
        # written_read_ids, threads_per_file)
        end3 = timer()
        logging.debug(f'Time to pair unpaired anonymized reads: {end3 - start3}')
        DEBUG_TOTAL_TIMES['unpaired_searches'] += end3 - start3
        # If there are still objects in the collection, write them as single end reads
        memory_usage = process.memory_info().rss / (1024 * 1024)  # in MB
        mem_debug_writer.write(
            f'Memory usage after pairing leftover reads: {memory_usage} MB\n')
        # written_remaining_reads = list()
        # for read_id, _ in to_pair_anonymized_reads:
        #     if read_id in written_read_ids:
        #         written_remaining_reads.append(read_id)
        # for read_id in written_remaining_reads:
        #     to_pair_anonymized_reads.pop(read_id)
        # Remove all reads that are repeated due to nearby windows"""
    # NOTE: Suboptimal solution if pileup cannot return unmapped pairs
    if to_pair_anonymized_reads:
        logging.info('Searching for remaining unpaired unmapped pairs')
        start3 = timer()
        pair_unmapped_mates(windows_in_sample, to_pair_anonymized_reads, tumor_bam_file, normal_bam_file,
                            tumor_output_fastq, normal_output_fastq, ref_genome_file, written_read_ids)
        end3 = timer()
        logging.debug(f'Time to pair unpaired anonymized reads: {end3 - start3}')
        DEBUG_TOTAL_TIMES['unpaired_searches'] += end3 - start3
        # If there are still objects in the collection, write them as single end reads
        memory_usage = process.memory_info().rss / (1024 * 1024)  # in MB
        mem_debug_writer.write(
            f'Memory usage after pairing leftover reads: {memory_usage} MB\n')
    for k in written_read_ids:
        removed_read = to_pair_anonymized_reads.pop(k, '')
        if removed_read != '':
            # logging.debug(f'Duplicated read object event: read {k} was removed from remaining reads after '
            #               f'being written, this may happen because of nearby windows')
            pass
    # [to_pair_anonymized_reads.pop(k, '') for k in written_read_ids]
    if to_pair_anonymized_reads:
        start4 = timer()
        write_single_end_reads(to_pair_anonymized_reads, tumor_output_fastq, normal_output_fastq)
        end4 = timer()
        logging.debug(f'Time to write single end reads: {end4 - start4}')
        DEBUG_TOTAL_TIMES['write_pairs'] += end4 - start4
    memory_usage = process.memory_info().rss / (1024 * 1024)  # in MB
    mem_debug_writer.write(
        f'Final memory usage: {memory_usage} MB\n')
    mem_debug_writer.close()
    for k, v in DEBUG_TOTAL_TIMES.items():
        logging.debug(f'{k}={v} s')
    logging.info(f'Anonymization complete for samples {tumor_output_fastq} and {normal_output_fastq}')
    if record_statistics:
        logging.info(f'Writing anonymized variant statistics to: {recorder.file_output}')
        recorder.write_statistics()
        logging.info(f'Statistics wrote for sample: {recorder.file_output}')


def generate_subsamples_from_file(input_file, subsample_input_files, window_subsets_per_subsample, ref_genome_file,
                                  threads_per_file):
    with pysam.AlignmentFile(input_file, threads=threads_per_file, mode='rb',
                             reference_filename=ref_genome_file) as sample_file_reader:
        for subsample_file in subsample_input_files:
            windows_in_subsample = window_subsets_per_subsample.get(subsample_file)
            with pysam.AlignmentFile(subsample_file, threads=threads_per_file, mode='wb',
                                     reference_filename=ref_genome_file,
                                     header=sample_file_reader.header) as subsample_file_writer:
                for i, window in enumerate(windows_in_subsample):
                    window_it = sample_file_reader.fetch(contig=window.sequence, start=window.first,
                                                         end=window.last)
                    for read_aln in window_it:
                        subsample_file_writer.write(read_aln)


def divide_samples(inputs: List, ref_genome_file, cpus) -> (List, Dict[str, List[str]]):
    """
    Divide samples into smaller samples based on available cpus and sample sizes,
     and return their output names in a dictionary by divided output names
    """
    input_sample_keys = dict()
    output_sample_keys = dict()
    window_subsets_per_sample = dict()
    new_inputs = []
    sorted_inputs_by_size = []
    remaining_cpus = cpus - len(inputs)
    cpus_per_sample = [1] * len(inputs)
    total_size = 0
    for sample_windows, sample_t_n_pair, output_t_n_pair in inputs:
        # Compute samples window total sizes
        # sample_windows = list(chain.from_iterable(sample_windows_dict.values()))
        sample_bp_size = np.sum([abs(w.last - w.first) for w in sample_windows])
        total_size += sample_bp_size
        sorted_inputs_by_size.append((sample_windows, sample_t_n_pair, output_t_n_pair, sample_bp_size))
        # As this maintains insertion order, a t_n pair corresponds to the i (tumor) and i+1 (normal) positions
        input_sample_keys[sample_t_n_pair[DATASET_IDX_TUMORAL]] = []
        input_sample_keys[sample_t_n_pair[DATASET_IDX_NORMAL]] = []
        output_sample_keys[output_t_n_pair[DATASET_IDX_TUMORAL]] = []
        output_sample_keys[output_t_n_pair[DATASET_IDX_NORMAL]] = []
        # window_subsets_per_sample[sample_t_n_pair[DATASET_IDX_TUMORAL]] = []
        # window_subsets_per_sample[sample_t_n_pair[DATASET_IDX_NORMAL]] = []
    # Sort inputs by bp size of all the windows in each sample
    sorted_inputs_by_size.sort(key=lambda x: x[-1], reverse=True)
    bp_per_cpu = total_size // remaining_cpus
    for i, sample_tuple in enumerate(sorted_inputs_by_size):
        sample_windows, sample_t_n_pair, output_t_n_pair, sample_bp_size = sample_tuple
        # BPs serve as the tokens to assign cpus, if the remaining cpus cost less than the total possible cpus for
        # this sample, the remainers are assigned and further samples will only get 1 cpu
        sample_cpus = min(remaining_cpus, sample_bp_size // bp_per_cpu)
        cpus_per_sample[i] += sample_cpus
        remaining_cpus -= sample_cpus
        if sample_cpus == 0:
            # Then, sample keys from this sample will have an empty list of subsamples in the output_sample_keys dict
            new_input = (sample_windows, sample_t_n_pair, output_t_n_pair)
            new_inputs.append(new_input)
            break
        # Begin subsetting process of files
        sample_cpus = cpus_per_sample[i]
        n_windows = len(sample_windows)
        largest_window = max(sample_windows, key=lambda x: abs(x.last - x.first))
        n_windows_per_subsample = n_windows // sample_cpus
        compute_largest_window_alone = False
        largest_window_threshold = 1_000_000
        if (largest_window.last - largest_window.first) > largest_window_threshold:
            sample_windows.remove(largest_window)
            n_windows_per_subsample = len(sample_windows) // (sample_cpus - 1)
            compute_largest_window_alone = True
        left_limit = 0
        # print(f'onw = {n_windows} - nw = {len(sample_windows)} - cpus = {sample_cpus}')
        last_batch = sample_cpus - 2 if compute_largest_window_alone else sample_cpus - 1
        windows_to_add = 0
        for j in range(sample_cpus):
            if compute_largest_window_alone and j == sample_cpus - 1:
                # As the last batch of windows contains only the last one if it has a huge length
                windows_in_subsample = [largest_window]
                # bisect.insort(windows_in_subsample, largest_window,
                #              key=lambda x: (ref_sequences_dict.get(x.sequence), x.first, x.last))
            else:
                windows_to_add = min(n_windows_per_subsample, n_windows - left_limit)
                right_limit = left_limit + windows_to_add
                windows_in_subsample = sample_windows[left_limit::] if j == last_batch else sample_windows[
                                                                                            left_limit:right_limit]
            # print(f'j: {j} - last_batch: {last_batch} - ll: {left_limit} - rl: {right_limit} - wn: {windows_to_add}\n'
            #      f'real_n_windows: {len(windows_in_subsample)}')
            subsample_t_sample_name = f'{sample_t_n_pair[DATASET_IDX_TUMORAL]}.{j}_temp'
            subsample_n_sample_name = f'{sample_t_n_pair[DATASET_IDX_NORMAL]}.{j}_temp'
            subsample_t_output_prefix = f'{output_t_n_pair[DATASET_IDX_TUMORAL]}.{j}_temp'
            subsample_n_output_prefix = f'{output_t_n_pair[DATASET_IDX_NORMAL]}.{j}_temp'
            input_sample_keys[sample_t_n_pair[DATASET_IDX_TUMORAL]].append(subsample_t_sample_name)
            input_sample_keys[sample_t_n_pair[DATASET_IDX_NORMAL]].append(subsample_n_sample_name)
            output_sample_keys[output_t_n_pair[DATASET_IDX_TUMORAL]].append(subsample_t_output_prefix)
            output_sample_keys[output_t_n_pair[DATASET_IDX_NORMAL]].append(subsample_n_output_prefix)
            new_input = (windows_in_subsample, (subsample_t_sample_name, subsample_n_sample_name),
                         (subsample_t_output_prefix, subsample_n_output_prefix))
            new_inputs.append(new_input)
            window_subsets_per_sample[subsample_t_sample_name] = windows_in_subsample
            window_subsets_per_sample[subsample_n_sample_name] = windows_in_subsample
            left_limit += windows_to_add
    threads_by_sample_for_io = max(cpus // len(input_sample_keys), 1)
    with ProcessPoolExecutor(max_workers=cpus) as executor:
        tasks = []
        for input_file, subsample_input_files in input_sample_keys.items():
            # window_subset = window_subsets_per_sample.get(input_file)
            if not subsample_input_files:
                continue
            tasks.append(executor.submit(generate_subsamples_from_file, input_file, subsample_input_files,
                                         window_subsets_per_sample, ref_genome_file, threads_by_sample_for_io))
        for task in as_completed(tasks):
            task.result()
    return new_inputs, input_sample_keys, output_sample_keys


def join_fastq_output_from_subsamples(final_output_sample, subsample_outputs):
    with (open(final_output_sample + '.1.fastq', 'wb') as output_file_pair_1,
          open(final_output_sample + '.2.fastq', 'wb') as output_file_pair_2):
        for subsample in subsample_outputs:
            with (open(subsample + '.1.fastq', 'rb') as subsample_1,
                  open(subsample + '.2.fastq', 'rb') as subsample_2):
                shutil.copyfileobj(subsample_1, output_file_pair_1)
                shutil.copyfileobj(subsample_2, output_file_pair_2)
                # This may not be needed, as the produced files already contain a newline character
                # output_file_pair_1.write(b'\n')
                # output_file_pair_2.write(b'\n')


def run_short_read_tumor_normal_anonymizer(vcf_variants_per_sample: List[str],
                                           tumor_normal_samples: List[Tuple[str, str]],
                                           ref_genome_file: str, anonymizer: Anonymizer,
                                           output_filenames: List[Tuple[str, str]], record_statistics: bool,
                                           cpus: int, enhance_parallelization: bool):
    """
    Anonymizes genomic sequencing from short read tumor-normal pairs, in the windows from each VCF variant
    Args:
        :param vcf_variants_per_sample: The VCF variants around which the sequencing data will be anonymized, for each sample
        :param tumor_normal_samples: The list of tumor-normal samples containing the reads to be anonymized. Each sample is a
            tuple containing the tumor and normal bam files in that order
        :param ref_genome_file: The reference genome to which the reads were mapped
        :param anonymizer: The specified anonymizing method
        :param output_filenames: The output filenames for the anonymized reads, in the same format as the input samples
        :param record_statistics: If true, record statistics about anonymized variants per window and type
        :param cpus: The number of cpus to use for the anonymization of each tumor-normal sample
        :param enhance_parallelization: Divide samples by their total content in bp to maximize the cpu usage
    """
    # TODO: Implement multithreading for each sample pair
    tasks = []
    inputs_per_sample = []
    ref_genome = pysam.FastaFile(ref_genome_file)
    ref_idx_sequences = get_ref_idxs(ref_genome)
    ref_genome.close()
    for sample_vcf_variants, sample_pairs, sample_pair_outputs in zip(vcf_variants_per_sample, tumor_normal_samples,
                                                                      output_filenames):
        vars_extractor = VariantExtractor(sample_vcf_variants)
        windows_in_sample = get_windows(vars_extractor, ref_idx_sequences)
        # DEBUG/
        # if sample_pairs[0] == 'input_PCAWG_0/mosaic_genome_PCAWG_0_T.bam_0a3dd00df8fa7d14d1df99450759c11a.bam':
        #    for window in windows_in_sample:
        #         logging.debug(str(window))
        # \DEBUG
        inputs_per_sample.append((windows_in_sample, sample_pairs, sample_pair_outputs))
        vars_extractor.close()
    extra_processors = cpus - len(tumor_normal_samples)
    logging.debug(f'# extra_processors = {extra_processors}')
    # enhance_parallelization = extra_processors > 1
    if enhance_parallelization:
        inputs_per_sample, input_sample_keys, output_sample_keys = divide_samples(inputs_per_sample,
                                                                                  ref_genome_file, cpus)
        inv_input_sample = dict()
        for (k, v) in input_sample_keys.items():
            if v:
                inv_input_sample[v] = k
            else:
                inv_input_sample[k] = k
        # DEBUG
        for input_sample, samples in input_sample_keys.items():
            print(f'# input_sample = {input_sample}')
            print(f'# input_subsamples: {samples}')
        # DEBUG
        # print(f'# output sample keys: {output_sample_keys}')
    # Test up until here
    # exit()
    with ProcessPoolExecutor(max_workers=cpus) as executor:
        processes_by_sample = 1 if cpus <= len(inputs_per_sample) or enhance_parallelization else cpus // len(
            inputs_per_sample)
        for windows_in_sample, samples, sample_output_files in inputs_per_sample:
            # TODO: Implement joining the result fastq files per original samples
            original_samples = samples
            if enhance_parallelization:
                original_samples = (inv_input_sample[samples[DATASET_IDX_TUMORAL]],
                                    inv_input_sample[samples[DATASET_IDX_NORMAL]])
            tasks.append(executor.submit(anonymize_genome,
                                         windows_in_sample, samples[DATASET_IDX_TUMORAL], samples[DATASET_IDX_NORMAL],
                                         ref_genome_file, anonymizer,
                                         sample_output_files[DATASET_IDX_TUMORAL],
                                         sample_output_files[DATASET_IDX_NORMAL],
                                         original_samples[DATASET_IDX_TUMORAL], original_samples[DATASET_IDX_NORMAL],
                                         record_statistics, processes_by_sample))
        for task in as_completed(tasks):
            task.result()
        if enhance_parallelization:
            tasks = []
            for final_output_sample, subsample_outputs in output_sample_keys.items():
                tasks.append(executor.submit(join_fastq_output_from_subsamples, final_output_sample, subsample_outputs))
            for task in as_completed(tasks):
                task.result()
