# @author: Nicolas Gaitan
import gc
import gzip
import logging
import os
import sys
from timeit import default_timer as timer
import psutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Tuple, Iterable, Dict
import pysam
from variant_extractor import VariantExtractor
from variant_extractor.variants import VariantType
from src.GenomeAnonymizer.anonymizer_methods import Anonymizer, AnonymizedRead, \
    add_anonymized_read_pair_to_collection_from_alignment, \
    add_or_update_anonymized_read_from_other, anonymized_read_pair_is_writeable
from src.GenomeAnonymizer.variation_classifier import DATASET_IDX_TUMORAL, DATASET_IDX_NORMAL, PAIR_1_IDX, PAIR_2_IDX, \
    DEBUG_TOTAL_TIMES

"""Module for anonymizing short read tumor-normal pair genomes, using any object that implements the Anonymizer 
protocol"""

# DEBUG
process = psutil.Process()


def get_windows(variants, ref_genome, window_size=2000):
    half_window = int(window_size / 2)
    windows = dict()
    for seq in ref_genome.references:
        windows[seq] = []
    for variant in variants:
        windows_in_seq = windows[variant.contig]
        if variant.variant_type == VariantType.INV:
            if variant.pos + half_window > variant.end - half_window:
                windows_in_seq.append((variant.pos - half_window, variant.end + half_window, variant))
            else:
                windows_in_seq.append((variant.pos - half_window, variant.pos + half_window, variant))
                windows_in_seq.append((variant.end - half_window, variant.end + half_window, variant))
        elif variant.variant_type == VariantType.TRA:
            windows_in_seq.append((variant.pos - half_window, variant.pos + half_window, variant))
            windows_in_seq.append((variant.end - half_window, variant.end + half_window, variant))
        elif variant.variant_type == VariantType.SNV:
            windows_in_seq.append((variant.pos - half_window, variant.pos + half_window, variant))
        else:
            windows_in_seq.append((variant.pos - half_window, variant.end + half_window, variant))
    for seq_name, seq_windows in windows.items():
        seq_windows.sort(key=lambda x: (x[0], x[1]))
    return windows


def write_pair(indexed_writer_streams, anonymized_read_pair1, anonymized_read_pair2):
    fastq_record_pair1 = str(anonymized_read_pair1.get_anonymized_fastq_record())
    fastq_record_pair2 = str(anonymized_read_pair2.get_anonymized_fastq_record())
    dataset_idx = anonymized_read_pair1.dataset_idx
    indexed_writer_streams[dataset_idx][PAIR_1_IDX].write(f'{fastq_record_pair1}\n')
    indexed_writer_streams[dataset_idx][PAIR_2_IDX].write(f'{fastq_record_pair2}\n')


def close_paired_streams(indexed_pair_writer_streams):
    for pair_writer in indexed_pair_writer_streams:
        pair_writer[PAIR_1_IDX].close()
        pair_writer[PAIR_2_IDX].close()
    indexed_pair_writer_streams = []


# def anonymize_window(seq_name, window, tumor_bam, normal_bam, tumor_output_fastq, normal_output_fastq, ref_genome,
#                      anonymizer, to_pair_anonymized_reads)
def anonymize_window(seq_name, window, tumor_bam, normal_bam, tumor_output_fastq, normal_output_fastq, ref_genome,
                     anonymizer, to_pair_anonymized_reads, mem_debug_writer):
    tumor_reads_pileup = tumor_bam.pileup(contig=seq_name, start=window[0], stop=window[1],
                                          fastafile=ref_genome, min_base_quality=0)
    normal_reads_pileup = normal_bam.pileup(contig=seq_name, start=window[0], stop=window[1],
                                            fastafile=ref_genome, min_base_quality=0)
    start2 = timer()
    anonymized_reads_generator = anonymizer.anonymize(window[2], tumor_reads_pileup, normal_reads_pileup, ref_genome)
    end2 = timer()
    DEBUG_TOTAL_TIMES['anonymize_call'] += end2 - start2
    # Anonymized reads generated per window
    # TODO: Actually yield the anonymized reads without saving to memory for ever
    indexed_pair_writer_streams = [
        [open(tumor_output_fastq + '.1.fastq', 'a'), open(tumor_output_fastq + '.2.fastq', 'a')],
        [open(normal_output_fastq + '.1.fastq', 'a'), open(normal_output_fastq + '.2.fastq', 'a')]]
    start3 = timer()
    for anonymized_read_pair in anonymized_reads_generator:
        anonymized_read_pair1 = anonymized_read_pair[PAIR_1_IDX]
        anonymized_read_pair2 = anonymized_read_pair[PAIR_2_IDX]
        # gc.collect(2)
        if anonymized_read_pair_is_writeable(anonymized_read_pair1, anonymized_read_pair2):
            # TODO: Manage anonymized pairs that have a supplementary in another window
            write_pair(indexed_pair_writer_streams, anonymized_read_pair1, anonymized_read_pair2)
        else:
            if anonymized_read_pair1 is not None:
                dataset_idx = anonymized_read_pair1.dataset_idx
                add_or_update_anonymized_read_from_other(to_pair_anonymized_reads,
                                                         anonymized_read_pair1)
                read_id = anonymized_read_pair1.query_name
            if anonymized_read_pair2 is not None:
                dataset_idx = anonymized_read_pair2.dataset_idx
                add_or_update_anonymized_read_from_other(to_pair_anonymized_reads,
                                                         anonymized_read_pair2)
                read_id = anonymized_read_pair2.query_name
            # A read_aln pair is present in the collection, but may be the same pair, if it has supplementary alignments in other windows
            updated_anonymized_read_pair = to_pair_anonymized_reads.get(read_id)
            updated_anonymized_read_pair1 = updated_anonymized_read_pair[PAIR_1_IDX]
            updated_anonymized_read_pair2 = updated_anonymized_read_pair[PAIR_2_IDX]
            if anonymized_read_pair_is_writeable(updated_anonymized_read_pair1, updated_anonymized_read_pair2):
                if updated_anonymized_read_pair1.has_left_overs_to_mask:
                    updated_anonymized_read_pair1.mask_or_anonymize_left_over_variants()
                if updated_anonymized_read_pair2.has_left_overs_to_mask:
                    updated_anonymized_read_pair2.mask_or_anonymize_left_over_variants()
                write_pair(indexed_pair_writer_streams, updated_anonymized_read_pair1,
                           updated_anonymized_read_pair2)
                to_pair_anonymized_reads.pop(read_id)
                # anonymized_read_pair = updated_anonymized_read_pair[pair_not_missing]
                # pair_chr
                # limits_in_chr = fetch_limits_per_chr[seq_name]
    end3 = timer()
    logging.debug(
        f'Time to write anonymized reads in window {seq_name} {window[0]}-{window[1]} type {window[2].variant_type.name}: {end3 - start3}')
    DEBUG_TOTAL_TIMES['mask_germlines_left_overs_in_window'] += end3 - start3
    close_paired_streams(indexed_pair_writer_streams)
    memory_usage = process.memory_info().rss / (1024 * 1024)  # in MB
    mem_debug_writer.write(
        f'Memory usage after window {seq_name}-{window[0]}-{window[1]}'
        f': {memory_usage} MB\n')
    # Reset the anonymizer for the next window
    # anonymizer.reset()


def pair_unpaired_or_supplementaries(to_pair_anonymized_reads, tumor_bam_file, normal_bam_file,
                                     tumor_output_fastq, normal_output_fastq, threads_per_file):
    with pysam.AlignmentFile(tumor_bam_file) as tumor_reads_file, \
            pysam.AlignmentFile(normal_bam_file) as normal_reads_file:
        tumor_reads_stream = tumor_reads_file.fetch(until_eof=True)
        normal_reads_stream = normal_reads_file.fetch(until_eof=True)
        indexed_pair_writer_streams = [
            [open(tumor_output_fastq + '.1.fastq', 'a'), open(tumor_output_fastq + '.2.fastq', 'a')],
            [open(normal_output_fastq + '.1.fastq', 'a'), open(normal_output_fastq + '.2.fastq', 'a')]]
        for dataset_idx, current_reads_stream in enumerate((tumor_reads_stream, normal_reads_stream)):
            for read_aln in current_reads_stream:
                read_id = read_aln.query_name
                if read_id in to_pair_anonymized_reads:
                    add_anonymized_read_pair_to_collection_from_alignment(to_pair_anonymized_reads, read_aln,
                                                                          dataset_idx)
                    updated_anonymized_read_pair = to_pair_anonymized_reads.get(read_id)
                    updated_anonymized_read_pair1 = updated_anonymized_read_pair[PAIR_1_IDX]
                    updated_anonymized_read_pair2 = updated_anonymized_read_pair[PAIR_2_IDX]
                    if anonymized_read_pair_is_writeable(updated_anonymized_read_pair1, updated_anonymized_read_pair2):
                        if updated_anonymized_read_pair1.has_left_overs_to_mask:
                            updated_anonymized_read_pair1.mask_or_anonymize_left_over_variants()
                        if updated_anonymized_read_pair2.has_left_overs_to_mask:
                            updated_anonymized_read_pair2.mask_or_anonymize_left_over_variants()
                        start4 = timer()
                        write_pair(indexed_pair_writer_streams, updated_anonymized_read_pair1,
                                   updated_anonymized_read_pair2)
                        end4 = timer()
                        DEBUG_TOTAL_TIMES['write_pairs'] += end4 - start4
                        to_pair_anonymized_reads.pop(read_id)
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


def anonymize_genome(vcf_variant_file: str, tumor_bam_file: str, normal_bam_file: str,
                     ref_genome_file: str, anonymizer: Anonymizer, tumor_output_fastq: str,
                     normal_output_fastq: str, available_threads: int):
    """
    Anonymizes genomic data using the provided VCF variants, normal and tumor BAM files, reference genome,
    classifier, and anonymizer object.
    """
    mem_debug_writer = open(f'{tumor_output_fastq.split("/")[-1]}_{normal_output_fastq.split("/")[-1]}.mem_debug', 'w')
    vcf_variants = VariantExtractor(vcf_variant_file)
    threads_per_file = available_threads  # max(available_threads  - 1 // 2, 1)
    # remaining_thread = available_threads - threads_per_file
    # tumor_bam = pysam.AlignmentFile(tumor_bam_file, threads=threads_per_file+remaining_threads)
    # normal_bam = pysam.AlignmentFile(normal_bam_file, threads=threads_per_file)
    ref_genome = pysam.FastaFile(ref_genome_file)
    start1 = timer()
    windows = get_windows(vcf_variants, ref_genome)
    end1 = timer()
    vcf_variants.close()
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
    for seq_name, seq_windows in windows.items():
        for window in seq_windows:
            with pysam.AlignmentFile(tumor_bam_file) as tumor_bam, \
                    pysam.AlignmentFile(normal_bam_file) as normal_bam:
                # Window is a tuple (start, end, VariantRecord)
                start2 = timer()
                # anonymize_window(seq_name, window, tumor_bam, normal_bam, tumor_output_fastq, normal_output_fastq,
                #                 ref_genome, anonymizer, to_pair_anonymized_reads)
                #DEBUG CALL
                anonymize_window(seq_name, window, tumor_bam, normal_bam, tumor_output_fastq, normal_output_fastq,
                                 ref_genome, anonymizer, to_pair_anonymized_reads, mem_debug_writer)
                end2 = timer()
                logging.debug(
                    f'Time to anonymize reads in window {seq_name} {window[0]}-{window[1]} type {window[2].variant_type.name}: {end2 - start2}')
                DEBUG_TOTAL_TIMES['anonymize_windows'] += end2 - start2
    # Deal with any remaining anonymized reads, which have pairs or primary mappings in non-window regions
    memory_usage = process.memory_info().rss / (1024 * 1024)  # in MB
    mem_debug_writer.write(
        f'Memory usage after windows: {memory_usage} MB\n')
    ref_genome.close()
    logging.info('Searching for remaining unpaired anonymized reads')
    if to_pair_anonymized_reads:
        start3 = timer()
        pair_unpaired_or_supplementaries(to_pair_anonymized_reads, tumor_bam_file, normal_bam_file,
                                         tumor_output_fastq, normal_output_fastq, threads_per_file)
        end3 = timer()
        logging.debug(f'Time to pair unpaired anonymized reads: {end3 - start3}')
        DEBUG_TOTAL_TIMES['unpaired_searches'] += end3 - start3
        # If there are still objects in the collection, write them as single end reads
        memory_usage = process.memory_info().rss / (1024 * 1024)  # in MB
        mem_debug_writer.write(
            f'Memory usage after pairing leftover reads: {memory_usage} MB\n')
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


def run_short_read_tumor_normal_anonymizer(vcf_variants_per_sample: List[str],
                                           tumor_normal_samples: List[Tuple[str, str]],
                                           ref_genome: str, anonymizer: Anonymizer,
                                           output_filenames: List[Tuple[str, str]], cpus):
    """
    Anonymizes genomic sequencing from short read tumor-normal pairs, in the windows from each VCF variant
    Args:
        :param vcf_variants_per_sample: The VCF variants around which the sequencing data will be anonymized, for each sample
        :param tumor_normal_samples: The list of tumor-normal samples containing the reads to be anonymized. Each sample is a
            tuple containing the tumor and normal bam files in that order
        :param ref_genome: The reference genome to which the reads were mapped
        :param anonymizer: The specified anonymizing method
        :param output_filenames: The output filenames for the anonymized reads, in the same format as the input samples
        :param cpus: The number of cpus to use for the anonymization of each tumor-normal sample
    """
    # TODO: Implement multithreading for each sample pair
    tasks = []
    with ProcessPoolExecutor(max_workers=cpus) as executor:
        processes_by_sample = max(cpus // len(tumor_normal_samples), 1)
        for vcf_variants, samples, sample_output_files in zip(vcf_variants_per_sample, tumor_normal_samples,
                                                              output_filenames):
            tasks.append(executor.submit(anonymize_genome, vcf_variants, samples[0], samples[1], ref_genome, anonymizer,
                                         sample_output_files[0], sample_output_files[1], processes_by_sample))
        for task in as_completed(tasks):
            task.result()
