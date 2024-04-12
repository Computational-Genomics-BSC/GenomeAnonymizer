# @author: Nicolas Gaitan
import gzip
import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Tuple, Iterable, Dict
import pysam
from variant_extractor import VariantExtractor
from variant_extractor.variants import VariantType
from src.GenomeAnonymizer.anonymizer_methods import Anonymizer, AnonymizedRead, add_anonymized_read_pair_to_collection
from src.GenomeAnonymizer.variation_classifier import DATASET_IDX_TUMORAL, DATASET_IDX_NORMAL, PAIR_1_IDX, PAIR_2_IDX

"""Module for anonymizing short read tumor-normal pair genomes, using any object that implements the Anonymizer 
protocol"""


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


def anonymized_read_pair_is_writeable(anonymized_read_pair1: AnonymizedRead,
                                      anonymized_read_pair2: AnonymizedRead) -> bool:
    if anonymized_read_pair1 is None or anonymized_read_pair2 is None:
        return False
    if anonymized_read_pair1.has_only_supplementary or anonymized_read_pair2.has_only_supplementary:
        return False
    return True


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


def anonymize_genome(vcf_variant_file: str, tumor_bam_file: str, normal_bam_file: str,
                     ref_genome_file: str, anonymizer: Anonymizer, tumor_output_fastq: str,
                     normal_output_fastq: str):
    """
    Anonymizes genomic data using the provided VCF variants, normal and tumor BAM files, reference genome,
    classifier, and anonymizer object.
    """
    vcf_variants = VariantExtractor(vcf_variant_file)
    tumor_bam = pysam.AlignmentFile(tumor_bam_file)
    normal_bam = pysam.AlignmentFile(normal_bam_file)
    ref_genome = pysam.FastaFile(ref_genome_file)
    windows = get_windows(vcf_variants, ref_genome)
    # This collection contains anonymized reads that are unpaired, because their pairs are not present in the same window
    to_pair_anonymized_reads: Dict[str, List[AnonymizedRead]] = dict()
    open(tumor_output_fastq + '.1.fastq', 'w').close()
    open(tumor_output_fastq + '.2.fastq', 'w').close()
    open(normal_output_fastq + '.1.fastq', 'w').close()
    open(normal_output_fastq + '.2.fastq', 'w').close()
    for seq_name, seq_windows in windows.items():
        for window in seq_windows:
            # Window is a tuple (start, end, VariantRecord)
            tumor_reads_pileup = tumor_bam.pileup(contig=seq_name, start=window[0], stop=window[1],
                                                  fastafile=ref_genome, min_base_quality=0)
            normal_reads_pileup = normal_bam.pileup(contig=seq_name, start=window[0], stop=window[1],
                                                    fastafile=ref_genome, min_base_quality=0)
            anonymizer.anonymize(window[2], tumor_reads_pileup, normal_reads_pileup, ref_genome)
            # Anonymized reads generated per window
            anonymized_reads_generator: Iterable[
                Tuple[AnonymizedRead, AnonymizedRead]] = anonymizer.yield_anonymized_reads()
            indexed_pair_writer_streams = [
                [open(tumor_output_fastq + '.1.fastq', 'a'), open(tumor_output_fastq + '.2.fastq', 'a')],
                [open(normal_output_fastq + '.1.fastq', 'a'), open(normal_output_fastq + '.2.fastq', 'a')]]
            for anonymized_read_pair in anonymized_reads_generator:
                anonymized_read_pair1 = anonymized_read_pair[PAIR_1_IDX]
                anonymized_read_pair2 = anonymized_read_pair[PAIR_2_IDX]
                if anonymized_read_pair_is_writeable(anonymized_read_pair1, anonymized_read_pair2):
                    write_pair(indexed_pair_writer_streams, anonymized_read_pair1, anonymized_read_pair2)
                else:
                    if anonymized_read_pair1 is not None:
                        anonymized_read_pair1 = anonymized_read_pair[PAIR_1_IDX]
                        dataset_idx = anonymized_read_pair1.dataset_idx
                        add_anonymized_read_pair_to_collection(to_pair_anonymized_reads,
                                                               anonymized_read_pair1.read_alignment,
                                                               dataset_idx)
                        read_id = anonymized_read_pair1.read_id
                    if anonymized_read_pair2 is not None:
                        anonymized_read_pair2 = anonymized_read_pair[PAIR_2_IDX]
                        dataset_idx = anonymized_read_pair2.dataset_idx
                        add_anonymized_read_pair_to_collection(to_pair_anonymized_reads,
                                                               anonymized_read_pair2.read_alignment,
                                                               dataset_idx)
                        read_id = anonymized_read_pair2.read_id
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
                    """if anonymized_read_pair1 is not None:
                        dataset_idx = anonymized_read_pair1.dataset_idx
                        fastq_record = str(anonymized_read_pair1.get_anonymized_fastq_record())
                    elif anonymized_read_pair2 is not None:
                        dataset_idx = anonymized_read_pair2.dataset_idx
                        fastq_record = str(anonymized_read_pair2.get_anonymized_fastq_record())
                    if dataset_idx == DATASET_IDX_TUMORAL:
                        tumor_fastq_writer_single_end.write(f'{fastq_record}\n')
                    elif dataset_idx == DATASET_IDX_NORMAL:
                        normal_fastq_writer_single_end.write(f'{fastq_record}\n')"""
            close_paired_streams(indexed_pair_writer_streams)
            # Reset the anonymizer for the next window
            anonymizer.reset()
    # Deal with any remaining anonymized reads, which have pairs or primary mappings in non-window regions
    if to_pair_anonymized_reads:
        tumor_reads_stream = pysam.AlignmentFile(tumor_bam_file).fetch(until_eof=True)
        normal_reads_stream = pysam.AlignmentFile(normal_bam_file).fetch(until_eof=True)
        indexed_pair_writer_streams = [
            [open(tumor_output_fastq + '.1.fastq', 'a'), open(tumor_output_fastq + '.2.fastq', 'a')],
            [open(normal_output_fastq + '.1.fastq', 'a'), open(normal_output_fastq + '.2.fastq', 'a')]]
        for dataset_idx, current_reads_stream in enumerate((tumor_reads_stream, normal_reads_stream)):
            for read_aln in current_reads_stream:
                read_id = read_aln.query_name
                if read_id in to_pair_anonymized_reads:
                    add_anonymized_read_pair_to_collection(to_pair_anonymized_reads, read_aln, dataset_idx)
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
        close_paired_streams(indexed_pair_writer_streams)
        # If there are still objects in the collection, write them as single end reads
        if to_pair_anonymized_reads:
            with (open(tumor_output_fastq + '.single_end.fastq', 'w') as tumor_fastq_writer_single_end,
                  open(normal_output_fastq + '.single_end.fastq', 'w') as normal_fastq_writer_single_end):
                for read_id, anonymized_read_pair in to_pair_anonymized_reads.items():
                    if anonymized_read_pair[PAIR_1_IDX] is not None:
                        single_anonymized_read = anonymized_read_pair[PAIR_1_IDX]
                        logging.warning('Single pair2 read not found for read id: ' + read_id)
                    elif anonymized_read_pair[PAIR_2_IDX] is not None:
                        single_anonymized_read = anonymized_read_pair[PAIR_2_IDX]
                        logging.warning('Single pair1 read not found for read id: ' + read_id)
                    if single_anonymized_read.has_left_overs_to_mask:
                        single_anonymized_read.mask_or_anonymize_left_over_variants()
                    dataset_idx = single_anonymized_read.dataset_idx
                    fastq_record = str(single_anonymized_read.get_anonymized_fastq_record())
                    if dataset_idx == DATASET_IDX_TUMORAL:
                        tumor_fastq_writer_single_end.write(f'{fastq_record}\n')
                    elif dataset_idx == DATASET_IDX_NORMAL:
                        normal_fastq_writer_single_end.write(f'{fastq_record}\n')


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
        for vcf_variants, samples, sample_output_files in zip(vcf_variants_per_sample, tumor_normal_samples,
                                                              output_filenames):
            tasks.append(executor.submit(anonymize_genome, vcf_variants, samples[0], samples[1], ref_genome, anonymizer,
                                         sample_output_files[0], sample_output_files[1]))
        for task in as_completed(tasks):
            task.result()
