# @author: Nicolas Gaitan
import gzip
import logging
from typing import List, Tuple, Iterable
import pysam
from variant_extractor import VariantExtractor
from variant_extractor.variants import VariantType
from src.GenomeAnonymizer.anonymizer_methods import Anonymizer, AnonymizedRead
from src.GenomeAnonymizer.variation_classifier import DATASET_IDX_TUMORAL, DATASET_IDX_NORMAL
from timeit import default_timer as timer

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


def anonymize_genome(vcf_variants: VariantExtractor, tumor_bam: pysam.AlignmentFile, normal_bam: pysam.AlignmentFile,
                     ref_genome: pysam.FastaFile, anonymizer: Anonymizer, tumor_output_fastq: str,
                     normal_output_fastq: str):
    """
    Anonymizes genomic data using the provided VCF variants, normal and tumor BAM files, reference genome,
    classifier, and anonymizer object.
    """
    start1 = timer()
    windows = get_windows(vcf_variants, ref_genome)
    end1 = timer()
    print(f'Time to get windows: {end1 - start1} s')
    open(tumor_output_fastq + '.1.fastq', 'w').close()
    open(normal_output_fastq + '.1.fastq', 'w').close()
    open(tumor_output_fastq + '.2.fastq', 'w').close()
    open(normal_output_fastq + '.2.fastq', 'w').close()
    for seq_name, seq_windows in windows.items():
        for window in seq_windows:
            # Window is a tuple (start, end, VariantRecord)
            tumor_reads_pileup = tumor_bam.pileup(contig=seq_name, start=window[0], stop=window[1],
                                                  fastafile=ref_genome, min_base_quality=0)
            normal_reads_pileup = normal_bam.pileup(contig=seq_name, start=window[0], stop=window[1],
                                                    fastafile=ref_genome, min_base_quality=0)
            start2 = timer()
            anonymizer.anonymize(window[2], tumor_reads_pileup, normal_reads_pileup, ref_genome)
            end2 = timer()
            print(f'Time to anonymize reads in window {seq_name} {window[0]}-{window[1]} type {window[2].variant_type.name}: {end2 - start2} s')
            # Anonymized reads generated per window
            anonymized_reads_generator: Iterable[AnonymizedRead] = anonymizer.yield_anonymized_reads()
            start3 = timer()
            with (open(tumor_output_fastq + '.1.fastq', 'a') as tumor_fastq_writer_pair1,
                  open(normal_output_fastq + '.1.fastq', 'a') as normal_fastq_writer_pair1,
                  open(tumor_output_fastq + '.2.fastq', 'a') as tumor_fastq_writer_pair2,
                  open(normal_output_fastq + '.2.fastq', 'a') as normal_fastq_writer_pair2
                  ):
                for anonymized_read in anonymized_reads_generator:
                    fastq_record = anonymized_read.get_anonymized_read()
                    string_record = str(fastq_record)
                    dataset_idx = anonymized_read.dataset_idx
                    if dataset_idx == DATASET_IDX_TUMORAL:
                        if anonymized_read.is_read1:
                            tumor_fastq_writer_pair1.write(string_record)
                        elif anonymized_read.is_read2:
                            tumor_fastq_writer_pair2.write(string_record)
                    elif dataset_idx == DATASET_IDX_NORMAL:
                        if anonymized_read.is_read1:
                            normal_fastq_writer_pair1.write(string_record)
                        elif anonymized_read.is_read2:
                            normal_fastq_writer_pair2.write(string_record)
            end3 = timer()
            print(f'Time to write anonymized reads in window {seq_name} {window[0]}-{window[1]} type {window[2].variant_type.name}: {end3 - start3} s')
            anonymizer.reset()


def run_short_read_tumor_normal_anonymizer(vcf_variants_per_sample: List[VariantExtractor],
                                           tumor_normal_samples: List[Tuple[pysam.AlignmentFile, pysam.AlignmentFile]],
                                           ref_genome: pysam.FastaFile, anonymizer: Anonymizer,
                                           output_filenames: List[Tuple[str, str]]):
    """
    Anonymizes genomic sequencing from short read tumor-normal pairs, in the windows from each VCF variant
    Args:
        :param vcf_variants_per_sample: The VCF variants around which the sequencing data will be anonymized, for each sample
        :param tumor_normal_samples: The list of tumor-normal samples containing the reads to be anonymized. Each sample is a
            tuple containing the tumor and normal bam files in that order
        :param ref_genome: The reference genome to which the reads were mapped
        :param anonymizer: The specified anonymizing method
        :param output_filenames: The output filenames for the anonymized reads, in the same format as the input samples
    """
    # TODO: Implement multithreading for each sample pair
    for vcf_variants, samples, sample_output_files in zip(vcf_variants_per_sample, tumor_normal_samples,
                                                          output_filenames):
        anonymize_genome(vcf_variants, samples[0], samples[1], ref_genome, anonymizer, sample_output_files[0],
                         sample_output_files[1])
