# @author: Nicolas Gaitan
import gzip
import logging
from typing import Union, List, Tuple

import pysam
import os
from variant_extractor import VariantExtractor
from variant_extractor.variants import VariantRecord
from anonymizer_methods import Anonymizer
from variation_classifier import DATASET_IDX_TUMORAL, DATASET_IDX_NORMAL

"""Module for anonymizing short read tumor-normal pair genomes, using any object that implements the Anonymizer protocol"""


def get_window(variant_record, window_size=2000):
    """
    A function that calculates a window around a variant record position.
    Returns:
    tuple: A tuple containing contig, start position, and end position of the window.
    """
    half_window = int(window_size / 2)
    return variant_record.contig, variant_record.pos - half_window, variant_record.end + half_window


def anonymize_genome(vcf_variants: VariantExtractor, tumor_bam: pysam.AlignmentFile, normal_bam: pysam.AlignmentFile,
                     ref_genome: pysam.FastaFile, anonymizer: Anonymizer, tumor_output_fastq: str,
                     normal_output_fastq: str):
    """
    Anonymizes genomic data using the provided VCF variants, normal and tumor BAM files, reference genome,
    classifier, and anonymizer object.
    """
    try:
        for variant_record in vcf_variants:
            window = get_window(variant_record)
            tumor_reads_pileup = tumor_bam.pileup(contig=window[0], start=window[1], stop=window[2],
                                                  fastafile=ref_genome, min_base_quality=0)
            normal_reads_pileup = normal_bam.pileup(contig=window[0], start=window[1], stop=window[2],
                                                    fastafile=ref_genome, min_base_quality=0)
            anonymizer.anonymize(variant_record, tumor_reads_pileup, normal_reads_pileup)  # , ref_genome)
            anonymized_reads_generator = anonymizer.yield_anonymized_reads()
            with (gzip.open(tumor_output_fastq + '.1.fastq.gz', 'wb') as tumor_fastq_writer_pair1,
                  gzip.open(normal_output_fastq + '.1.fastq.gz', 'wb') as normal_fastq_writer_pair1,
                  gzip.open(tumor_output_fastq + '.2.fastq.gz', 'wb') as tumor_fastq_writer_pair2,
                  gzip.open(normal_output_fastq + '.2.fastq.gz', 'wb') as normal_fastq_writer_pair2
                  ):
                for dataset_idx, fastq_record in anonymized_reads_generator:
                    byte_coded_record = str(fastq_record).encode()
                    if dataset_idx == DATASET_IDX_TUMORAL:
                        if fastq_record.is_read1:
                            tumor_fastq_writer_pair1.write(byte_coded_record)
                        elif fastq_record.is_read2:
                            tumor_fastq_writer_pair2.write(byte_coded_record)
                        # tumor_fastq_writer.write(byte_coded_record)
                    elif dataset_idx == DATASET_IDX_NORMAL:
                        if fastq_record.is_read1:
                            normal_fastq_writer_pair1.write(byte_coded_record)
                        elif fastq_record.is_read2:
                            normal_fastq_writer_pair2.write(byte_coded_record)
                        # normal_fastq_writer.write(byte_coded_record)
    except Exception as e:
        logging.error(e)


def run_short_read_tumor_normal_anonymizer(vcf_variants: VariantExtractor,
                                           tumor_normal_samples: List[Tuple[pysam.AlignmentFile, pysam.AlignmentFile]],
                                           ref_genome: pysam.FastaFile, anonymizer: Anonymizer,
                                           output_filenames: List[Tuple[str, str]]):
    """
    Anonymizes genomic sequencing from short read tumor-normal pairs, in the windows from each VCF variant
    Args:
        :param vcf_variants: The VCF variants around which the sequencing data will be anonymized
        :param tumor_normal_samples: The list of tumor-normal samples containing the reads to be anonymized. Each sample is a
            tuple containing the tumor and normal bam files in that order
        :param ref_genome: The reference genome to which the reads were mapped
        :param anonymizer: The specified anonymizing method
        :param output_filenames: The output filenames for the anonymized reads, in the same format as the input samples
    """
    # TODO: Implement multithreading for each sample pair
    for sample, sample_output_files in zip(tumor_normal_samples, output_filenames):
        try:
            anonymize_genome(vcf_variants, sample[0], sample[1], ref_genome, anonymizer, sample_output_files[0],
                             sample_output_files[1])
        except Exception as e:
            logging.error(e)
