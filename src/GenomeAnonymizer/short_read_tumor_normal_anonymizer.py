# @author: Nicolas Gaitan
import gzip

import pysam
import os
from variant_extractor import VariantExtractor
from variant_extractor.variants import VariantRecord
from anonymizer_methods import Anonymizer
from variation_classifier import DATASET_IDX_TUMORAL, DATASET_IDX_NORMAL

"""Module for anonymizing short read tumor-normal pair genomes, using any anonymize function"""


def get_window(variant_record, window_size=2000):
    """
    A function that calculates a window around a variant record position.
    Returns:
    tuple: A tuple containing contig, start position, and end position of the window.
    """
    half_window = int(window_size / 2)
    return variant_record.contig, variant_record.pos - half_window, variant_record.end + half_window


def anonymize_genome(vcf_variants: VariantExtractor, tumor_bam: pysam.AlignmentFile, normal_bam: pysam.AlignmentFile,
                     ref_genome: pysam.Fastafile, anonymizer: Anonymizer, tumor_output_fastq: str,
                     normal_output_fastq: str):
    """
    Anonymizes genomic data using the provided VCF variants, normal and tumor BAM files, reference genome,
    classifier, and anonymizer object.
    """
    try:
        for variant_record in vcf_variants:
            window = get_window(variant_record)
            tumor_reads_pileup = tumor_bam.pileup(contig=window[0], start=window[1], stop=window[2],
                                                  fastafile=ref_genome)
            normal_reads_pileup = normal_bam.pileup(contig=window[0], start=window[1], stop=window[2],
                                                    fastafile=ref_genome)
            anonymizer.anonymize(variant_record, tumor_reads_pileup, normal_reads_pileup)
            anonymized_reads_generator = anonymizer.yield_anonymized_reads()
            # TODO: modify to account for read pairs identification
            with (gzip.open(tumor_output_fastq + '.gz', 'wb') as tumor_fastq_writer,
                  gzip.open(normal_output_fastq + '.gz', 'wb') as normal_fastq_writer):
                for dataset_idx, fastq_record in anonymized_reads_generator:
                    byte_coded_record = str(fastq_record).encode()
                    if dataset_idx == DATASET_IDX_TUMORAL:
                        tumor_fastq_writer.write(byte_coded_record)
                    elif dataset_idx == DATASET_IDX_NORMAL:
                        normal_fastq_writer.write(byte_coded_record)
    except Exception as e:
        print(e)


