import re
import sys
import os
from argparse import ArgumentParser
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Tuple, List
from src.GenomeAnonymizer.anonymizer_methods import CompleteGermlineAnonymizer, Anonymizer
# from src.GenomeAnonymizer.short_read_tumor_normal_anonymizer import run_short_read_tumor_normal_anonymizer
from timeit import default_timer as timer

from src.GenomeAnonymizer.short_read_tumor_normal_anonymizer import anonymize_genome

# Anonymizer algorithm options
COMPLETE_GERMLINE_ANONYMIZER_ALGORITHM = 'complete_germline'

ANONYMIZER_ALGORITHMS = set()
ANONYMIZER_ALGORITHMS.add(COMPLETE_GERMLINE_ANONYMIZER_ALGORITHM)


def exec_parser():
    parser = ArgumentParser(
        prog='GenomeAnonymizer',
        description='Anonymization of sequencing data by removing germline variation',
        epilog='')
    parser.add_argument('-d', '--directory', type=str, help='Directory in which the tumor-normal sample pairs '
                                                            'and the samples text file are stored', required=True)
    parser.add_argument('-s', '--samples', type=str,
                        help='Text file with two columns with the tumor <first column>, normal <second column> '
                             'sample pair file names, and the vcf files <third column> from each sample, separated by tab </\t>',
                        required=True)
    parser.add_argument('-r', '--reference', type=str, help='reference genome to which the reads are mapped',
                        required=True)
    parser.add_argument('-m', '--method', type=str, help='anonymization method to apply on the samples '
                                                         'complete_germline: Mask all SNVs in the reads in the variant windows',
                        required=False,
                        default='complete_germline', choices=['complete_germline'])
    parser.add_argument('-c', '--cpu', help='Number of CPUs available for the execution', type=int,
                        required=False,
                        default=1)
    parser.add_argument('-v', '--verbose', help='Verbosity of logging', type=int,
                        required=False,
                        default=2)
    config = parser.parse_args()
    return config


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


def name_output(sample):
    output_suffix = '.anonymized'
    sample_name = re.sub('.bam|.sam|.cram', output_suffix, sample)
    return sample_name


def join_dir_file(directory, param):
    return ''.join((directory, '/', param)) if not directory.endswith('/') else ''.join(
        (directory, param))


def run_anonymizer():
    config = exec_parser()
    verbosity = config.verbose * 10
    logging.basicConfig(level=verbosity)
    # logging.basicConfig(level=logging.INFO)
    start1 = timer()
    # logging.basicConfig(level=logging.DEBUG)
    logging.info('Beginning execution of GenomeAnonymizer v 0.0.2')
    variants_per_sample = []
    directory = config.directory
    ref_genome = config.reference
    anonymizer_algorithm_name = config.method
    if anonymizer_algorithm_name not in ANONYMIZER_ALGORITHMS:
        logging.error('Anonymizer algorithm %s is not a valid option', anonymizer_algorithm_name)
        sys.exit(1)
    anonymizer_algorithm: CompleteGermlineAnonymizer = CompleteGermlineAnonymizer()  # Default anonymizer algorithm
    path_to_samples = join_dir_file(directory, config.samples)
    logging.info('Reading inputs from %s', path_to_samples)
    samples = []
    output_samples = []
    try:
        with open(path_to_samples) as samples_file:
            for line in samples_file:
                if line.startswith('#'):
                    continue
                sample_files = line.strip().split('\t')
                tumor_sample = join_dir_file(directory, sample_files[0])
                normal_sample = join_dir_file(directory, sample_files[1])
                vcf_sample = join_dir_file(directory, sample_files[2])
                logging.info('Reading sample files %s and %s', tumor_sample, normal_sample)
                samples_aln_files: Tuple[str, str] = (tumor_sample, normal_sample)
                samples.append(samples_aln_files)
                logging.info('Reading vcf sample %s', vcf_sample)
                variants_per_sample.append(vcf_sample)
                tumor_output_prefix = name_output(tumor_sample)
                normal_output_prefix = name_output(normal_sample)
                logging.info('Anonymized samples will be written as %s and %s', tumor_output_prefix,
                             normal_output_prefix)
                samples_output_files = (tumor_output_prefix, normal_output_prefix)
                output_samples.append(samples_output_files)
        logging.info('Beginning execution of the anonymizer algorithm')
        run_short_read_tumor_normal_anonymizer(variants_per_sample, samples, ref_genome, anonymizer_algorithm,
                                               output_samples, config.cpu)
        logging.info('Finished execution of GenomeAnonymizer successfully')
    except Exception as e:
        logging.error(e)
        raise
    end1 = timer()
    logging.debug(f'Total execution time: {end1 - start1} s')


if __name__ == "__main__":
    run_anonymizer()
