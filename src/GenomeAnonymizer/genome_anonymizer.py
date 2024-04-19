import re
import sys
import os
from argparse import ArgumentParser
import logging
from typing import Tuple
from src.GenomeAnonymizer.anonymizer_methods import CompleteGermlineAnonymizer
from src.GenomeAnonymizer.short_read_tumor_normal_anonymizer import run_short_read_tumor_normal_anonymizer
from timeit import default_timer as timer


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
    parser.add_argument('-c', '--cpu', help='Number of CPUs available for the execution', type=int, default=1)
    config = parser.parse_args()
    return config


def name_output(sample):
    output_suffix = '.anonymized'
    sample_name = re.sub('.bam|.sam|.cram', output_suffix, sample)
    return sample_name


def join_dir_file(directory, param):
    return ''.join((directory, '/', param)) if not directory.endswith('/') else ''.join(
        (directory, param))


def run_anonymizer():
    # logging.basicConfig(level=logging.INFO)
    start1 = timer()
    logging.basicConfig(level=logging.DEBUG)
    logging.info('Beginning execution of GenomeAnonymizer v 0.0.1')
    config = exec_parser()
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
