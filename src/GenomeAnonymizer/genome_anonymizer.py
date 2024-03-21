import os
import re
import sys
from argparse import ArgumentParser
import logging
import pysam
from variant_extractor import VariantExtractor
from src.GenomeAnonymizer.anonymizer_methods import Anonymizer
from src.GenomeAnonymizer.short_read_tumor_normal_anonymizer import run_short_read_tumor_normal_anonymizer


def exec_parser():
    parser = ArgumentParser(
        prog='GenomeAnonymizer',
        description='Anonymization of sequencing data by removing germline variation',
        epilog='')
    parser.add_argument('-i', '--input', type=str, help='Input vcf file (.vcf, .vcf.gz)', required=True)
    parser.add_argument('-d', '--directory', type=str, help='Directory in which the tumor-normal sample pairs '
                                                            'and the samples text file are stored', required=True)
    parser.add_argument('-s', '--samples', type=str,
                        help='Text file with two columns with the tumor <first column> normal <second column> '
                             'sample pair file names, separated by tab', required=True)
    parser.add_argument('-r', '--reference', type=str, help='reference genome to which the mappings were done',
                        required=True)
    parser.add_argument('-m', '--method', type=str, help='anonymization method to apply on the samples ' 
                                                         'complete_germline: Mask all SNVs in the reads in the variant windows',
                        required=False,
                        default='complete_germline', choices=['complete_germline'])
    config = parser.parse_args()
    return config


def name_output(sample):
    output_suffix = '.anonymized'
    sample_name = re.sub('.bam | .sam | .cram', output_suffix, sample)
    return sample_name


def run_anonymizer():
    logging.basicConfig(level=logging.INFO)
    logging.info('Beginning execution')
    config = exec_parser()
    vcf_path = config.input
    variants = VariantExtractor(vcf_path)
    directory = config.directory
    ref_genome = pysam.FastaFile(config.reference)
    anonymizer_algorithm_name = config.method
    anonymizer_algorithm: Anonymizer = Anonymizer()
    if anonymizer_algorithm_name == 'complete_germline':
        from anonymizer_methods import CompleteGermlineAnonymizer
        anonymizer_algorithm: CompleteGermlineAnonymizer = CompleteGermlineAnonymizer()
    path_to_samples = ''.join((directory, '/', config.samples))
    logging.info('Reading inputs from %s', path_to_samples)
    samples = []
    output_samples = []
    with open(path_to_samples) as samples_file:
        for line in samples_file:
            sample_files = line.strip().split('\t')
            logging.info('Reading sample files %s and %s', sample_files[0], sample_files[1])
            samples_aln_files = (pysam.AlignmentFile(sample_files[0]), pysam.AlignmentFile(sample_files[1]))
            samples.append(samples_aln_files)
            logging.info('Anonymized samples will be written as %s and %s', sample_files[0], sample_files[1])
            samples_output_files = (name_output(sample_files[0]), name_output(sample_files[1]))
            output_samples.append(samples_output_files)
    logging.info('Beginning execution of the anonymizer algorithm')
    run_short_read_tumor_normal_anonymizer(variants, samples, ref_genome, anonymizer_algorithm, output_samples)
    logging.info('Finished execution of GenomeAnonymizer successfully')


if __name__ == "__main__":
    run_anonymizer()
