import pathlib
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
import logging
import numpy as np
import pandas as pd
import pysam
from variant_extractor import VariantExtractor
from variant_extractor.variants import VariantRecord, VariantType
from src.GenomeAnonymizer.short_read_tumor_normal_anonymizer import Window, get_ref_idxs


def read_variation_windows(stats_files_list: list[str], ref_idx_dict) -> (dict[str, dict[str, list]],
                                                                          dict[str, list[Window]]):
    variation_per_window_by_seq: dict[str, dict[str, list]] = {k: {} for k in ref_idx_dict.keys()}
    # print(variation_per_window_by_seq.keys())
    window_order_lists: dict[str, list[Window]] = {k: [] for k in ref_idx_dict.keys()}
    for file in stats_files_list:
        with open(file, 'r') as f:
            line = f.readline()
            while line != '':
                if not line.startswith('#'):
                    elems = line.strip().split('\t')
                    # window_key = ','.join((elems[0], elems[1], elems[2]))
                    window = Window(sequence=elems[0], first=int(elems[1]), last=int(elems[2]))
                    counts_per_type = list(map(int, elems[3:]))
                    # print(window)
                    variation_per_window_in_seq = variation_per_window_by_seq.get(window.sequence)
                    variation_per_window_in_seq[str(window)] = counts_per_type
                    window_order_list = window_order_lists.get(window.sequence)
                    window_order_list.append(window)
                line = f.readline()
    # windows.sort(key=lambda x: (ref_sequences_dict.get(x.sequence), x.first, x.last))
    for _, window_order_list in window_order_lists.items():
        window_order_list.sort(key=lambda x: (ref_idx_dict.get(x.sequence), x.first, x.last))
    return variation_per_window_by_seq, window_order_lists


def compare(seq_idx1: int, first1: int, last1: int, seq_idx2: int, first2: int, last2: int) -> int:
    overlap = first2 <= last1 and last2 >= first1
    if seq_idx1 < seq_idx2:
        return -3
    if seq_idx1 > seq_idx2:
        return 3
    # For these cases seq_idx1 == seq_idx2
    if last1 < last2:
        return -1 if overlap else -2
    if last2 < last1:
        return 1 if overlap else 2
    # For these cases last1 == last2
    if first1 < first2:
        return -1
    if first2 < first1:
        return 1
    return 0


def process_variation_from_seq(file: str, window_order_lists, ref_idxs, min_af) -> (str, list):
    panel_variant_counts_in_seq: dict[str, list] = {}
    extractor = VariantExtractor(file)
    it = extractor.__iter__()
    variant: VariantRecord = next(it, None)
    if variant is None:
        raise ValueError(f'No variants in this file: {file}')
    sequence = variant.contig
    # windows_in_seq = anonymized_variation_per_window.get(sequence)
    windows_in_seq = window_order_lists.get(sequence)
    for window in windows_in_seq:
        window_panel_counts = [0] * len(VariantType)
        while variant is not None:
            cmp = compare(ref_idxs[variant.contig], variant.pos, variant.end, ref_idxs[window.sequence], window.first,
                          window.last)
            # logging.info(f'Variant chr: {variant.contig} idx: {ref_idxs[variant.contig]} {variant.pos} {variant.end} {ref_idxs[window.sequence]} {window.first} {window.last}'
            #             f' Window {str(window)} cmp={cmp}')
            if cmp < -1:
                variant = next(it, None)
            elif cmp > 1:
                break
            else:
                # if variant is None:
                #     logging.info(f'BUG: Variant in window: {str(window)} is None, in file: {file}, with'
                #                  f'cmp={cmp}')
                var_type_idx = variant.variant_type.value-1
                var_allele_freq = variant.info.get('AF', None)[0]
                # logging.info(f' pre-VAF: {var_allele_freq}')
                if var_allele_freq is None:
                    var_allele_freq = 0
                    logging.warning(f'Variant {str(variant)} does not have allele frequency (AF) field')
                # print(var_allele_freq)
                if var_allele_freq > min_af:
                    # logging.info(f'{len(window_panel_counts)}: {var_type_idx} from {variant.variant_type.name}')
                    window_panel_counts[var_type_idx] += 1
                    # logging.info(
                    #     f'Variant chr: {variant.contig} idx: {ref_idxs[variant.contig]} {variant.pos} {variant.end} {ref_idxs[window.sequence]} {window.first} {window.last}'
                    #     f' Window {str(window)} cmp={cmp}'
                    #     f' VAF: {var_allele_freq} - count: {window_panel_counts[var_type_idx]}')
                    # logging.info(f' VAF: {var_allele_freq} - count: {window_panel_counts[var_type_idx]}')
                variant = next(it, None)
        panel_variant_counts_in_seq[str(window)] = window_panel_counts
        # logging.info(f'Changing window {str(window)}')
    logging.info(f'Finished processing variants in panel sequence {sequence} from file: {file}')
    return sequence, panel_variant_counts_in_seq


def read_panel_variation(panel_files, window_order_lists, ref_idx_dict, min_af, cpus):
    panel_variant_counts: dict[str, dict[str, list]] = {k: {} for k in ref_idx_dict.keys()}
    executor = ProcessPoolExecutor(max_workers=min(len(panel_files), cpus))
    tasks = []
    for panel_file in panel_files:
        tasks.append(
            executor.submit(process_variation_from_seq, panel_file, window_order_lists, ref_idx_dict, min_af))
    for task in tasks:
        seq, panel_variant_counts_in_seq = task.result()
        # if panel_variant_counts[seq]:
        panel_variant_counts[seq] |= panel_variant_counts_in_seq
        # else:
        #     panel_variant_counts[seq] = panel_variant_counts_in_seq
    return panel_variant_counts


def results_to_dataframe(window_lists, anon_variation, panel_variation):
    data_rows = []
    columns = ['window_seq', 'window_first', 'window_last',
               'anon_SNV', 'anon_DEL', 'anon_INS', 'anon_DUP', 'anon_INV', 'anon_CNV', 'anon_TRA', 'anon_SGL',
               'panel_SNV', 'panel_DEL', 'panel_INS', 'panel_DUP', 'panel_INV', 'panel_CNV', 'panel_TRA', 'panel_SGL']
    for seq, windows_in_seq in window_lists.items():
        anon_windows_in_seq = anon_variation.get(seq)
        panel_windows_in_seq = panel_variation.get(seq)
        for window in windows_in_seq:
            data_row = [window.sequence, window.first, window.last]
            anon_counts = anon_windows_in_seq.get(str(window))
            panel_window_counts = panel_windows_in_seq.get(str(window))
            if panel_window_counts is None:
                panel_window_counts = [None] * len(VariantType)
            data_row += anon_counts + panel_window_counts
            data_rows.append(data_row)
    return pd.DataFrame(data_rows, columns=columns)


if __name__ == '__main__':
    parser = ArgumentParser(
        prog='Analysis script: Anonymized Variation in windows',
        description='Compare germinal variation in windows from anonymization to a reference panel', )

    parser.add_argument('-d', '--directory', type=str, help='Directory in which the statistic files'
                                                            'are located', required=True)
    parser.add_argument('-pd', '--panel_directory', type=str, help='Directory in which the reference panel files'
                                                                   'are located',
                        default='panel',
                        required=False)
    parser.add_argument('-c', '--cpu', help='Number of CPUs available for the execution', type=int,
                        required=False,
                        default=1)
    parser.add_argument('-r', '--reference', type=str, help='reference genome to which the reads are mapped',
                        required=True)
    parser.add_argument('--min_AF', type=float, help='minimum MAF to consider variants in the panel for quantification',
                        default=0.0,
                        required=False)
    parser.add_argument('--ignore_seqs', )
    logging.basicConfig(level=20)
    logging.info(f'Beginning analysis of anonymized variants')
    config = parser.parse_args()
    try:
        ref_genome = pysam.FastaFile(config.reference)
        ref_idx_sequences = get_ref_idxs(ref_genome)
        stats_files = list(map(str, pathlib.Path(config.directory).glob('*.statistics.txt')))
        for stats_file in stats_files:
            logging.info(f'Analyzing {stats_file}')
        panel_chr_files = list(
            map(str, pathlib.Path(config.directory).glob(f'{config.panel_directory}/*.haplotypes.vcf.gz')))
        for panel_chr_file in panel_chr_files:
            logging.info(f'Reading from panel {panel_chr_file}')
        # Dictionaries hashed by Windows in each sequence, with values as a np.array with the counts per variation type,
        # indexed by the VariationType enum value
        anonymized_variation, window_order = read_variation_windows(stats_files, ref_idx_sequences)
        window_germline_variants_from_panel: dict[str, dict[str, list]] = read_panel_variation(
            panel_chr_files, window_order, ref_idx_sequences, config.min_AF, config.cpu)
        results_df = results_to_dataframe(window_order, anonymized_variation, window_germline_variants_from_panel)
        results_df.to_csv(f'{config.directory}/anonymized_variation_vs_panel.csv', sep='\t', index=False)
    except Exception as e:
        logging.error(f'Error while analyzing anonymized variants: {e}')
        raise
