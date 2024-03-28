# @author: Nicolas Gaitan

from typing import Protocol, Dict, List
import pysam
from variant_extractor.variants import VariantType
from src.GenomeAnonymizer.variants import CalledGenomicVariant, SomaticVariationType
from src.GenomeAnonymizer.variation_classifier import classify_variation_in_pileup_column
from src.GenomeAnonymizer.variation_classifier import DATASET_IDX_TUMORAL, DATASET_IDX_NORMAL
from timeit import default_timer as timer


class Anonymizer(Protocol):
    """
    Protocol that dictates the methods an anonymizer class must implement.
    """

    def anonymize(self, variant_record, tumor_reads_pileup, normal_reads_pileup, ref_genome) -> None:
        pass

    def mask_germline_snvs(self, pileup_column, called_genomic_variants,
                           variant_record) -> None:  # , ref_genome) -> None:
        pass

    def get_anonymized_reads(self):
        pass

    def yield_anonymized_reads(self):
        pass

    def reset(self):
        pass


class CompleteGermlineAnonymizer:
    def __init__(self):
        self.anonymized_reads: Dict[str, AnonymizedRead] = dict()
        self.has_anonymized_reads = False

    def reset(self):
        self.anonymized_reads = dict()
        self.has_anonymized_reads = False

    def anonymize(self, variant_record, tumor_reads_pileup, normal_reads_pileup, ref_genome):
        called_genomic_variants = {}
        seen_read_alns = set()
        for dataset_idx, current_pileup in enumerate((tumor_reads_pileup, normal_reads_pileup)):
            for pileup_column in current_pileup:
                start1 = timer()
                classify_variation_in_pileup_column(pileup_column, dataset_idx, seen_read_alns, ref_genome, called_genomic_variants)
                end1 = timer()
                print(f"Classify variation time: {end1 - start1}")
                # TODO: mask indels
                start2 = timer()
                self.mask_germline_snvs(pileup_column, dataset_idx, called_genomic_variants, variant_record)
                end2 = timer()
                print(f"Mask germline snvs time: {end2 - start2}")
                # , ref_genome)
        self.has_anonymized_reads = True

    def mask_germline_snvs(self, pileup_column, dataset_idx, called_genomic_variants, variant_record):  # , ref_genome):
        pos = pileup_column.reference_pos
        # ref_base = -1
        # ref_base = pileup_column.reference_base ref_base = ref_genome.fetch(pileup_column, pos - 1, pos) -> Might
        # cause issues due to accesing the same file again, after classification
        for pileup_read in pileup_column.pileups:
            aln = pileup_read.alignment
            # if ref_base == -1:  # Requires the MD tag to be present in the alignment
            #    if pos > len(aln.get_reference_sequence()):
            #        raise ValueError(f"Position {pos} is out of bounds of reference sequence {aln.get_reference_sequence()}")
            #    ref_base = aln.get_reference_sequence()[pos]
            if aln.query_name not in self.anonymized_reads:
                self.anonymized_reads[aln.query_name] = AnonymizedRead(aln, dataset_idx)
        variants_in_column: List[CalledGenomicVariant] = called_genomic_variants.get(pos, None)
        if variants_in_column is None:
            return
        variant_to_keep = CalledGenomicVariant.from_variant_record(variant_record)
        for called_variant in variants_in_column:
            if (called_variant.somatic_variation_type == SomaticVariationType.TUMORAL_NORMAL_VARIANT
                    and called_variant != variant_to_keep):
                if called_variant.variant_type == VariantType.SNV:
                    # TODO: mask indels
                    for read_id, var_read_pos in called_variant.supporting_reads.items():
                        self.anonymized_reads[read_id].mask_or_modify_base_pair(var_read_pos, called_variant.ref_allele)

    def get_anonymized_reads(self):
        fastx_records = ([], [])
        if not self.has_anonymized_reads:
            raise ValueError("No reads have been anonymized, call anonymize() first")
        for read_id, anonymized_read_obj in self.anonymized_reads.items():
            dataset_records = fastx_records[anonymized_read_obj.dataset_idx]
            dataset_records.append(anonymized_read_obj)
        return fastx_records

    def yield_anonymized_reads(self):
        if not self.has_anonymized_reads:
            raise ValueError("No reads have been anonymized, call anonymize() first")
        for read_id, anonymized_read_obj in self.anonymized_reads.items():
            yield anonymized_read_obj


class AnonymizedRead:
    def __init__(self, read_alignment: pysam.AlignedSegment, dataset_idx: int):
        self.read_id = read_alignment.query_name
        self.is_read1 = read_alignment.is_read1
        self.is_read2 = read_alignment.is_read2
        self.read_alignment = read_alignment
        self.original_sequence = read_alignment.query_sequence
        self.original_qualities = read_alignment.query_qualities
        self.anonymized_sequence_list = list(self.original_sequence)
        self.anonymized_qualities_list = list(self.original_qualities)
        self.dataset_idx = dataset_idx

    def mask_or_modify_base_pair(self, pos_in_read: int, new_base: str, modify_qualities: bool = False, new_quality: int = 0):
        self.anonymized_sequence_list[pos_in_read] = new_base
        if modify_qualities:
            self.anonymized_qualities_list[pos_in_read] = new_quality

    def get_anonymized_read(self) -> pysam.FastxRecord:
        anonymized_read_seq = ''.join(self.anonymized_sequence_list)
        anonimized_read_qual = ''.join(map(lambda x: chr(x + 33), self.anonymized_qualities_list))
        # anonimized_read_qual = ''.join(self.anonymized_qualities_list)
        if len(anonymized_read_seq) != len(self.original_sequence):
            raise ValueError("Anonymized read length does not match original read length")
        if len(anonimized_read_qual) != len(self.original_qualities):
            raise ValueError("Anonymized qualities length does not match original qualities length")
        anonymized_read: pysam.FastxRecord = pysam.FastxRecord(name=self.read_id, sequence=anonymized_read_seq,
                                                               quality=anonimized_read_qual)
        return anonymized_read
