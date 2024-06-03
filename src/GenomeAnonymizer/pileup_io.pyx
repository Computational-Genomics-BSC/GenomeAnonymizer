cimport cython
import pysam
from pysam.libcalignmentfile cimport AlignmentFile
from pysam.libcalignedsegment cimport pysam_bam_get_cigar, pysam_bam_get_seq, \
    pysam_bam_get_qual, PileupColumn


def iter_pileups(alignment_file1, alignment_file2, ref_genome, seq_name, start, stop):
    cdef:
        PileupColumn p1
        PileupColumn p2
    it1 = alignment_file1.pileup(reference=seq_name, start=start, end=stop, fastafile=ref_genome, min_base_quality=0,
                                 min_mapping_quality=0, max_depth=100000, stepper='nofilter')
    it2 = alignment_file2.pileup(reference=seq_name, start=start, end=stop, fastafile=ref_genome, min_base_quality=0,
                                 min_mapping_quality=0, max_depth=100000, stepper='nofilter')
    p1 = next(it1, None)
    p2 = next(it2, None)
    while p1 is not None and p2 is not None:
        if p1.reference_pos < p2.reference_pos or p2 is None:
            yield p1, None
            p1 = next(it1, None)
        elif p1.reference_pos > p2.reference_pos or p1 is None:
            yield None, p2
            p2 = next(it2, None)
        else:
            yield p1, p2
            p1 = next(it1, None)
            p2 = next(it2, None)