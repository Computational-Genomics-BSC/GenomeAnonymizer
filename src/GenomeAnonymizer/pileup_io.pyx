cimport cython
import pysam
from pysam.libcalignmentfile cimport AlignmentFile, IteratorRow
from pysam.libcalignedsegment cimport pysam_bam_get_cigar, pysam_bam_get_seq, \
    pysam_bam_get_qual, PileupColumn, AlignedSegment
from pysam.libcalignedsegment import AlignedSegment

def iter_pileups(alignment_file1, alignment_file2, ref_genome, seq_name, start, stop):
    cdef:
        PileupColumn p1
        PileupColumn p2
    it1 = alignment_file1.pileup(reference=seq_name, start=start, end=stop, fastafile=ref_genome, min_base_quality=0,
                                 min_mapping_quality=0, max_depth=1000000, stepper='nofilter', ignore_overlaps=False,
                                 ignore_orphans=False)
    it2 = alignment_file2.pileup(reference=seq_name, start=start, end=stop, fastafile=ref_genome, min_base_quality=0,
                                 min_mapping_quality=0, max_depth=1000000, stepper='nofilter', ignore_overlaps=False,
                                 ignore_orphans=False)
    p1 = next(it1, None)
    p2 = next(it2, None)
    # while p1 is not None and p2 is not None:
    while True:
        if p1 is not None and p2 is not None:
            if p1.reference_pos < p2.reference_pos:
                yield p1, None
                p1 = next(it1, None)
            elif p1.reference_pos > p2.reference_pos:
                yield None, p2
                p2 = next(it2, None)
            else:
                yield p1, p2
                p1 = next(it1, None)
                p2 = next(it2, None)
        elif p1 is None and p2 is None:
            break
        else:
            if p2 is None:
                yield p1, None
                p1 = next(it1, None)
            elif p1 is None:
                yield None, p2
                p2 = next(it2, None)


cdef compare(seq_idx1: int, first1: int, last1: int, seq_idx2: int, first2: int, last2: int):
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


cdef compare_read_alignments_intersection(r1, r2):
    if r1 is None:
        return -4
    if r2 is None:
        return 4
    seq1 = r1.reference_id
    seq2 = r2.reference_id
    first1 = r1.reference_start
    first2 = r2.reference_start
    last1 = r1.reference_end if r1.is_mapped else first1
    last2 = r2.reference_end if r2.is_mapped else first2
    return compare(seq1, first1, last1, seq2, first2, last2)


# def collect_intersecting_reads(it, r_aln, arr, unmapped, p_right, read_ids):
def collect_intersecting_reads(it, r_aln, arr, unmapped, p_right):
    right = p_right
    next_aln = r_aln
    # print(f'Before cycle')
    while True:
        next_aln = next(it, None)
        # if next_aln is not None:
            # print(f'aln: {next_aln.to_string()}')
        if next_aln is None:
            # print(f'Breaks for None status')
            break
        # print(f'{r_aln.to_string()}')
        # if next_aln.query_name not in read_ids:
            # print(f'Read not in read ids')
        #     continue
        if next_aln.is_unmapped:
            # print('unmapped read appends')
            unmapped.append(next_aln)
            continue
        intersects = -1 <= compare_read_alignments_intersection(arr[-1], next_aln) <= 1
        # if next_aln.query_name == 'HWI-ST1133:223:C1940ACXX:4:2314:12959:99156':
        #     print(f'Tracked read pair 2={next_aln.is_read2}: intersects={intersects} ----- {next_aln.to_string()}')
        if not intersects:
            # print(f'Breaks for not intersecting')
            return next_aln
        arr.append(next_aln)
        # if next_aln.query_name == 'HWI-ST1133:223:C1940ACXX:4:2314:12959:99156' and next_aln.is_read2:
        #     print(f'Tracked read appends ----- {next_aln.to_string()}')
    return next_aln


cdef get_righmost_pos(arr, prev_right):
    right = prev_right
    if right is None:
        right = 0
    for aln in arr:
        if aln.is_mapped:
            # if aln.reference_end is None:
            #     print(f'Aln reference end is None for aln: {aln.to_string()}')
            # if right is None:
            #     print(f'Right value is None compared to aln: {aln.to_string()}')
            right = max(right, aln.reference_end)
    return right


# def iter_fetch_pair(alignment_file1, alignment_file2, read_ids, ref_genome, seq=None, first=None, last=None):
def iter_fetch_pair(alignment_file1, alignment_file2, ref_genome, seq=None, first=None, last=None):
    cdef:
        IteratorRow it1
        IteratorRow it2
        AlignedSegment r1
        AlignedSegment r2
        list r1_array
        list r2_array
        list r1_unmapped
        list r2_unmapped
    # track = seq == '3' and first == 146666756 and last == 198022429
    until_eof = True
    if seq is not None:
        until_eof = False
    it1 = alignment_file1.fetch(reference=seq, start=first, stop=last, until_eof=until_eof)
    it2 = alignment_file2.fetch(reference=seq, start=first, stop=last, until_eof=until_eof)
    r1_array = []
    r2_array = []
    r1_unmapped = []
    r2_unmapped = []
    r1 = next(it1, None)
    r2 = next(it2, None)
    r1_yielded = True
    r2_yielded = True
    if r1 is None and r2 is None:
        # if track:
        #     print('init r1 and 2 are none')
        return None
    if r1 is not None:
        # if track:
        #     print('init r1 is not none')
        seq1 = r1.reference_id
        seq_name1 = r1.reference_name
        left1 = r1.reference_start
        right1 = r1.reference_end
        # while r1.query_name not in read_ids:
        #     r1 = next(it1, None)
        r1_array.append(r1)
    if r2 is not None:
        # if track:
        #     print('init r2 is not none')
        seq2 = r2.reference_id
        seq_name2 = r2.reference_name
        left2 = r2.reference_start
        right2 = r2.reference_end
        # while r2.query_name not in read_ids:
        #     r2 = next(it2, None)
        r2_array.append(r2)
    while True:
        if r1_yielded and r1 is not None:
            # r1 = collect_intersecting_reads(it1, r1, r1_array, r1_unmapped, right1, read_ids)
            r1 = collect_intersecting_reads(it1, r1, r1_array, r1_unmapped, right1)
            right1 = get_righmost_pos(r1_array, right1)
            r1_yielded = False
            # if track:
            #     print(f'r1 colected')
        if r2_yielded and r2 is not None:
            # r2 = collect_intersecting_reads(it2, r2, r2_array, r2_unmapped, right2, read_ids)
            r2 = collect_intersecting_reads(it2, r2, r2_array, r2_unmapped, right2)
            right2 = get_righmost_pos(r2_array, right2)
            r2_yielded = False
            # if track:
            #     print(f'r2 colected')
        if r1 is None and r2 is None:
            yield r1_array, None, None
            yield None, r2_array, None
            break
        elif r1 is not None and r2 is not None:
            inter_cmp = compare(seq1, left1, right1, seq2, left2, right2)
            # if track:
            #     print(f'INTER_CMP={inter_cmp}')
            # if r2 is None or inter_cmp < -1:
            if inter_cmp < -1:
                # if r1_array[0].query_name == 'HWI-ST142:601:D1D1NACXX:1:1307:16993:85790' and r1_array[0].is_read1:
                #     print(f'Yielded list of Tracked read in T with len: {len(r1_array)}')
                # if r1 is None:
                #     continue
                # if track:
                #     print(f'r1 yielded')
                yield r1_array, None, None
                # if r1.query_name == 'HWI-ST142:601:D1D1NACXX:1:1307:16993:85790' and r1.is_read1:
                #     print('Yielded previous list to Tracked read in T')
                r1_yielded = True
                r1_array = [r1]
                # if r1_array[0].query_name == 'HWI-ST142:601:D1D1NACXX:1:1307:16993:85790' and r1_array[0].is_read1:
                #     print(f'New list is tracked read {r1_array[0].to_string()}')
                seq1 = r1.reference_id
                seq_name1 = r1.reference_name
                left1 = r1.reference_start
                right1 = r1.reference_end
            # elif r1 is None or inter_cmp > 1:
            elif inter_cmp > 1:
                # if r2_array[0].query_name == 'HWI-ST1133:217:D1D4WACXX:1:2103:3953:77371' and r2_array[0].is_read2:
                #     print(f'Yielded list of Tracked read in N with len: {len(r2_array)}')
                # if r2 is None:
                #     continue
                # if track:
                #     print(f'r2 yielded')
                yield None, r2_array, None
                # if r2.query_name == 'HWI-ST1133:217:D1D4WACXX:1:2103:3953:77371' and r2.is_read2:
                #     print('Yielded previous list to Tracked read in N')
                r2_yielded = True
                r2_array = [r2]
                # if r2_array[0].query_name == 'HWI-ST1133:217:D1D4WACXX:1:2103:3953:77371' and r2_array[0].is_read2:
                #     print(f'New list is tracked read {r2_array[0].to_string()}')
                seq2 = r2.reference_id
                seq_name2 = r2.reference_name
                left2 = r2.reference_start
                right2 = r2.reference_end
            else:
                # if r1_array[0].query_name == 'HWI-ST142:601:D1D1NACXX:1:1307:16993:85790' and r1_array[0].is_read1:
                #     print(f'Yielded list of Tracked read in T-N pair with len: {len(r1_array)} with pileup coordinates: '
                #           f'{(seq_name1, min(left1, left2), max(right1, right2))}')
                # if len(r2_array) < 5:
                #     for rt2 in r2_array:
                #         if rt2.query_name == 'HWI-ST142:601:D1D1NACXX:1:1307:16993:85790' and rt2.is_read1:
                #             print(
                #                 f'Yielded list of Tracked read in T-N pair with len: {len(r2_array)} with pileup coordinates: '
                #                 f'{(seq_name1, min(left1, left2), max(right1, right2))}')
                # if track:
                #     print(f'both yielded')
                yield r1_array, r2_array, (seq_name1, min(left1, left2), max(right1, right2))
                r1_yielded = True
                r2_yielded = True
                # if r1 is not None:
                    # if r1.query_name == 'HWI-ST142:601:D1D1NACXX:1:1307:16993:85790' and r1.is_read1:
                    #     print('Yielded previous list to Tracked read in T-N pair')
                r1_array = [r1]
                seq1 = r1.reference_id
                seq_name1 = r1.reference_name
                left1 = r1.reference_start
                right1 = r1.reference_end
                    # if r1_array[0].query_name == 'HWI-ST142:601:D1D1NACXX:1:1307:16993:85790' and r1_array[0].is_read1:
                    #     print('New list is tracked read')
                # if r2 is not None:
                r2_array = [r2]
                seq2 = r2.reference_id
                seq_name2 = r2.reference_name
                left2 = r2.reference_start
                right2 = r2.reference_end
        else:
            if r1 is not None:
                # if track:
                #     print(f'r1 yielded')
                yield r1_array, None, None
                # if r1.query_name == 'HWI-ST142:601:D1D1NACXX:1:1307:16993:85790' and r1.is_read1:
                #     print('Yielded previous list to Tracked read in T')
                r1_yielded = True
                r1_array = [r1]
                # if r1_array[0].query_name == 'HWI-ST142:601:D1D1NACXX:1:1307:16993:85790' and r1_array[0].is_read1:
                #     print(f'New list is tracked read {r1_array[0].to_string()}')
                seq1 = r1.reference_id
                seq_name1 = r1.reference_name
                left1 = r1.reference_start
                right1 = r1.reference_end
            if r2 is not None:
                yield None, r2_array, None
                # if r2.query_name == 'HWI-ST1133:217:D1D4WACXX:1:2103:3953:77371' and r2.is_read2:
                #     print('Yielded previous list to Tracked read in N')
                r2_yielded = True
                r2_array = [r2]
                # if r2_array[0].query_name == 'HWI-ST1133:217:D1D4WACXX:1:2103:3953:77371' and r2_array[0].is_read2:
                #     print(f'New list is tracked read {r2_array[0].to_string()}')
                seq2 = r2.reference_id
                seq_name2 = r2.reference_name
                left2 = r2.reference_start
                right2 = r2.reference_end
        """if r1 is not None:
            if track:
                print(f'r1 is not none {r1.to_string()}')
        if r2 is not None:
            if track:
                print(f'r2 is not none {r2.to_string()}')"""
    # print(f'Unmapped being yielded: r1={len(r1_unmapped)} and r2={len(r2_unmapped)}')
    yield None, None, (r1_unmapped, r2_unmapped)