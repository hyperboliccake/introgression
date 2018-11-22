import sys
import os
import gzip
import gene_predictions
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc/')
import read_fasta

cen_starts = [151465, 238207, 114385, 449711, 151987, 148510, 496920, 105586, 355629, 436307, 440129, 150828, 268031, 628758, 326584, 555957]
cen_starts = [x-1 for x in cen_starts]

cen_ends =[151582,238323,114501,449821,152104,148627,497038,105703,355745,436425,440246,150947,268149,628875,326702,556073]
cen_ends = [x-1 for x in cen_ends]

tel_coords = [1,801,229411,230218,1,6608,812379,813184,1,1098,315783,316620,1,904,1524625,1531933,1,6473,569599,576874,1,5530,269731,270161,1,781,1083635,1090940,1,5505,556105,562643,1,7784,439068,439888,1,7767,744902,745751,1,807,665904,666816,1,12085,1064281,1078177,1,6344,923541,924431,1,7428,783278,784333,1,847,1083922,1091291,1,7223,942396,948010]
tel_coords = [x-1 for x in tel_coords]
tel_left_starts = [tel_coords[i] for i in range(0, len(tel_coords), 4)]
tel_left_ends = [tel_coords[i] for i in range(1, len(tel_coords), 4)]
tel_right_starts = [tel_coords[i] for i in range(2, len(tel_coords), 4)]
tel_right_ends = [tel_coords[i] for i in range(3, len(tel_coords), 4)]

def distance_from_telomere(start, end, chrm):

    assert start <= end, str(start) + ' ' + str(end)

    i = gp.chrms.index(chrm)
    # region entirely on left arm
    if end <= cen_starts[i]:
        return start - tel_left_ends[i]
    # region entirely on right arm
    if start >= cen_ends[i]:
        return tel_right_starts[i] - end
    # region overlaps centromere: return minimum distance from either telomere
    return min(start - tel_left_ends[i], tel_right_starts[i] - end)

def distance_from_centromere(start, end, chrm):

    assert start <= end, str(start) + ' ' + str(end)

    i = gp.chrms.index(chrm)
    # region entirely on left arm
    if end <= cen_starts[i]:
        return cen_starts[i] - end
    # region entirely on right arm
    if start >= cen_ends[i]:
        return start - cen_ends[i]
    # region overlaps centromere: return 0
    return 0

def write_region_summary_plus(fn, regions, fields):
    f = open(fn, 'w')
    f.write('region_id\t' + '\t'.join(fields) + '\n')
    keys = sorted(regions.keys(), key=lambda x: int(x[1:]))
    for region_id in keys:
        f.write(region_id + '\t')
        f.write('\t'.join([str(regions[region_id][field]) for field in fields]))
        f.write('\n')
    f.close()


def gap_columns(seqs):
    g = 0
    for i in range(len(seqs[0])):
        for seq in seqs:
            if seq[i] == gp.gap_symbol:
                g += 1
                break
    return g

def longest_consecutive(s, c):
    max_consecutive = 0
    current_consecutive = 0
    in_segment = False
    for i in range(len(s)):
        if s[i] == c:
            current_consecutive += 1
            in_segment = True
        else:
            if in_segment:
                max_consecutive = max(max_consecutive, current_consecutive)
                current_consecutive = 0
            in_segment = False
    return max_consecutive

#def in_intervals(i, intervals):
#    left = [x[0] for x in intervals]
#    right = [x[1] for x in intervals]
#    ind = bisect.bisect_right(left, i) - 1
#    if start < 0:
#        return False
#    start = left[ind]
#    end = right[ind]
#    assert i >= start
#    if i <= end:
#        return True

def masked_columns(seqs):
    # return two things:
    # - number of columns that are masked in any sequence
    # - above, but excluding columns with gaps
    num_seqs = len(seqs)
    num_sites = len(seqs[0])
    mask_total = 0
    mask_non_gap_total = 0
    for ps in range(num_sites):
        mask = False
        gap = False
        for s in range(num_seqs):
            if seqs[s][ps] == gp.gap_symbol:
                gap = True
            elif seqs[s][ps] == gp.masked_symbol:
                mask = True
        if mask:
            mask_total += 1
            if not gap:
                mask_non_gap_total += 1
    return mask_total, mask_non_gap_total

def index_by_reference(ref_seq, seq):
    # return dictionary keyed by reference index, with value the
    # corresponding index in non-reference sequence

    d = {}
    ri = 0
    si = 0
    for i in range(len(ref_seq)):
        if ref_seq[i] != gp.gap_symbol:
            d[ri] = si
            ri += 1
        if seq[i] != gp.gap_symbol:
            si += 1
    return d

def index_alignment_by_reference(ref_seq):
    # want a way to go from reference sequence coordinate to index in
    # alignment
    l = []
    for i in range(len(ref_seq)):
        if ref_seq[i] != gp.gap_symbol:
            l.append(i)
    return l
    

#def slice_alignment_by_reference(seq, ref_seq, ref_start, ref_end):
    
    
