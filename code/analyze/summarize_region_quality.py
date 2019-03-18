import bisect
import gzip
import global_params as gp
from misc import binary_search
import numpy as np

cen_starts = [151465, 238207, 114385, 449711, 151987, 148510,
              496920, 105586, 355629, 436307, 440129, 150828,
              268031, 628758, 326584, 555957]
cen_starts = [x-1 for x in cen_starts]

cen_ends = [151582, 238323, 114501, 449821, 152104, 148627,
            497038, 105703, 355745, 436425, 440246, 150947,
            268149, 628875, 326702, 556073]
cen_ends = [x-1 for x in cen_ends]

tel_coords = [1, 801, 229411, 230218,
              1, 6608, 812379, 813184,
              1, 1098, 315783, 316620,
              1, 904, 1524625, 1531933,
              1, 6473, 569599, 576874,
              1, 5530, 269731, 270161,
              1, 781, 1083635, 1090940,
              1, 5505, 556105, 562643,
              1, 7784, 439068, 439888,
              1, 7767, 744902, 745751,
              1, 807, 665904, 666816,
              1, 12085, 1064281, 1078177,
              1, 6344, 923541, 924431,
              1, 7428, 783278, 784333,
              1, 847, 1083922, 1091291,
              1, 7223, 942396, 948010]
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
    return np.where(ref_seq != gp.gap_symbol)[0]


def num_sites_between(sites, start, end):
    # sites are sorted
    i = bisect.bisect_left(sites, start)
    j = bisect.bisect_right(sites, end)
    return j - i, sites[i:j]


def read_masked_intervals(fn):
    with open(fn, 'r') as reader:
        reader.readline()  # header
        ints = []
        for line in reader:
            line = line.split()
            ints.append((int(line[0]), int(line[2])))

    return ints


def convert_intervals_to_sites(ints):
    sites = []
    for start, end in ints:
        sites += range(start, end + 1)
    return np.array(sites)


def seq_id_hmm(seq1, seq2, offset, include_sites):
    sites = np.array(include_sites) - offset

    info_gap = np.logical_or(seq1 == gp.gap_symbol,
                             seq2 == gp.gap_symbol)
    info_unseq = np.logical_or(seq1 == gp.unsequenced_symbol,
                               seq2 == gp.unsequenced_symbol)
    info_match = seq1 == seq2
    info_hmm = np.zeros(info_match.shape, bool)
    sites = sites[np.logical_and(sites < len(info_match), sites >= 0)]
    info_hmm[sites] = True

    total_sites = np.sum(info_hmm)
    total_moatch = np.sum(np.logical_and(info_hmm, info_match))

    # check all included are not gapped or skipped
    include_in_skip = np.logical_and(
        info_hmm, np.logical_or(
            info_unseq, info_gap))
    if np.any(include_in_skip):
        ind = np.where(include_in_skip)[0][0]
        raise AssertionError(f'{seq1[ind]} {seq2[ind]} {ind}')

    return total_match, total_sites, \
        {'gap_flag': info_gap, 'unseq_flag': info_unseq,
         'hmm_flag': info_hmm, 'match': info_match}


def seq_id_unmasked(seq1, seq2, offset, exclude_sites1, exclude_sites2):
    # total_sites is number of sites at which neither sequence is
    # masked or has a gap or unsequenced character; total_match is the
    # number of those sites at which the two sequences match
    # gapped and unsequenced locations
    info_gap = np.logical_or(seq1 == gp.gap_symbol,
                             seq2 == gp.gap_symbol)
    info_unseq = np.logical_or(seq1 == gp.unsequenced_symbol,
                               seq2 == gp.unsequenced_symbol)

    # convert offset excluded sites to boolean array
    info_mask = np.zeros(seq1.shape, bool)
    if exclude_sites1 != []:
        sites1 = np.array(exclude_sites1) - offset
        sites1 = sites1[np.logical_and(sites1 < len(info_gap),
                                       sites1 >= 0)]
        info_mask[sites1] = True
    if exclude_sites2 != []:
        sites2 = np.array(exclude_sites2) - offset
        sites2 = sites2[np.logical_and(sites2 < len(info_gap),
                                       sites2 >= 0)]
        info_mask[sites2] = True

    # find sites that are not masked, gapped, or unsequenced
    sites = np.logical_not(
        np.logical_or(
            info_mask,
            np.logical_or(
                info_gap, info_unseq)))

    # determine totals
    total_sites = np.sum(sites)
    total_match = np.sum(
        np.logical_and(
            seq1 == seq2,
            sites))

    return total_match, total_sites, {'mask_flag': info_mask}

    n = len(seq1)
    total_sites = 0
    total_match = 0

    skip = [gp.gap_symbol, gp.unsequenced_symbol]
    info_mask = [False for i in range(n)]
    for i in range(n):

        if binary_search.present(exclude_sites1, i + offset) or \
           binary_search.present(exclude_sites2, i + offset):
            info_mask[i] = True
            continue
        if seq1[i] not in skip and seq2[i] not in skip:
            total_sites += 1
            if seq1[i] == seq2[i]:
                total_match += 1

    # TODO: keep track of gapped/masked sites for master/predicted to
    # incorporate into info string later
    return total_match, total_sites, {'mask_flag': info_mask}


def make_info_string_unknown(info, master_ind):

    # used with indices to decode result
    decoder = np.array(list('Xx._-'))
    indices = np.zeros(info['gap_any_flag'].shape, int)

    indices[info['match_flag'][:, master_ind]] = 1  # x
    matches = np.all(info['match_flag'], axis=1)
    indices[matches] = 2  # .
    indices[info['mask_any_flag']] = 3  # _
    indices[info['gap_any_flag']] = 4  # -

    return ''.join(decoder[indices])


def make_info_string(info, master_ind, predict_ind):

    if predict_ind >= info['match_flag'].shape[1]:
        return make_info_string_unknown(info, master_ind)

    decoder = np.array(list('xXpPcCbB._-'))
    indices = np.zeros(info['match_flag'].shape[0], int)

    indices[info['match_flag'][:, predict_ind]] += 2  # x to p if true
    indices[info['match_flag'][:, master_ind]] += 4  # x to c, p to b
    indices[info['hmm_flag']] += 1  # to upper

    matches = np.all(info['match_flag'], axis=1)
    indices[matches] = 8  # .
    indices[np.any(
        info['mask_flag'][:, [master_ind, predict_ind]],
        axis=1)] = 9  # _
    indices[np.any(
        info['gap_flag'][:, [master_ind, predict_ind]],
        axis=1)] = 10  # -

    return ''.join(decoder[indices])

def read_region_file(fn):
    f = gzip.open(fn, 'rb')
    d = {}
    line = f.readline().decode()
    while line != '':
        region_id = line[1:-1]
        line = f.readline().decode()
        seqs = {}
        while line[0] != '#':
            line = line[:-1].split(' ')
            strain = line[1]
            seqs[strain] = {}
            if len(line) > 2:
                seqs[strain]['start'] = int(line[2])
                seqs[strain]['end'] = int(line[3])
            seqs[strain]['seq'] = f.readline().decode()[:-1]
            line = f.readline().decode()
            if line == '':
                break
        d[region_id] = seqs

    f.close()
    return d
    
