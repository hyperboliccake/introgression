
r = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'a':'t', 't':'a', 'g':'c', 'c':'g'}

def reverse_complement(s):
    rs = ''
    for b in s[::-1]:
        if b in r:
            rs += r[b]
        else:
            rs += b
    return rs


def index_ignoring_gaps(s, i, s_start, gap_symbol):
    '''returns the index of the ith (starting at 0) non-gap character in
    s, given that the start index of s is s_start (instead of needing
    to be 0); for example, s='-A-AA-A', i=2, s_start=0 => returns 4
    '''

    assert s_start >= 0

    x = 0
    non_gap_count = 0
    i -= s_start
    if i < 0:
        return -1
    while x < len(s):
        if s[x] != gap_symbol and non_gap_count >= i:
            return x
        if s[x] != gap_symbol:
            non_gap_count += 1
        x += 1
    return x

def seq_id(ref_seq, seq):
    n = len(ref_seq)
    total_sites = 0
    total_match = 0
    for i in range(n):
        if ref_seq[i] in r and seq[i] in r:
            total_sites += 1
            if ref_seq[i] == seq[i]:
                total_match += 1
    if total_sites == 0:
        return 'NA'
    return float(total_match) / total_sites