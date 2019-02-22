
r = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'a':'t', 't':'a', 'g':'c', 'c':'g'}

codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

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
    return total_match, total_sites

def seq_id_windowed(seq1, seq2, window):
    n = len(seq1)
    total_sites = 0
    total_match = 0
    matches = []
    for i in range(n):
        if seq1[i] in r and seq2[i] in r:
            total_sites += 1
            if seq1[i] == seq2[i]:
                total_match += 1
            if total_sites == window:
                matches.append(float(total_match)/window)
                total_sites = 0
                total_match = 0
    return matches

def translate(seq):
    if len(seq) % 3 != 0:
        return ''
    a = ''
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        a += codon_table[codon]
    return a

