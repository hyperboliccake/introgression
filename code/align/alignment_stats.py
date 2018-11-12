import sys
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc')
import read_fasta

# count sites where n, ..., 3, 2, 1 genomes aligned, etc.
def num_strains_aligned_by_site(seqs):
    nseqs = len(seqs)
    nsites = len(seqs[0])
    # number of sites with 0-n strains aligned
    num_strains_hist = [0] * (nseqs + 1)
    for i in range(nsites):
        column = [seq[i] for seq in seqs]
        num_aligned = nseqs - column.count(gp.gap_symbol)
        num_strains_hist[num_aligned] += 1

    return num_strains_hist

# fraction of each strain's sequence contained in alignment 
# (should be 1)
def fraction_strains_aligned(headers, seqs):
    nseqs = len(seqs)
    nsites = len(seqs[0])
    seq_lengths = []
    fracs_aligned = []
    for i in range(nseqs):
        h = headers[i].split(' ')
        actual = nsites - seqs[i].count(gp.gap_symbol)
        seq_lengths.append(actual)
        s = read_fasta.read_fasta(h[-1])
        expected = len(s[1][0])
        fracs_aligned.append(float(actual)/expected)

    return fracs_aligned, seq_lengths

# using each genome as reference, percentage of other genomes aligned
def frac_aligned_to_reference(seqs, seq_lengths):
    nseqs = len(seqs)
    nsites = len(seqs[0])
    fracs_aligned_to_ref = []
    for ref in range(nseqs):
        r = []
        for other in range(nseqs):
            if other == ref:
                r.append(1)
            else:
                total = 0
                for i in range(nsites):
                    if seqs[ref][i] != gp.gap_symbol and seqs[other][i] != gp.gap_symbol:
                        total += 1
                r.append(float(total) / seq_lengths[other])
        fracs_aligned_to_ref.append(r)
    return fracs_aligned_to_ref
