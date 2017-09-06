import os
from alignment_stats import *
from align_helpers import *
sys.path.insert(0, '..')
import global_params as gp

# gives info related to how good an alignment is:
# - number of sites where 3, 2, 1, genomes aligned
# - percent of each genome contained in alignment (hopefully 100%)
# - using each genome as reference, percentage of other genomes aligned

# get all non-reference strains of cerevisiae and paradoxus
s = get_strains(flatten(gp.non_ref_dirs.values()))
strain, d = s[int(sys.argv[1])]
gp_dir = '../'

fn_start = gp_dir + gp.alignments_dir + '_'.join(gp.alignment_ref_order) + '_' + strain + '_chr'

for chrm in gp.chrms:
    print chrm
    sys.stdout.flush()

    headers, seqs = read_fasta.read_fasta(fn_start + chrm + '_mafft.maf')
    a = dict(zip(headers, seqs))
    
    f_out = open(fn_start + chrm + '_mafft.stats', 'w')

    # number of sites where 3,2,1 genomes aligned
    num_strains_by_site = num_strains_aligned_by_site(seqs)
    f_out.write(\
        '# histogram of number of strains aligned across all alignment columns\n')
    for n in range(len(num_strains_by_site)):
        f_out.write(str(n) + ',' + str(num_strains_by_site[n]) + '\n')
    f_out.write('\n')

    # fraction of genomes aligned (should all be 1)
    fracs_aligned, seq_lengths = fraction_strains_aligned(headers, seqs)
    for frac in fracs_aligned:
        assert frac == 1, fracs_aligned

    # length of chromosomes
    f_out.write('chromosome aligned lengths\n')
    for n in range(len(seqs)):
        f_out.write(headers[n].split(' ')[0][1:] + ',' + str(seq_lengths[n]) + '\n')
    f_out.write('\n')

    # using each genome as reference, fraction of other genomes aligned
    f_out.write('fraction aligned to reference\n')
    frac_aligned_to_ref = frac_aligned_to_reference(seqs, seq_lengths)
    for ref in range(len(seqs)):
        f_out.write(headers[ref].split(' ')[0][1:])
        for other in range(len(seqs)): 
            f_out.write(',' + str(frac_aligned_to_ref[ref][other]))
        f_out.write('\n')
    f_out.write('\n')

    f_out.close()
