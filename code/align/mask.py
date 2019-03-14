import sys
import os
from mask_helpers import *
import align_helpers
sys.path.insert(0, '..')
import global_params as gp

only_ref = True

s = []

if not only_ref:
    # get all non-reference strains of cerevisiae and paradoxus
    s = align_helpers.get_strains(align_helpers.flatten(gp.non_ref_dirs.values()))

gp_dir = '../'
a = []

for r in gp.alignment_ref_order:
    s.append((gp.ref_fn_prefix[r], gp.ref_dir[r]))

strain_fn = '*_chr?' + gp.fasta_suffix
strain_masked_fn = '*_chr?_masked' + gp.fasta_suffix
intervals_fn = '*_chr?_intervals.txt'
intervals_d = gp_dir + gp.mask_dir

for i in range(len(s)):
    strain, d = s[i]

    print strain
    sys.stdout.flush()

    current_strain_fn = d + strain_fn.replace('*', strain)
    current_strain_masked_fn = d + strain_masked_fn.replace('*', strain)
    current_strain_intervals_fn = intervals_d + intervals_fn.replace('*', strain)

    for chrm in gp.chrms:

        in_fn = current_strain_fn.replace('?', chrm)
        out_fn = current_strain_intervals_fn.replace('?', chrm)
    
        # get dustmasker intervals
        cmd_string = gp.blast_install_path + 'dustmasker' + \
                     ' -in ' + in_fn + \
                     ' -out ' + out_fn + \
                     ' -outfmt interval'
        
        os.system(cmd_string)

        # replace those intervals with Ns and write to masked fasta file
        masked_fn = current_strain_masked_fn.replace('?', chrm)
        mask(in_fn, masked_fn, out_fn)
        



