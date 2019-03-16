import sys
import os
from mask_helpers import *
import align_helpers
from analyze import read_args
import global_params as gp

args = read_args.read_setup_args(sys.argv[1])

only_ref = True

s = []

if not only_ref:
    # get all non-reference strains of cerevisiae and paradoxus
    s = args['strain_dirs']

a = []

for r in args['references']:
    s.append((r, args['reference_directories'][r]))

strain_fn = '*_chr?' + gp.fasta_suffix
strain_masked_fn = '*_chr?_masked' + gp.fasta_suffix
intervals_fn = '*_chr?_intervals.txt'
intervals_d = gp.mask_dir

for i in range(len(s)):
    strain, d = s[i]

    print(strain)
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
        



