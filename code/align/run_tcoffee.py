import sys
import os
sys.path.insert(0, '..')
from align_helpers import *
import global_params as gp

# get all non-reference strains of cerevisiae and paradoxus
s = get_strains(flatten(gp.non_ref_dirs.values()))

gp_dir = '../'
a = []
if gp.resume_alignment:
    a = os.listdir(gp_dir + gp.alignments_dir)


ref_prefix = '_'.join(gp.alignment_ref_order) + '_'
ref_fns = [gp.ref_dir[r] + gp.ref_fn_prefix[r] + '_chr' + '?' + gp.fasta_suffix \
               for r in gp.alignment_ref_order]

strain_fn = '*_chr?' + gp.fasta_suffix

for strain, d in s:
    print strain

    # building up one command string so that we don't create a new
    # shell instance every time (I think there's a limit on the
    # command character count or something which is why we're not
    # making a single string for all strains)
    cmd_string = ''

    current_strain_fn = d + strain_fn.replace('*', strain)
        
    for chrm in gp.chrms[:2]:
        print chrm
        align_fn = ref_prefix + strain + '_chr' + chrm + \
            '_tcoffee' + gp.alignment_suffix
        # if we don't already have an alignment for this
        # strain/chromosome, then make one
        if align_fn not in a:
            # first put all sequences in same (temporary) file
            ref_fns_chrm = [x.replace('?', chrm) for x in ref_fns]
            current_strain_fn_chrm = current_strain_fn.replace('?', chrm)
            combined_fn = 'run_tcoffee_' + strain + chrm + '.temp'

            concatenate_fasta(ref_fns_chrm + [current_strain_fn_chrm], combined_fn)

            cmd_string += gp.tcoffee_install_path + '/t_coffee ' + \
                combined_fn + '; '

            #cmd_string += 'rm ' + combined_fn + ';'

    # commands can only be up to a certain length so break it up this way
    print cmd_string
    os.system(cmd_string)
    sys.exit()
