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

# need to add this on the start of each command because os.system()
# creates a new shell instance every time
cmd_string_start = 'export MUGSY_INSTALL=' + gp.mugsy_install_path + '; '
cmd_string_start += 'export PATH=$PATH:$MUGSY_INSTALL:$MUGSY_INSTALL/mapping; '
cmd_string_start += 'export PERL5LIB=$MUGSY_INSTALL/perllibs; '

ref_prefix = '_'.join(gp.alignment_ref_order) + '_'
ref_dirs = [gp.ref_dir[ref] for ref in gp.alignment_ref_order]

for strain, d in s:
    print strain

    cmd_string = cmd_string_start
        
    for chrm in [gp.chrms[-1]]:
        align_fn = ref_prefix + strain + '_chr' + chrm + gp.alignment_suffix
        # if we don't already have an alignment for this strain/chromosome, then make one
        if align_fn not in a:
            cmd_string += gp.mugsy_install_path + '/mugsy ' + \
                '--directory ' + gp_dir + gp.alignments_dir + ' ' + \
                '--prefix ' + ref_prefix + strain + '_chr' + chrm
            for ref in gp.alignment_ref_order:
                cmd_string += ' ' + gp.ref_dir[ref] + '/' + \
                gp.ref_fn_prefix[ref] + '_chr' + chrm + gp.fasta_suffix
            cmd_string += ' ' + d + '/' + strain + '_chr' + chrm + gp.fasta_suffix + '; '


    # commands can only be up to a certain length so break it up this way
    print cmd_string
    os.system(cmd_string)
