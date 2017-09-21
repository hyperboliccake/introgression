import sys
import os
from align_helpers import *
sys.path.insert(0, '..')
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

strain, d = s[int(sys.argv[1])]

print strain

# building up one command string so that we don't create a new
# shell instance every time (I think there's a limit on the
# command character count or something which is why we're not
# making a single string for all strains)
#cmd_string = ''

current_strain_fn = d + strain_fn.replace('*', strain)

for chrm in gp.chrms:
    print chrm
    sys.stdout.flush()

    align_fn = ref_prefix + strain + '_chr' + chrm + \
        '_mafft' + gp.alignment_suffix
    align_fn_abs = gp_dir + gp.alignments_dir + align_fn
    # if we don't already have an alignment for this strain/chromosome
    # (or that alignment file is empty), then make one
    if (align_fn not in a) or (os.stat(align_fn_abs).st_size == 0):
        cmd_string = ''

        # first put all sequences in same (temporary) file
        ref_fns_chrm = [x.replace('?', chrm) for x in ref_fns]
        current_strain_fn_chrm = current_strain_fn.replace('?', chrm)
        combined_fn = 'run_mafft_' + strain + chrm + '.temp'
        
        concatenate_fasta(ref_fns_chrm + [current_strain_fn_chrm], combined_fn)
        
        cmd_string += gp.mafft_install_path + '/mafft ' + \
            combined_fn + ' > ' + align_fn_abs + '; '
        
        cmd_string += 'rm ' + combined_fn + ';'

        print cmd_string
        sys.stdout.flush()
        os.system(cmd_string)

        # want some kind of indication if alignment fails (due to
        # running out of memory probably)
        if os.stat(align_fn_abs).st_size == 0:
            print 'alignment failed: ' + strain + ' chromosome ' + chrm
            sys.stdout.flush()
            sys.exit()
    else:
        print "already did this alignment: " + strain + ' chromosome ' + chrm

#print cmd_string
#sys.stdout.flush()
#os.system(cmd_string)
