import sys
import os
from align_helpers import *
from analyze import read_args
import global_params as gp

args = read_args.read_setup_args(sys.argv[2])

masked = False
mask_suffix = ''
if masked:
    mask_suffix = '_masked'

ref_only = False

a = []
if gp.resume_alignment:
    for fn in os.listdir(args['alignments_directory']):
        if os.stat(args['alignments_directory'] + fn).st_size != 0:
            a.append(fn)
ref_prefix = '_'.join(args['references']) + '_'
ref_fns = [args['reference_directories'][r] + r + '_chr' + '?' + \
           mask_suffix + gp.fasta_suffix \
           for r in args['references']]

if ref_only:

    chrm = gp.chrms[int(sys.argv[1])]

    ref_fns_chrm = [x.replace('?', chrm) for x in ref_fns]
    combined_fn = 'run_mafft_' + chrm + '.temp'

    concatenate_fasta(ref_fns_chrm, \
                      args['references'], combined_fn)
    
    align_fn = ref_prefix + 'chr' + chrm + \
        '_mafft' + gp.alignment_suffix
    align_fn_abs = args['alignments_directory'] + align_fn

    cmd_string = gp.mafft_install_path + '/mafft ' + \
                 combined_fn + ' > ' + align_fn_abs + '; '
        
    cmd_string += 'rm ' + combined_fn + ';'

    print(cmd_string)
    os.system(cmd_string)
    sys.stdout.flush()

    sys.exit()

# get all non-reference strains of cerevisiae and paradoxus
s = get_strains(args['strain_dirs'])

strain_fn = '*_chr?' + mask_suffix + gp.fasta_suffix

strain, d = s[int(sys.argv[1])]

print(strain)

# building up one command string so that we don't create a new
# shell instance every time (I think there's a limit on the
# command character count or something which is why we're not
# making a single string for all strains)
#cmd_string = ''

current_strain_fn = d + strain_fn.replace('*', strain)

for chrm in gp.chrms:
    print(chrm)
    sys.stdout.flush()

    align_fn = ref_prefix + strain + '_chr' + chrm + \
        '_mafft' + gp.alignment_suffix
    align_fn_abs = args['alignments_directory'] + align_fn
    # if we don't already have an alignment for this strain/chromosome
    # (or that alignment file is empty), then make one
    #if (align_fn not in a) or (os.stat(align_fn_abs).st_size == 0):
    if align_fn not in a:
        cmd_string = ''

        # first put all sequences in same (temporary) file
        ref_fns_chrm = [x.replace('?', chrm) for x in ref_fns]
        current_strain_fn_chrm = current_strain_fn.replace('?', chrm)
        combined_fn = 'run_mafft_' + strain + chrm + '.temp'
        
        concatenate_fasta(ref_fns_chrm + [current_strain_fn_chrm], \
                          args['references'] + [strain], combined_fn)
        
        # add --ep 0.123 to maybe get shorter alignment
        #cmd_string += gp.mafft_install_path + '/mafft --ep 0.123 ' + \
        #    combined_fn + ' > ' + align_fn_abs + '; '
        #cmd_string += gp.mafft_install_path + '/mafft --retree 1 ' + \
        #    combined_fn + ' > ' + align_fn_abs + '; '
        cmd_string += gp.mafft_install_path + '/mafft ' + \
            combined_fn + ' > ' + align_fn_abs + '; '
        
        cmd_string += 'rm ' + combined_fn + ';'

        print(cmd_string)
        sys.stdout.flush()
        os.system(cmd_string)

        # want some kind of indication if alignment fails (due to
        # running out of memory probably)
        if os.stat(align_fn_abs).st_size == 0:
            print('alignment failed: ' + strain + ' chromosome ' + chrm)
            sys.stdout.flush()
            sys.exit()
    else:
        print("already did this alignment: " + strain + ' chromosome ' + chrm)
        sys.stdout.flush()

#print cmd_string
#sys.stdout.flush()
#os.system(cmd_string)
