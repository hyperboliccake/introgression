# extract open reading frame sequences from all strains x chromosomes
# do this using orffinder, and take all non-overlapping above a certain length

import sys
import os
#from orf import *
sys.path.insert(0, '../align')
import align_helpers
sys.path.insert(0, '..')
import global_params as gp

ref_fns = [gp.ref_dir[r] + gp.ref_fn_prefix[r] + '_chr' + '?' + \
           gp.fasta_suffix \
           for r in gp.alignment_ref_order]

# get all non-reference strains of cerevisiae and paradoxus
s = align_helpers.get_strains(align_helpers.flatten(gp.non_ref_dirs.values()))
# and get paradoxus reference as well
s.append((gp.ref_fn_prefix[gp.alignment_ref_order[1]], gp.ref_dir[gp.alignment_ref_order[1]]))

strain_fn = '*_chr?' + gp.fasta_suffix

f = open('orfs.sh', 'w')
for i in range(78,94):

    strain, d = s[i]
    print strain
    current_strain_fn = strain_fn.replace('*', strain)
    for chrm in gp.chrms:
        print chrm
        sys.stdout.flush()
        
        current_strain_chrm_fn = current_strain_fn.replace('?', chrm)
        orf_fn = strain + '_chr' + chrm + \
                 '_orfs' + gp.fasta_suffix
        orf_d = d + '/orfs/'
        if not os.path.isdir(orf_d):
            os.makedirs(orf_d)

        cmd_string = gp.orffinder_install_path + '/ORFfinder' + \
                     ' -in ' + d + current_strain_chrm_fn + \
                     ' -s 0' + \
                     ' -out ' + orf_d + orf_fn + \
                     ' -outfmt 1 -n true; \n'
        #print cmd_string
        os.system(cmd_string)
        f.write(cmd_string)
f.close()


# "../../../../software/ORFfinder -in /tigress/AKEY/akey_vol2/aclark4/nobackup/100_genomes/genomes_gb/yjm248_chrI.fa -out a.txt -outfmt 1 -n true"
