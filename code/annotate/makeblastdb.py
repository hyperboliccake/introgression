import sys
import os
#from orf import *
sys.path.insert(0, '../align')
import align_helpers
sys.path.insert(0, '..')
import global_params as gp

#d = '/tigress/AKEY/akey_vol2/aclark4/nobackup/100_genomes/genomes_gb/orfs/'
d = '../../data/CBS432/orfs/'
#d = '/tigress/AKEY/akey_vol2/aclark4/nobackup/100_genomes/genomes_gb/orfs/'
fns = os.listdir(d)
for fn in fns:
    cmd_string = gp.blast_install_path + 'makeblastdb' + \
                 ' -dbtype nucl' + \
                 ' -in ' + d + fn
    print cmd_string
    os.system(cmd_string)
