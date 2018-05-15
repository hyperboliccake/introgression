import sys
import os
#from orf import *
sys.path.insert(0, '../align')
import align_helpers
sys.path.insert(0, '..')
import global_params as gp

#d = '/tigress/AKEY/akey_vol2/aclark4/nobackup/100_genomes/genomes_gb/orfs/'
d = '../../data/CBS432/orfs/'
fns = os.listdir(d)
for fn in fns:
    new_fn = fn[len('S288c_CBS432_'):]
    os.system('mv ' + d + fn + ' ' + d + new_fn)
