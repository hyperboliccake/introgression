from sim_analyze_hmm_bw import *
from concordance_functions import *
sys.path.insert(0, '..')
import global_params as gp

def read_list(l):
    return [x.strip() for x in l[1:-1].split(',')]

gp_dir = '../'
outfilename = gp.sim_out_prefix + tag + '.txt'
results_filename = gp.sim_out_prefix + tag + '_summary.txt'
