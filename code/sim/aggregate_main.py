import sys
from aggregate import *
import process_args
sys.path.append('..')
import global_params as gp

# read in all summary files and generate one summary file with
# averages, std devs for all stats

arg_lines = open(sys.argv[1], 'r').readlines()
tags = [l.split(' ')[0] for l in arg_lines]
# summary, introgressed_actual, introgressed_predicted, etc
output_type  = sys.argv[2] 
sim_id = sys.argv[3]

gp_dir = '../'
fns = [gp_dir + gp.sim_out_dir + gp.sim_out_prefix + tag + '_' + output_type + '.txt' \
       for tag in tags]

fn_all = gp_dir + gp.sim_out_dir + gp.sim_out_prefix + 'all_' + output_type + \
         '_' + sim_id + '.txt'

aggregate_summary_files(fns, fn_all, tags)
