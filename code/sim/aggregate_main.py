import sys
from aggregate import *
import process_args
sys.path.append('..')
import global_params as gp

# read in all summary files and generate one summary file with
# averages, std devs for all stats

arg_lines = open(sys.argv[1], 'r').readlines()
tags = [l.split(' ')[0] for l in arg_lines]

gp_dir = '../'
fns = [gp_dir + gp.sim_out_dir + gp.sim_out_prefix + tag + '_' + sys.argv[2] + '.txt' \
       for tag in tags]

fn_all = gp_dir + gp.sim_out_dir + gp.sim_out_prefix + 'all_' + sys.argv[2] + '.txt'

aggregate_summary_files(fns, fn_all, tags)
