import sys
from sim_analyze_hmm_bw import *
from concordance_functions import *
sys.path.insert(0, '..')
import global_params as gp

def process_list(l):
    r = [x.strip() for x in l[1:-1].split(',')]
    return [process_entry(x) for x in r]

def process_entry(l):
    if len(l) == 0:
        return None
    if l[0] == '[':
        return process_list(l)
    try:
        return float(l)
    except:
        return l

def process_line(l, d, labels):
    l = l.strip().split('\t')
    for i in range(len(l)):
        r = process_entry(l[i])
        d[labels[i]].append(r)

tag = sys.argv[1]

gp_dir = '../'
results_filename = gp_dir + gp.sim_out_dir + gp.sim_out_prefix + tag + '_summary.txt'
processed_filename = gp_dir + gp.sim_out_dir + gp.sim_out_prefix + tag + '_summary_processed.txt'

#####
# read results file into a table
#####

f = open(results_filename, 'r')

# header
line = f.readline()
labels = line.split('\t')

results_table = dict(zip(labels, [[] for x in range(len(labels))]))
num_sims = 0

line = f.readline()
while line != '':
    process_line(line, results_table, labels)
    num_sims += 1
    line = f.readline()
f.close()

#####
# do things with it...and write them to summary file in ggplot form
# (one row per thing)
#####

f = open(processed_filename, 'w')
f.write('variable value\n')

# power (to detect individual introgressed bases)

# num_bases_actual_par_predicted_par / num_bases_actual_par_predicted_*

# average power across all simulations

power = []
for i in range(num_sims):
    num = sum(results_table['num_bases_actual_par_predicted_par'][i])
    den = sum(results_table['num_bases_actual_par_predicted_par'][i])
    den += sum(results_table['num_bases_actual_par_predicted_cer'][i])
    if den == 0:
        power.append('NA')
    else:
        power.append(num / den)

for x in power:
    f.write('power ' + str(x) + '\n')

f.close()
