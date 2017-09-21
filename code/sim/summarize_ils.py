import os
import math
import numpy.random
import sys
sys.path.insert(0, '../misc/')
import mystats

params = [line.strip().split(' ') for line in open('sim_multi_model_args.txt', 'r').readlines()]
param_names = ['tag', 'model', 'N0', 'num_samples_par', 'num_samples_cer', 'par_cer_migration', 't_par_cer', 'num_sites', 'rho', 'outcross_rate', 'num_reps']
assert len(param_names) == len(params[0]), str(len(param_names)) + ' != ' + str(len(params[0]))
out_dir = '../../results/sim/run_3/'
prefix = 'sim_out_'
suffix = '_summary.txt'
ids = range(1, len(params) + 1)
f_out = open(out_dir + prefix + 'ils' + suffix, 'w')
for i in ids:
    print i
    f = open(out_dir + prefix + str(i) + suffix, 'r')
    col_names = f.readline().strip().split('\t') # header
    print col_names[22]
    if i == ids[0]:
        for p in param_names[:-1]:
            f_out.write(p + '\t')
        f_out.write(p[-1])
        for x in range(100):
            f_out.write('\t' + col_names[22] + '.' + str(x))
        f_out.write('\n')
    
    agg = []
    line = f.readline()
    while line != '':
        line = line.strip().split('\t')
        agg.append(line[22])
        line = f.readline()
    # mean row
    for p in params[int(i) - 1][:-1]:
        f_out.write(p + '\t')
    f_out.write(params[int(i) - 1][-1])
    for item in agg:
        f_out.write('\t' + str(item))
    f_out.write('\n')
f_out.close()
