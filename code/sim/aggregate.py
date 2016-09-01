import os

params = [line.split(' ') for line in open('sim_multi_model_args.txt', 'r').readlines()]
param_names = ['tag', 'model', 'num_samples_par', 'num_samples_cer', 'par_cer_migration', 't_par_cer', 'num_sites', 'rho', 'outcross_rate', 'mu', 'theta', 'num_reps']
assert len(param_names) == len(params]0])
out_dir = '../../results/sim/run3/'
prefix = 'sim_out_'
suffix = '_summary.txt'
ids = range(1, len(args) + 1)
f_all = open(out_dir + prefix + 'all' + suffix, 'w')
types = ['list'] * 20 + 
for i in ids:
    f = open(out_dir + prefix + str(i) + suffix, 'r')
    col_names = f.readline() # header
    if i == ids[0]:
        for p in param_names:
            f_all.write(p + '\t')
        f_all.write('row_type\t')
        for c in col_names[:-1]:
            f_all.write(c + '\t')
        f_all.write(col_names[-1] + '\n')
    line = f.readline()
    means = []
    sds = []
    while line != '':
        line = line.split('\t')
        for item in line:
            if item[0] == '[':
                item = parse_list(item)
            m = mean(item)
            sd = std_dev(item)
        line = f.readline()
    f.close()
f_all.close()
