import os
import math
import numpy.random

def parse_list(l):
    assert len(l) >= 2, l
    if len(l) == 2:
        return []
    return [float(x) for x in l[1:-1].split(',')]
    
def mean(l):
    if len(l) == 0:
        return 'NA'
    return float(sum(l)) / len(l)

def std_dev(l):
    if len(l) == 0:
        return 'NA'
    if len(l) == 1:
        return 0
    m = mean(l)
    return math.sqrt(sum([(x - m)**2 for x in l]) / (len(l) - 1))

def std_err(l):
    if len(l) == 0:
        return 'NA'
    return std_dev(l) / math.sqrt(len(l))

def bootstrap(l, n = 100, alpha = .05):
    x = len(l)
    if x == 0:
        return 'NA', 'NA'
    a = []
    for i in range(n):
        a.append(mean(numpy.random.choice(l, size = x, replace = True)))
    a.sort()
    print len(a), a.count(0)
    print mean(a)
    return a[int(alpha * n * .5)], a[int((1 - alpha * .5) * n)]


params = [line.strip().split(' ') for line in open('sim_multi_model_args.txt', 'r').readlines()]
param_names = ['tag', 'model', 'N0', 'num_samples_par', 'num_samples_cer', 'par_cer_migration', 't_par_cer', 'num_sites', 'rho', 'outcross_rate', 'num_reps']
assert len(param_names) == len(params[0]), str(len(param_names)) + ' != ' + str(len(params[0]))
out_dir = '../../results/sim/run_3/'
prefix = 'sim_out_'
suffix = '_summary.txt'
ids = range(1, len(params) + 1)
f_all = open(out_dir + prefix + 'all' + suffix, 'w')
for i in ids:
    print i
    f = open(out_dir + prefix + str(i) + suffix, 'r')
    col_names = f.readline().strip().split('\t') # header
    if i == ids[0]:
        for p in param_names:
            f_all.write(p + '\t')
        f_all.write('row_type\t')
        for c in col_names[:-1]:
            f_all.write(c + '\t')
        f_all.write(col_names[-1] + '\n')
    line = f.readline()
    agg = [[] for item in range(len(col_names))]
    while line != '':
        line = line.strip().split('\t')
        for x in range(len(line)):
            item = line[x]
            if item[0] == '[':
                item = parse_list(item)
                agg[x] += item
            else:
                agg[x].append(float(item))
        line = f.readline()
    # mean row
    for p in params[int(i) - 1]:
        f_all.write(p + '\t')
    f_all.write('mean')
    for item in agg:
        m = mean(item)
        f_all.write('\t' + str(m))
    f_all.write('\n')
    # standard error row
    for p in params[int(i) - 1]:
        f_all.write(p + '\t')
    f_all.write('std_err')
    for item in agg:
        se = std_err(item)
        f_all.write('\t' + str(se))
    f_all.write('\n')
    # bootstrap
    bls = []
    bus = []
    for item in agg:
        bl, bu = bootstrap(item)
        print bl, bu
        bls.append(bl)
        bus.append(bu)
    print bls
    print bus
    # lower bootstrap row
    for p in params[int(i) - 1]:
        f_all.write(p + '\t')
    f_all.write('lower')
    for x in range(len(agg)):
        f_all.write('\t' + str(bls[x]))
    f_all.write('\n')
    # upper bootstrap row
    for p in params[int(i) - 1]:
        f_all.write(p + '\t')
    f_all.write('upper')
    for x in range(len(agg)):
        f_all.write('\t' + str(bus[x]))
    f_all.write('\n')
    f.close()
f_all.close()
