# lol because i'm so bad at R

import sys
sys.path.insert(0, '..')
import global_params as gp

tag = 'u3_i.001_tv_l1000_f.01'
fn = gp.analysis_out_dir_absolute + tag + '/' + \
     '/polymorphism/polymorphism.txt'
f = open(fn, 'r')
lines = [line[:-1].split('\t') for line in f.readlines()]
f.close()

d = {}
d_sums = {}
d2_sums = {}
for line in lines[1:]:
    chrm = line[0]
    if not d_sums.has_key(chrm):
        d_sums[chrm] = 0
        d2_sums[chrm] = 0
        d[chrm] = {}
    c = int(line[-1])
    if line[1] != '1':
        d_sums[chrm] += c
    if line[1] == '2':
        d2_sums[chrm] += c
    d[chrm][tuple(line[1:5])] = c


fn = gp.analysis_out_dir_absolute + tag + '/' + \
     '/polymorphism/polymorphism_summary.txt'
f = open(fn, 'w')
f.write('chromosome\tmatch\tsites\tfrac\tcount\n')
for chrm in gp.chrms + ['all']:
    # for frac one, want * false true *
    # for frac one biallelic, want 2 false true *
    # for frac all, want * false true true
    # for frac all biallelic, want 2 false true true
    fo = 0
    fob = 0
    fa = 0
    fab = 0
    for key in d[chrm].keys():
        c = d[chrm][key]
        if key[1] == 'False' and key[2] == 'True' and key[0] != '1':
            fo += c
            if key[3] == 'True':
                fa += c
            if key[0] == '2':
                fob += c
                if key[3] == 'True':
                    fab += c
    try:
        fo = str(float(fo)/d_sums[chrm])
    except:
        fo = 'NaN'
    try:
        fob = str(float(fob)/d2_sums[chrm])
    except:
        fob = 'NaN'
    try:
        fa = str(float(fa)/d_sums[chrm])
    except:
        fa = 'NaN'
    try:
        fab = str(float(fab)/d2_sums[chrm])
    except:
        fab = 'NaN'

    f.write(chrm + '\tone\tpolymorphic\t' + fo + '\t' + str(d_sums[chrm]) + '\n')
    f.write(chrm + '\tone\tbiallelic\t' + fob + '\t' + str(d2_sums[chrm]) + '\n')
    f.write(chrm + '\tall\tpolymorphic\t' + fa + '\t' + str(d_sums[chrm]) + '\n')
    f.write(chrm + '\tall\tbiallelic\t' + fab + '\t' + str(d2_sums[chrm]) + '\n')

f.close()
