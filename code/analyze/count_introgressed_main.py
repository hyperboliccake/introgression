# counts total amount of sites introgressed on each chromosome

import re
import sys
import os
import math
import gzip
import itertools
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../align/')
import align_helpers
sys.path.insert(0, '../misc/')
import overlap
import read_table
import read_fasta
import write_fasta
import mystats


chrm_sizes = [230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751, 666816, 1078177, 924431, 784333, 1091291, 948066]

tag = 'u3_i.001_tv_l1000_f.01'

regions_by_chrm = dict(zip(gp.chrms, [[] for i in range(len(gp.chrms))]))
fn_regions = gp.analysis_out_dir_absolute + tag + '/' + \
             'introgressed_blocks_filtered_par_' + tag + '_summary_plus.txt'
d, labels = read_table.read_table_rows(fn_regions, '\t')

for region in d:
    chrm = d[region]['chromosome']
    strain = d[region]['strain']
    regions_by_chrm[chrm].append((strain, \
                                  int(d[region]['start']), \
                                  int(d[region]['end'])))

hist = {}
for chrm in gp.chrms:
    print chrm
    chrm_size = chrm_sizes[gp.chrms.index(chrm)]
    x = [0 for i in range(chrm_size)]
    for ri in range(len(regions_by_chrm[chrm])):
        strain, start, end = regions_by_chrm[chrm][ri]
        for i in range(start, end+1):
            x[i] += 1
    max_count = max(x)
    hist_chrm = []
    for count in range(max_count + 1):
        hist_chrm.append(x.count(count))
    hist[chrm] = hist_chrm

fn = gp.analysis_out_dir_absolute + tag + '/' + \
     'count_introgressed.txt'
f = open(fn, 'w')
f.write('chromosome\tat_least_one\tat_least_one_frac\tchromosome_size\thist\n')
total = 0
for chrm in gp.chrms:
    f.write(chrm + '\t')
    chrm_size = chrm_sizes[gp.chrms.index(chrm)]    
    at_least_one = chrm_size - hist[chrm][0]
    total += at_least_one
    f.write(str(at_least_one) + '\t')
    f.write(str(float(at_least_one)/chrm_size) + '\t')
    f.write(str(chrm_size) + '\t')
    f.write(','.join([str(c) for c in hist[chrm]]) + '\n')
f.write('all\t' + str(total) + '\t')
f.write(str(float(total)/sum(chrm_sizes)) + '\t')
f.write(str(sum(chrm_sizes)) + '\t')
f.write('\n')
f.close()
