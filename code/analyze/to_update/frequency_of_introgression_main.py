import re
import sys
import os
import copy
import itertools
from collections import defaultdict
import gene_predictions
import predict
from filter_helpers import *
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../sim/')
import sim_analyze_hmm_bw as sim
sys.path.insert(0, '../misc/')
import mystats
import read_table
import read_fasta

tag = 'u3_i.001_tv_l1000_f.01'
species_from = 'par'

#strains3 = ['yjm1252', 'yjm1078', 'yjm248']

fn = gp.analysis_out_dir_absolute + tag + '/' + \
     'introgressed_blocks_filtered_' + species_from + \
     '_' + tag + '_summary_plus.txt'
regions, fields = read_table.read_table_rows(fn, '\t')

bases_by_strains = defaultdict(lambda: defaultdict(int))
strains = set([])
for region_id in regions:
    strain = regions[region_id]['strain']
    #if strain not in strains3:
    strains.add(strain)
    chrm = regions[region_id]['chromosome']
    start = int(regions[region_id]['start'])
    end = int(regions[region_id]['end'])
    for base in range(start, end + 1):
        bases_by_strains[chrm][base] += 1

counts = defaultdict(int)
for chrm in bases_by_strains.keys():
    for base in bases_by_strains[chrm]:
        counts[bases_by_strains[chrm][base]] += 1

fn = gp.analysis_out_dir_absolute + tag + '/' + \
     'introgressed_frequency.txt'
f = open(fn, 'w')
f.write('num_strains\tcount\n')
for i in range(len(strains)):
    f.write(str(i) + '\t' + str(counts[i]) + '\n')
f.close()

