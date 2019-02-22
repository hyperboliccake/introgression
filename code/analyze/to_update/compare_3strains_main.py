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

strains = ['yjm1252', 'yjm1078', 'yjm248']

tag = 'u3_i.001_tv_l1000_f.01'
species_from = 'par'

fn = gp.analysis_out_dir_absolute + tag + '/' + \
     'introgressed_blocks_filtered_' + species_from + \
     '_' + tag + '_summary_plus' + '.txt'
regions, fields = read_table.read_table_rows(fn, '\t')

bases_by_strains = defaultdict(lambda: defaultdict(list))
for region_id in regions:
    strain = regions[region_id]['strain']
    if strain in strains:
        chrm = regions[region_id]['chromosome']
        start = int(regions[region_id]['start'])
        end = int(regions[region_id]['end'])
        for base in range(start, end + 1):
            bases_by_strains[chrm][base].append(strain)

#for base in sorted(bases_by_strains['I'].keys()):
#    print base, bases_by_strains['I'][base]
        
categories = []
for i in range(1,len(strains) + 1):
    categories += [tuple(sorted(x)) for x in itertools.combinations(strains, i)]

cat_counts = defaultdict(int)
for chrm in bases_by_strains.keys():
    for base in bases_by_strains[chrm]:
        cat_counts[tuple(sorted(bases_by_strains[chrm][base]))] += 1

fn = gp.analysis_out_dir_absolute + tag + '/' + \
     'base_overlap_' + '_'.join(strains) + '.txt'
f = open(fn, 'w')
f.write('group\tcount\n')
for cat in categories:
    f.write(','.join(cat) + '\t' + str(cat_counts[cat]) + '\n')
f.close()

