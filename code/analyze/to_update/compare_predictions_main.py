import re
import sys
import os
import copy
import itertools
import gene_predictions
import predict
from collections import defaultdict
from filter_helpers import *
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../sim/')
import sim_analyze_hmm_bw as sim
sys.path.insert(0, '../misc/')
import mystats
import read_table
import read_fasta

# similar to find_pops function in structure_3_main.py
def overlap_with_any(start, end, blocks):
    count = 0
    for block in blocks:
        if block[0] <= start and start <= block[1]:
            count += min(end, block[1]) - start + 1
            if end <= block[1]:
                break
        if start < block[0] and end > block[1]:
            count += block[1] - block[0] + 1
            continue
        if block[0] <= end and end <= block[1]:
            count += end - block[0] + 1
            break
    return count

args = predict.process_predict_args(sys.argv[1:])

## comparing to other prediction run; e.g. comparing using just one
## introgressed reference state to using multiple; this is a little
## janky because some of the file names and formatting have changed
other_region_fn = gp.analysis_out_dir_absolute + 'u3_i.001_tv_l1000_f.01/' + \
                  'introgressed_blocks_filtered_par_u3_i.001_tv_l1000_f.01_summary_plus.txt'
rt_other, fields_other = read_table.read_table_rows(other_region_fn, '\t')
regions_other = defaultdict(lambda: defaultdict(list))
for region_id in rt_other:
    chrm = rt_other[region_id]['chromosome']
    strain = rt_other[region_id]['strain']
    regions_other[chrm][strain].append((int(rt_other[region_id]['start']), \
                                        int(rt_other[region_id]['end'])))
for chrm in gp.chrms:
    for strain in regions_other[chrm].keys():
        regions_other[chrm][strain].sort(key = lambda x: x[0])

regions = defaultdict(lambda: defaultdict(list))

for species_from in args['known_states'][1:]:

    fn_filtered2i = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
                    'blocks_' + species_from + \
                    '_' + args['tag'] + '_filtered2intermediate.txt'

    rt, fields = read_table.read_table_rows(fn_filtered2i, '\t')

    for region_id in rt:
        chrm = rt[region_id]['chromosome']
        strain = rt[region_id]['strain']
        regions[chrm][strain].append((int(rt[region_id]['start']), \
                                      int(rt[region_id]['end']),
                                      rt[region_id]['alternative_states']))
for chrm in gp.chrms:
    for strain in regions[chrm].keys():
        regions[chrm][strain].sort(key = lambda x: x[0])

# count bases found in every possible combination of species_from +
# presence/absence in regions_other
# e.g. d[(other, ref2, ref3)] = x or d[(ref1,)] = y
d = defaultdict(lambda: defaultdict(int))

for chrm in gp.chrms:
    # current predictions
    for strain in regions[chrm].keys():
        for region in regions[chrm][strain]:
            x = overlap_with_any(region[0], region[1], regions_other[chrm][strain])
            length = region[1] - region[0] + 1
            alt_states = region[2].split(',')
            d[strain][tuple(['other'] + alt_states)] += x
            d[strain][tuple(alt_states)] += length - x
            assert x <= length


    # other predictions
    for strain in regions_other[chrm].keys():
        for region in regions_other[chrm][strain]:
            x = overlap_with_any(region[0], region[1], regions[chrm][strain])
            length = region[1] - region[0] + 1
            d[strain][('other', 'any')] += x
            d[strain][('other',)] += length - x
            assert x <= length
            

fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + 'state_counts_comparison.txt'
f = open(fn, 'w')

f.write('strain\tlabel\tcount\n')
for strain in d.keys():
    for label in d[strain].keys():
        f.write(strain + '\t' + ','.join(label) + '\t' + str(d[strain][label]) + '\n')
f.close()
