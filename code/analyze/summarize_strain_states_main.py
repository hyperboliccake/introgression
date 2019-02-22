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
sys.path.insert(0, '../misc/')
import read_table
import read_fasta

args = predict.process_predict_args(sys.argv[1:])

d = defaultdict(lambda: defaultdict(int))
for species_from in args['known_states'][1:]:

    print species_from

    fn_filtered1i = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
                    'blocks_' + species_from + \
                    '_' + args['tag'] + '_filtered1intermediate.txt'
    fn_filtered2i = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
                    'blocks_' + species_from + \
                    '_' + args['tag'] + '_filtered2intermediate.txt'

    regions1, fields1 = read_table.read_table_rows(fn_filtered1i, '\t')
    regions2, fields2 = read_table.read_table_rows(fn_filtered2i, '\t')

    for region_id in regions1:

        strain = regions1[region_id]['strain']
        length = int(regions1[region_id]['end']) - int(regions1[region_id]['start']) + 1
        d[strain]['num_regions_' + species_from] += 1
        d[strain]['num_regions_total'] += 1
        d[strain]['num_bases_' + species_from] += length
        d[strain]['num_bases_total'] += length
        if regions1[region_id]['reason'] == '':
            d[strain]['num_regions_' + species_from + '_filtered1'] += 1
            d[strain]['num_regions_total_filtered1'] += 1
            d[strain]['num_bases_' + species_from + '_filtered1'] += length
            d[strain]['num_bases_total_filtered1'] += length

            alt_states = regions2[region_id]['alternative_states'].split(',')
            for species_from_alt in alt_states:
                d[strain]['num_regions_' + species_from_alt + \
                          '_filtered2_inclusive'] += 1
                d[strain]['num_bases_' + species_from_alt + \
                          '_filtered2_inclusive'] += length
                if species_from_alt == species_from:
                    d[strain]['num_regions_total_filtered2_inclusive'] += 1
                    d[strain]['num_bases_total_filtered2_inclusive'] += length
                
            if len(alt_states) == 1:
                d[strain]['num_regions_' + species_from + \
                          '_filtered2'] += 1
                d[strain]['num_regions_total_filtered2'] += 1
                d[strain]['num_bases_' + species_from + \
                          '_filtered2'] += length
                d[strain]['num_bases_total_filtered2'] += length


            else:
                d[strain]['num_bases_' + '_or_'.join(sorted(alt_states)) + '_filtered2i'] += length

            d[strain]['num_bases_' + str(len(alt_states)) + '_filtered2i'] += length


strain_info = [line[:-1].split('\t') for line in open('../../100_genomes_info.txt', 'r')]
strain_origins = dict(zip([x[0].lower() for x in strain_info], \
                          [(x[5], x[3], x[4]) for x in strain_info]))
for strain in d.keys():
    d[strain]['population'] = strain_origins[strain][0]
    d[strain]['geographic_origin'] = strain_origins[strain][1]
    d[strain]['environmental_origin'] = strain_origins[strain][2]

fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + 'state_counts_by_strain.txt'
f = open(fn, 'w')
fields = []

fields += ['population', 'geographic_origin', 'environmental_origin']

fields += ['num_regions_' + x for x in args['known_states'][1:]]
fields += ['num_regions_total']
fields += ['num_regions_' + x + '_filtered1' for x in args['known_states'][1:]]
fields += ['num_regions_total_filtered1']
fields += ['num_regions_' + x + '_filtered2' for x in args['known_states'][1:]]
fields += ['num_regions_total_filtered2']
fields += ['num_regions_' + x + '_filtered2_inclusive' for x in args['known_states'][1:]]
fields += ['num_regions_total_filtered2_inclusive']

fields += ['num_bases_' + x for x in args['known_states'][1:]]
fields += ['num_bases_total']
fields += ['num_bases_' + x + '_filtered1' for x in args['known_states'][1:]]
fields += ['num_bases_total_filtered1']
fields += ['num_bases_' + x + '_filtered2' for x in args['known_states'][1:]]
fields += ['num_bases_total_filtered2']
fields += ['num_bases_' + x + '_filtered2_inclusive' for x in args['known_states'][1:]]
fields += ['num_bases_total_filtered2_inclusive']

r = sorted(gp.alignment_ref_order[1:])
for n in range(2, len(r)+1):
    x = itertools.combinations(r, n)
    for combo in x:
        fields += ['num_bases_' + '_or_'.join(combo) + '_filtered2i']
    fields += ['num_bases_' + str(n) + '_filtered2i']

f.write('strain' + '\t' + '\t'.join(fields) + '\n')

for strain in sorted(d.keys()):
    f.write(strain + '\t')
    f.write('\t'.join([str(d[strain][x]) for x in fields]))
    f.write('\n')
f.close()
