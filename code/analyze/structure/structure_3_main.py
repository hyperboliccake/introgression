## generate three files:

## 1. introgressed regions annotated by which population background(s)
## they overlap

## 2. population backgrounds annotated by how much introgression they
## have from each reference strain (or ambiguous strains)

## 3. counts of bases in for each strain x population background x
## introgresssing reference [or lack of introgression]

import sys
import os
import gzip
import predict
from collections import defaultdict
import gene_predictions
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc/')
import read_fasta
import read_table
import seq_functions

args = predict.process_predict_args(sys.argv[3:])

run_id = sys.argv[1]
num_pops = int(sys.argv[2])
out_dir = gp.analysis_out_dir_absolute + args['tag'] + "/structure/"
out_dir_run = out_dir + run_id + '/'
if run_id == '0':
    out_dir_run = out_dir

# TODO maybe getting strains should be simpler...at least make this
# not copy pasta
strains = [line.split('\t')[0] for line in \
           open(gp.analysis_out_dir_absolute + args['tag'] + \
                '/state_counts_by_strain.txt', 'r').readlines()[1:]]


# returns all the pops that overlap the range [start, end], as well as
# number of bases in each overlap
def find_pops(start, end, pop_ranges):
    pops = []
    bases = []
    in_region = False
    for r in pop_ranges:
        if end >= r[0] and end <= r[1]:
            pops.append(r[2])
            bases.append(end - r[0] + 1)
            break
        if start >= r[0] and start <= r[1]:
            pops.append(r[2])
            bases.append(start - r[0] + 1)
            in_region = True
        elif in_region:
            pops.append(r[2])
            bases.append(r[1] - r[0] + 1)
    return pops, bases

population_int_counts = defaultdict(lambda: defaultdict(int))
strain_population_int_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
population_totals = defaultdict(int)
strain_population_totals = defaultdict(lambda: defaultdict(int))
all_alternative_states = set([])
for ref in gp.alignment_ref_order[1:]:
    regions_fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
                 'blocks_' + ref + \
                 '_' + args['tag'] + '_filtered2intermediate.txt'
    regions, labels = read_table.read_table_rows(regions_fn, '\t')
    regions_strain_chrm = defaultdict(lambda: defaultdict(dict))
    for region_id in regions:
        chrm = regions[region_id]['chromosome']
        strain = regions[region_id]['strain']
        regions_strain_chrm[strain][chrm][region_id] = regions[region_id]
    new_regions_fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
                     'blocks_' + ref + \
                     '_' + args['tag'] + '_populations.txt'
    f = open(new_regions_fn, 'w')
    labels = labels[1:] + ['population']
    f.write('region_id' + '\t' + '\t'.join(labels) + '\n')

    #for chrm in regions_strain_chrm[strain]:
    for strain in strains:
        for chrm in gp.chrms:
            # TODO get rid of run_id in filenames?
            pop_ranges_fn = out_dir_run + 'population_ranges/' + \
                            'population_ranges_' + strain + \
                            '_chr' + chrm + '_run' + run_id + '.txt'
            pop_ranges = [line[:-1].split('\t') for line in \
                          open(pop_ranges_fn, 'r').readlines()]
            pop_ranges = [(int(x[0]), int(x[1]), x[2]) for x in pop_ranges]
            for pr in pop_ranges:
                population_totals[pr[2]] += pr[1] - pr[0] + 1                
                strain_population_totals[strain][pr[2]] += pr[1] - pr[0] + 1

            for region_id in regions_strain_chrm[strain][chrm]:
                r = regions_strain_chrm[strain][chrm][region_id]

                all_alternative_states.add(r['alternative_states'])

                # find the population ranges that the region start and end
                # coordinates fall within
                pops, overlaps = find_pops(int(r['start']), int(r['end']), pop_ranges)
                regions_strain_chrm[strain][chrm][region_id]['population'] = \
                    ','.join(pops)
                f.write(region_id + '\t' + \
                        '\t'.join([str(regions_strain_chrm[strain][chrm][region_id][x])\
                                   for x in labels]) + '\n')

                for i in range(len(pops)):
                    population_int_counts[pops[i]][r['alternative_states']] += \
                        overlaps[i]
                    strain_population_int_counts[strain][pops[i]]\
                        [r['alternative_states']] += overlaps[i]

f = open(out_dir_run + 'population_introgression_counts_run' + run_id + '.txt', 'w')
f.write('population\treference\tnum_bases_introgressed\tfrac_bases_introgressed\n')
for i in population_int_counts.keys():
    for ref in population_int_counts[i].keys():
        f.write(str(i) + '\t' + ref + '\t' + str(population_int_counts[i][ref]) + '\t' + \
                str(float(population_int_counts[i][ref])/population_totals[i]) + '\n')
f.close()


f = open(out_dir_run + 'strain_population_introgression_counts_run' + \
         run_id + '.txt', 'w')
f.write('strain\tpopulation\treference\tnum_bases_introgressed' + \
        '\tfrac_bases_introgressed\n')
for strain in strains:
    for i in strain_population_int_counts[strain].keys():
        for ref in all_alternative_states:
            count = strain_population_int_counts[strain][i][ref]
            total = strain_population_totals[strain][i]
            #frac = 0
            #if total > 0:
            frac = float(count)/total
            f.write(strain + '\t' + str(i) + '\t' + ref + '\t' + 
                    str(count) + '\t' + str(frac) + '\n')
f.close()
