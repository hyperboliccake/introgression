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

args = predict.process_predict_args(sys.argv[2:])

run_id = sys.argv[1]
out_dir = gp.analysis_out_dir_absolute + args['tag'] + "/structure/"
out_dir_run = out_dir + run_id + '/'
if run_id == '0':
    out_dir_run = out_dir
if not os.path.isdir(out_dir_run):
    os.makedirs(out_dir_run)
    os.makedirs(out_dir_run + '/population_ranges')

# maybe getting strains should be simpler
strains = [line.split('\t')[0] for line in \
           open(gp.analysis_out_dir_absolute + args['tag'] + \
                '/state_counts_by_strain.txt', 'r').readlines()[1:]]

gp_dir = '../'

nuc_to_int = {'a':1, 't':2, 'g':3, 'c':4}

##======
# use program structure to find population proportion using either
# unlinked tagsnps from ldselect, or just all snps
##======

use_all_snps = True

all_snps = {}
strains = set([])
# map snp indices from program to chromosome and position
ind_to_chrm_ps = []
for chrm in gp.chrms:
    f = open(out_dir + 'ldselect_input_chr' + chrm + '.tsv', 'r')
    # keep list of all strain alleles for each snp (these should be in
    # order in the file we're reading)
    snp_strains = defaultdict(dict)
    line = f.readline()
    prev_snp_id = -1
    while line != '':
        line = line[:-1].split('\t')
        strain = line[1]
        strains.add(strain)
        snp_id = line[0]
        # keep only the snps that we've identified as the tag snps we're
        # going to use
        # if snp_id in all_tag_snps:
        snp_strains[int(snp_id)][strain] = nuc_to_int[line[-1]]
        if snp_id != prev_snp_id:
            ind_to_chrm_ps.append((chrm, int(snp_id)))
            prev_snp_id = snp_id
        line = f.readline()
    f.close()
    all_snps[chrm] = snp_strains

map_distances = {}
for chrm in gp.chrms:
    map_distances[chrm] = {}
    snps = sorted(all_snps[chrm].keys())
    map_distances[chrm][snps[0]] = -1
    for i in range(1, len(snps)):
        map_distances[chrm][snps[i]] = snps[i] - snps[i-1]

strains = sorted(list(strains))

f = open('../../100_genomes_info.txt', 'r')
lines = [line[:-1].split('\t') for line in f.readlines()]
f.close()
strain_pops = defaultdict(int)
pops = sorted(list(set([line[-1] for line in lines])))
pops.remove('mosaic')
pop_inds = dict(zip(pops, [str(x) for x in range(1, len(pops) + 1)]))
pop_inds['mosaic'] = 0

for line in lines:
    strain_pops[line[0].lower()] = pop_inds[line[-1]]

# write input file for structure, with snps as columns and strains as
# rows (snp name format is chr_pos or TODO whatever we've made it
# above); also add population based on 100 genomes paper

f = open(out_dir_run + 'structure_input_run' + run_id + '.txt', 'w')
for chrm in gp.chrms:
    f.write('\t\t\t' + '\t'.join([chrm + '_' + str(x) \
                              for x in sorted(all_snps[chrm].keys())]))
f.write('\n')
for chrm in gp.chrms:
    f.write('\t\t\t' + '\t'.join([str(map_distances[chrm][x]) \
                              for x in sorted(map_distances[chrm].keys())]))
f.write('\n')

for strain in strains:
    flag = '1'
    if strain_pops[strain] == 0:
        flag = '0'
    f.write(strain + '\t' + str(strain_pops[strain]) + '\t' + flag)
    for chrm in gp.chrms:
        for snp in sorted(all_snps[chrm].keys()):
            f.write('\t' + str(str(all_snps[chrm][snp][strain])))
    f.write('\n')
f.close()

num_snps = sum([len(all_snps[chrm]) for chrm in gp.chrms])

# run structure (with varying k values? or just k=6?)
"""
os.system(gp.structure_install_path + 'structure -L ' + str(num_snps) + \
          ' -K 6 -i ' + out_dir_run + 'structure_input_run' + run_id + \
          '.txt -o ' + out_dir_run + 'structure_output_k6_run' + run_id + '.txt')

os.system('mv ' + out_dir_run + 'structure_output_k6_run' + \
          run_id + '.txt_ss ' + out_dir_run + \
          'structure_output_ss_k6_run' + run_id + '.txt')

os.system('mv ' + out_dir_run + 'structure_output_k6_run' + \
          run_id + '.txt_f ' + out_dir_run + \
          'structure_output_k6_run' + run_id + '.txt')
"""
# generate output file with population fractions for each strain (for
# making typical structure barplot)

f_out = open(out_dir_run + 'strain_pop_fracs_k6_run' + run_id + '.txt', 'w')
f = open(out_dir_run + 'structure_output_k6_run' + run_id + '.txt', 'r')
line = f.readline()
while line != "Inferred ancestry of individuals:\n":
    line = f.readline()
f.readline() # column headings
line = f.readline()
f_out.write('strain\tpopulation\tfraction\tindex\n')
while line != "\n":
    line = line.split()
    strain = line[1]
    fracs = [float(x) for x in line[line.index(':')+1:]]
    ind = len(fracs)
    for i in range(len(fracs)):
        if fracs[i] > .6:
            ind = i
            break
    for i in range(len(fracs)):
        f_out.write(strain + '\t' + str(i + 1) + '\t' + \
                    str(fracs[i]) + '\t' + str(ind + 1) + '\n')
    line = f.readline()
f.close()
f_out.close()

# generate file with structure population coordinates for each strain
# x chrm

f = open(out_dir_run + 'structure_output_ss_k6_run' + run_id + '.txt', 'r')
k = 6

# read in posterior probabilities for each strain locus being in each population
line = f.readline()
while line.strip() == '\n':
    line = f.readline()
strain_snp_pop = defaultdict(lambda: defaultdict(dict))
while line != '':
    if line == '\n':
        line = f.readline()
        continue

    line = line.split()
    strain = strains[int(line[0]) - 1]
    chrm, ps = ind_to_chrm_ps[int(line[1]) - 1]

    pop_fracs = [float(x) for x in line[2:]]
    max_frac = 0
    max_ind = 0
    for fi in range(k):
        if pop_fracs[fi] > max_frac:
            max_frac = pop_fracs[fi]
            max_ind = fi + 1

    strain_snp_pop[strain][chrm][ps] = str(max_ind)
    line = f.readline()
f.close()



# TODO at some point associate numbered populations with logical names
# (i.e. ones from strope et al)

# generate file that groups consecutive snps with same population into
# ranges

# population_ranges_strain_chrX.txt
# start end popx
# start end popx/popy
# start end 

chrm_lengths = [line[:-1].split('\t') for line in \
                     open(out_dir + 'chromosome_lengths.txt', 'r').readlines()]
chrm_lengths = dict(zip([x[0] for x in chrm_lengths], \
                        [int(x[1]) for x in chrm_lengths]))

for strain in strains:
    for chrm in gp.chrms:
        ranges = []
        snps =  sorted(strain_snp_pop[strain][chrm].keys())
        start = snps[0]
        end = start
        previous_pop = strain_snp_pop[strain][chrm][start]
        ranges.append((0, start - 1, 'start'))
        for snp in snps:
            current_pop = strain_snp_pop[strain][chrm][snp]

            if current_pop == previous_pop:
                end = snp

            else:
                ranges.append((start, end, previous_pop))
                ranges.append((end + 1, snp - 1, previous_pop + '/' + current_pop))
                start = snp
                end = snp
                previous_pop = current_pop

        ranges.append((start, end, previous_pop))
        # TODO get chromosome lengths
        ranges.append((end + 1, chrm_lengths[chrm], 'end'))

        # TODO file location
        f = open(out_dir_run + 'population_ranges/population_ranges_' + strain + '_chr' + chrm + '_run' + run_id + '.txt', 'w')
        for r in ranges:
            f.write('\t'.join([str(x) for x in r]) + '\n')
        f.close()
