import sys
import os
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../align/')
import align_helpers
sys.path.insert(0, '../misc/')
import read_fasta
import read_table

def pad(s, n=10):
    s = s.strip()
    return s[:n] + (n - len(s)) * ' '


tag = 'u3_i.001_tv_l1000_f.01'
suffix = '_filtered'

gp_dir = '../'
s = align_helpers.get_strains(align_helpers.flatten(gp.non_ref_dirs.values()))

# read in filtered regions
fn_regions = gp.analysis_out_dir_absolute + tag + '/' + \
             'introgressed_blocks' + suffix + '_par_' + tag + '_summary_plus.txt'
regions, l = read_table.read_table_rows(fn_regions, '\t')


site_to_strains_intd = dict(zip(gp.chrms, [{} for c in gp.chrms]))

chrm_lengths = dict(zip(gp.chrms, [0 for c in gp.chrms]))
for region in regions:
    chrm = regions[region]['chromosome']
    region_end = int(regions[region]['end'])
    if region_end >= chrm_lengths[chrm]:
        chrm_lengths[chrm] = region_end + 1

intd = {}
for chrm in gp.chrms:
    intd[chrm] = [set() for site in range(chrm_lengths[chrm])]

# looping through all regions, keep track of which strains are
# introgressed at each site
for region in regions:
    for site in range(int(regions[region]['start']), int(regions[region]['end'])+1):
        chrm = regions[region]['chromosome']
        strain = regions[region]['strain']
        intd[chrm][site].add(strain)

# collapse into consecutive sites with same pattern of introgression
# across the strains
shared_regions = {}
for chrm in gp.chrms:
    shared_regions[chrm] = []
    #if intd[chrm][0] != set():
    site_to_strains_intd[chrm][0] = intd[chrm][0]
    shared_regions[chrm].append([0, 0, intd[chrm][0]])
    for i in range(1, len(intd[chrm])):
        if intd[chrm][i] != intd[chrm][i-1]:
            site_to_strains_intd[chrm][i] = intd[chrm][i]
            shared_regions[chrm].append([i, i, intd[chrm][i]])
        else:
            shared_regions[chrm][-1][1] = i

print shared_regions['I']

strain_matrix = {}
for chrm in gp.chrms:
    for site in site_to_strains_intd[chrm]:
        strains = list(site_to_strains_intd[chrm][site])
        for s1 in range(len(strains)):
            for s2 in range(s1, len(strains)):
                strain1 = strains[s1]
                strain2 = strains[s2]
                if not strain_matrix.has_key(strain1):
                    strain_matrix[strain1] = {}
                if not strain_matrix.has_key(strain2):
                    strain_matrix[strain2] = {}
                if not strain_matrix[strain1].has_key(strain2):
                    strain_matrix[strain1][strain2] = 0
                if not strain_matrix[strain2].has_key(strain1):
                    strain_matrix[strain2][strain1] = 0
                strain_matrix[strain1][strain2] += 1
                if strain1 != strain2:
                    strain_matrix[strain2][strain1] += 1


f = open('strain_shared_introgression_matrix.txt', 'w')
f.write(str(len(strain_matrix.keys())) + '\n')
for strain in sorted(strain_matrix.keys()):
    f.write(pad(strain))
    for strain_other in sorted(strain_matrix.keys()):
        f.write(' ' + str(strain_matrix[strain][strain_other]))
    f.write('\n')
f.close()


f = open('shared_introgression_list.txt', 'w')
f.write('region_number\tchromosome\tstart\tend\tnum_strains\tstrain_list\n')
count = 1
for chrm in gp.chrms:
    for r in shared_regions[chrm]:
        if r[2] != set():
            f.write(str(count) + '\t')
            f.write(chrm + '\t')
            f.write(str(r[0]) + '\t')
            f.write(str(r[1]) + '\t')
            f.write(str(len(r[2])) + '\t')
            f.write(','.join(sorted(list(r[2]))) + '\n')
            count += 1
f.close()

f = open('shared_introgression_nonsingleton_list.txt', 'w')
f.write('region_number\tchromosome\tstart\tend\tnum_strains\tstrain_list\n')
count_ns = 1
for chrm in gp.chrms:
    for r in shared_regions[chrm]:
        if len(r[2]) > 1:
            f.write(str(count_ns) + '\t')
            f.write(chrm + '\t')
            f.write(str(r[0]) + '\t')
            f.write(str(r[1]) + '\t')
            f.write(str(len(r[2])) + '\t')
            f.write(','.join(sorted(list(r[2]))) + '\n')
            count_ns += 1

f = open('strain_shared_introgression.txt', 'w')
f.write(str(len(strain_matrix.keys())) + ' ' + str(count - 1) + '\n')
for strain in sorted(strain_matrix.keys()):
    f.write(pad(strain))
    for chrm in gp.chrms:
        for i in sorted(site_to_strains_intd[chrm].keys()):
            if site_to_strains_intd[chrm][i] != set():
                if strain in site_to_strains_intd[chrm][i]:
                    f.write('1')
                else:
                    f.write('0')
    f.write('\n')
f.write(pad('S288c') + '0' * (count - 1) + '\n')
f.write(pad('CBS432') + '0' * (count - 1) + '\n')
f.close()

f = open('strain_shared_introgression_nonsingleton.txt', 'w')
f.write(str(len(strain_matrix.keys())) + ' ' + str(count_ns - 1) + '\n')
for strain in sorted(strain_matrix.keys()):
    f.write(pad(strain))
    for chrm in gp.chrms:
        for i in sorted(site_to_strains_intd[chrm].keys()):
            if len(site_to_strains_intd[chrm][i]) > 1:
                if strain in site_to_strains_intd[chrm][i]:
                    f.write('1')
                else:
                    f.write('0')
    f.write('\n')
f.write(pad('S288c') + '0' * (count_ns - 1) + '\n')
f.write(pad('CBS432') + '0' * (count_ns - 1) + '\n')
f.close()


# for chromosomes separately
for chrm in gp.chrms:

    f = open('shared_introgression_chr' + chrm + '_nonsingleton_list.txt', 'w')
    f.write('region_number\tchromosome\tstart\tend\tnum_strains\tstrain_list\n')
    count_ns = 1

    for r in shared_regions[chrm]:
        if len(r[2]) > 1:
            f.write(str(count_ns) + '\t')
            f.write(chrm + '\t')
            f.write(str(r[0]) + '\t')
            f.write(str(r[1]) + '\t')
            f.write(str(len(r[2])) + '\t')
            f.write(','.join(sorted(list(r[2]))) + '\n')
            count_ns += 1

    f.close()

    f = open('strain_shared_introgression_chr' + chrm + '_nonsingleton.txt', 'w')
    f.write(str(len(strain_matrix.keys())) + ' ' + str(count_ns - 1) + '\n')
    for strain in sorted(strain_matrix.keys()):
        f.write(pad(strain))
        for i in sorted(site_to_strains_intd[chrm].keys()):
            if len(site_to_strains_intd[chrm][i]) > 1:
                if strain in site_to_strains_intd[chrm][i]:
                    f.write('1')
                else:
                    f.write('0')
        f.write('\n')
    f.write(pad('S288c') + '0' * (count_ns - 1) + '\n')
    f.write(pad('CBS432') + '0' * (count_ns - 1) + '\n')
    f.close()

