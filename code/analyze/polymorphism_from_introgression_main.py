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

def convert_and_get(i1, l, s, default_value='n'):
    try:
        i2 = int(l[i1])
        return s[i2]
    except:
        return default_value

# loop through each site in cerevisiae reference
# count sites at which there is a variant in at least one other strain
# count sites where the only variants present are in introgressed regions
# count sites where the only variants present are in introgressed regions and match paradoxus

# and for comparison
# count variable sites where not all variants are in introgressed regions
# count variable sites where there's no introgression
# count sites where only variants are in introgressed regions but don't match paradoxus

tag = 'u3_i.001_tv_l1000_f.01'

regions_by_chrm_and_strain = dict(zip(gp.chrms, [{} for i in range(len(gp.chrms))]))
fn_regions = gp.analysis_out_dir_absolute + tag + '/' + \
             'introgressed_blocks_filtered_par_' + tag + '_summary_plus.txt'
d, labels = read_table.read_table_rows(fn_regions, '\t')

for region in d:
    chrm = d[region]['chromosome']
    strain = d[region]['strain']
    if not regions_by_chrm_and_strain[chrm].has_key(strain):
        regions_by_chrm_and_strain[chrm][strain] = []
    regions_by_chrm_and_strain[chrm][strain].append((int(d[region]['start']), \
                                                     int(d[region]['end'])))

# storing a venn diagram as a dictionary, with the key being
# a tuple with a value for each of these categories:
# 1. number of alleles in cerevisiae (1, 2, 3+)
# 2. cer ref == par ref (true, false)
# 3. >= 1 strain with paradoxus variant introgressed (true, false)
# 4. all strains with paradoxus variant introgressed (true, false)
x = list(itertools.product([True, False], repeat=3))
keys = []
for num_alleles in ['1', '2', '3+']:
    for xi in x:
        k = (num_alleles,) + xi
        keys.append(k)

v_by_chrm = {}

f = open(gp.analysis_out_dir_absolute + tag + '/polymorphism/polymorphism.txt', 'w')
f.write('chromosome\tnumber_cer_alleles\tcer_ref_match_par_ref\tone_match_intd\tall_match_intd\tcount\n')

f_sites = open(gp.analysis_out_dir_absolute + tag + '/polymorphism/polymorphism_with_sites.txt', 'w')
f_sites.write('chromosome\tnumber_cer_alleles\tcer_ref_match_par_ref\tone_match_intd\tall_match_intd\tcount\tsites\n')

total = 0
total_r2_n = 0

for chrm in gp.chrms:
    
    fn = gp.ref_dir[gp.master_ref] + gp.ref_fn_prefix[gp.master_ref] + \
         '_chr' + chrm + gp.fasta_suffix
    ref_seq = read_fasta.read_fasta(fn)[1][0].lower()

    ref2 = filter(lambda x: x != gp.master_ref, gp.alignment_ref_order)[0]
    fn = gp.ref_dir[ref2] + gp.ref_fn_prefix[ref2] + \
         '_chr' + chrm + gp.fasta_suffix
    ref2_seq = read_fasta.read_fasta(fn)[1][0].lower()
    coord_fn = gp.analysis_out_dir_absolute + 'coordinates/' + \
               gp.master_ref + '_to_' + ref2 + \
               '_chr' + chrm + '.txt.gz'
    f_coord = gzip.open(coord_fn, 'rb')
    ref_ind_to_ref2_ind = [line[:-1] for line in f_coord.readlines()]
    # TODO
    #m = max(ref_ind_to_ref2_ind)
    #ref_ind_to_ref2_ind = [x if x < m else x - 1 for x in ref_ind_to_ref2_ind]

    alleles = [[] for base in ref_seq]
    strain_dirs = align_helpers.get_strains(gp.non_ref_dirs[gp.master_ref])
    num_strains = len(strain_dirs)

    for strain, d in strain_dirs:
        print chrm, strain
        coord_fn = gp.analysis_out_dir_absolute + 'coordinates/' + \
                   gp.master_ref + '_to_' + strain + \
                   '_chr' + chrm + '.txt.gz'
        f_coord = gzip.open(coord_fn, 'rb')
        ref_ind_to_strain_ind = [line[:-1] for line in f_coord.readlines()]
        # TODO
        #m = max(ref_ind_to_strain_ind)
        #ref_ind_to_strain_ind = [x if x < m else x - 1 for x in ref_ind_to_strain_ind]

        strain_fn = d + strain + '_chr' + chrm + gp.fasta_suffix
        strain_seq = read_fasta.read_fasta(strain_fn)[1][0].lower()

        for i in range(len(ref_seq)):
            alleles[i].append(convert_and_get(i, ref_ind_to_strain_ind, strain_seq))
            
    v = {}
    sites = {}
    for k in keys:
        v[k] = 0
        sites[k] = []

    # throw out sites where cer ref unsequenced (but leave sites where just par is)
    introgressed_regions = regions_by_chrm_and_strain[chrm]
    for i in range(len(ref_seq)):
        r = ref_seq[i]
        if r == 'n':
            continue
        r2 = convert_and_get(i, ref_ind_to_ref2_ind, ref2_seq)
        if r2 == 'n':
            total_r2_n += 1
        total += 1

        # 1
        a = set([r] + alleles[i]).difference(set(['n']))
        num_alleles = '2'
        if len(a) == 1:
            num_alleles = '1'
        if len(a) > 2:
            num_alleles = '3+'
           
        # 2
        refs_match = False
        if r == r2:
            refs_match = True

        # 3 & 4
        one_match = False
        one_match_intd = False
        all_match_intd = True
        for s in range(num_strains):
            if alleles[i][s] == r2:
                one_match = True
                if introgressed_regions.has_key(strain_dirs[s][0]) and \
                   overlap.contained_any(i, introgressed_regions[strain_dirs[s][0]]):
                    one_match_intd = True
                else:
                    all_match_intd = False
        if not one_match:
            all_match_intd = False

        k = (num_alleles, refs_match, one_match_intd, all_match_intd)
        v[k] += 1
        sites[k].append(i)

    for k in keys:
        f.write(chrm + '\t')
        f.write('\t'.join([str(x) for x in k]))
        f.write('\t' + str(v[k]) + '\n')
        
        f_sites.write(chrm + '\t')
        f_sites.write('\t'.join([str(x) for x in k]))
        f_sites.write('\t' + str(v[k]) + '\t')
        f_sites.write(','.join([str(x) for x in sites[k]]) + '\n')

    v_by_chrm[chrm] = v

print 'total sites considered:', total
print 'of sites considered, number at which paradoxus not sequenced/aligned:', total_r2_n

for k in keys:
    v_all = 0
    for chrm in gp.chrms:
        v_all += v_by_chrm[chrm][k]
    f.write('all' + '\t')
    f.write('\t'.join([str(x) for x in k]))
    f.write('\t' + str(v_all) + '\n')

f.close()
f_sites.close()
