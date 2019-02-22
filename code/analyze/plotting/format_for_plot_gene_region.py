# outputs file with all variants in the region
# for each variant site, a symbol for each stain:
# p for matches only par
# c for matches only cer
# n for matches neither
# blank for matches both
# - for gap


import re
import sys
import os
import copy
import gzip
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../align/')
import align_helpers
sys.path.insert(0, '../misc/')
import read_fasta
import read_table

# copy pasta

def try_int(s, default=-1):
    try:
        i = int(s)
        return i
    except:
        return default

def referize(strain_seq, ref_ind_to_strain_ind, skip_char = 'N'):
    s = [skip_char for r in ref_ind_to_strain_ind]
    for i in range(len(ref_ind_to_strain_ind)):
        si = ref_ind_to_strain_ind[i]
        if si == -1:
            continue
        if strain_seq[si] in 'atgc':
            s[i] = strain_seq[si]
    return s

#region_start = 787000
#region_end = 794000
#chrm = 'II'
region_start = 917571 - 100
region_end = 921647 + 100
chrm = 'IV'
region_length = region_end - region_start + 1

##======
# get strains
##======

strain_dirs = align_helpers.get_strains(align_helpers.flatten(gp.non_ref_dirs.values()))
num_strains = len(strain_dirs)


##======
# loop through all strains, getting appropriate sequence
##======

# master reference and other reference seqs
master_ref = gp.alignment_ref_order[0]
master_fn = gp.ref_dir[master_ref] + gp.ref_fn_prefix[master_ref] + '_chr' + \
            chrm + gp.fasta_suffix
master_seq = read_fasta.read_fasta(master_fn)[1][0][region_start:region_end+1].lower()


other_ref = gp.alignment_ref_order[1]
coord_fn = gp.analysis_out_dir_absolute + 'coordinates/' + \
           gp.master_ref + '_to_' + other_ref + \
           '_chr' + chrm + '.txt.gz'
f_coord = gzip.open(coord_fn, 'rb')
ref_ind_to_strain_ind = [try_int(line[:-1]) for line in f_coord.readlines()]
other_ref_fn = gp.ref_dir[other_ref] + gp.ref_fn_prefix[other_ref] + \
               '_chr' + chrm + gp.fasta_suffix
other_ref_seq = referize(read_fasta.read_fasta(other_ref_fn)[1][0].lower(), \
                         ref_ind_to_strain_ind)[region_start:region_end+1]

# other strains
seqs = {}
for i in range(num_strains):
    strain, d = strain_dirs[i]
    print strain
    coord_fn = gp.analysis_out_dir_absolute + 'coordinates/' + \
               gp.master_ref + '_to_' + strain + \
               '_chr' + chrm + '.txt.gz'
    f_coord = gzip.open(coord_fn, 'rb')
    ref_ind_to_strain_ind = [try_int(line[:-1]) for line in f_coord.readlines()]
    strain_fn = d + strain + '_chr' + chrm + gp.fasta_suffix
    seqs[strain] = referize(read_fasta.read_fasta(strain_fn)[1][0].lower(), \
                            ref_ind_to_strain_ind)[region_start:region_end+1]

# write file
fn = 'gene_region_variants.txt'
f = open(fn, 'w')

f.write('ps\t' + '\t'.join([x[0] for x in strain_dirs]) + '\n') 
for i in range(region_length):
    
    f.write(str(region_start + i))
    for strain, d in strain_dirs:
        x = seqs[strain][i]
        f.write('\t')
        if x not in 'atgc' or \
           master_seq[i] not in 'atgc' or \
           other_ref_seq[i] not in 'atgc':
            f.write('-')
        else:
            if x == master_seq[i]:
                if x == other_ref_seq[i]:
                    f.write('b')
                else:
                    f.write('c')
            else:
                if x == other_ref_seq[i]:
                    f.write('p')
                else:
                    f.write('n')
    f.write('\n')
f.close()

