## calculate nucleotide diversity for all sites and for all sites
## excluding introgression; also calculate the same but only in coding
## regions

import re
import sys
import os
import copy
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

def try_int(s, default=-1):
    try:
        i = int(s)
        return i
    except:
        return default

def count_diffs(s, t, skip_char = 'N'):
    assert len(s) == len(t)
    num = 0
    den = 0
    for i in range(len(s)):
        if s[i] == skip_char or t[i] == skip_char:
            continue
        if s[i] != t[i]:
            num += 1
        den += 1
    return num, den

## generate a sequence that has the current strain's base for each
## site in the reference sequence, and skip_char for any site where
## the base is a gap/unknown (this is all based on the alignment)
def referize(strain_seq, ref_ind_to_strain_ind, skip_char = 'N'):
    s = [skip_char for r in ref_ind_to_strain_ind]
    for i in range(len(ref_ind_to_strain_ind)):
        si = ref_ind_to_strain_ind[i]
        if si == -1:
            continue
        if strain_seq[si] in 'atgc':
            s[i] = strain_seq[si]
    return s

def mark_excluded(seq, regions, fill='N'):
    seqi = copy.deepcopy(seq)
    for start, end in regions:
        for i in range(start, end+1):
            seqi[i] = fill
    return seqi

def mark_included(seq, regions, fill='N'):
    s = [fill for r in seq]
    for start, end in regions:
        for i in range(start, end+1):
            s[i] = seq[i]
    return s

tag = 'u3_i.001_tv_l1000_f.01'

########
## read in introgressed regions, as well as strains and reference genes
########

## dictionary of introgressed regions keyed by chromosome and then
## strain
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
## read in all strains
strain_dirs = align_helpers.get_strains(gp.non_ref_dirs[gp.master_ref])
num_strains = len(strain_dirs)

## read in genes in reference sequence into dictionary keyed by
## chromosome
ref_genes = {}
for chrm in gp.chrms:
    ref_genes[chrm] = []
    f = open(gp.analysis_out_dir_absolute + gp.master_ref + \
             '_chr' + chrm + '_genes.txt', 'r')
    line = f.readline()
    while line != '':
        line = line[:-1].split('\t')
        ref_genes[chrm].append((int(line[1]), int(line[2])))
        line = f.readline()
    f.close()

########
## calculate nucleotide diversity
########

# all sites
total_frac = 0
total_fracs = dict(zip(gp.chrms, [0 for c in gp.chrms]))

# all sites excluding introgression
total_frac_nonint = 0
total_fracs_nonint = dict(zip(gp.chrms, [0 for c in gp.chrms]))

# all coding sites
total_frac_coding = 0
total_fracs_coding = dict(zip(gp.chrms, [0 for c in gp.chrms]))

# all coding sites excluding introgression
total_frac_coding_nonint = 0
total_fracs_coding_nonint = dict(zip(gp.chrms, [0 for c in gp.chrms]))

# total number of strain pairs
num_comparisons = 0

## loop through all strains
for i in range(num_strains):
    strain_i, d_i = strain_dirs[i]
    strain_i_seqs = {}
    strain_i_seqs_nonint = {}
    strain_i_seqs_coding = {}
    strain_i_seqs_coding_nonint = {}

    ## for each 
    for chrm in gp.chrms:
        ## coordinate conversion between reference and current strain
        coord_fn = gp.analysis_out_dir_absolute + 'coordinates/' + \
                   gp.master_ref + '_to_' + strain_i + \
                   '_chr' + chrm + '.txt.gz'
        f_coord = gzip.open(coord_fn, 'rb')
        ref_ind_to_strain_i_ind = [try_int(line[:-1]) for line in f_coord.readlines()]

        ## current strain fasta file for current chromosome
        strain_fn = d_i + strain_i + '_chr' + chrm + gp.fasta_suffix
        print strain_i, chrm
        
        ## get chromosome sequence for this strain relative to
        ## reference strain (the base for this strain at each site in
        ## the reference, based on original alignment);
        ## gaps/unsequenced sites/etc marked as 'N'
        strain_i_seqs[chrm] = referize(read_fasta.read_fasta(strain_fn)[1][0].lower(),\
                                       ref_ind_to_strain_i_ind)

        ## get version of sequence where everything that doesn't fall
        ## within gene is replaced by 'N'
        strain_i_seqs_coding[chrm] = mark_included(strain_i_seqs[chrm],\
                                                   ref_genes[chrm])

        ## also get version of above sequences where introgressed sites are
        ## replaced by 'N'
        strain_i_seqs_nonint[chrm] = copy.deepcopy(strain_i_seqs[chrm])
        strain_i_seqs_coding_nonint[chrm] = copy.deepcopy(strain_i_seqs_coding[chrm])
        if regions_by_chrm_and_strain[chrm].has_key(strain_i):
            strain_i_seqs_nonint[chrm] = mark_excluded(strain_i_seqs[chrm],\
                                            regions_by_chrm_and_strain[chrm][strain_i])
            strain_i_seqs_coding_nonint[chrm] = \
                mark_excluded(strain_i_seqs_coding[chrm],\
                              regions_by_chrm_and_strain[chrm][strain_i])

    ## loop through all strains to get second strain for current pair
    for j in range(i+1, num_strains):
        strain_j, d_j = strain_dirs[j]

        print strain_i, strain_j
        ## keep track of total number of strain pairs we're looking
        ## at, so we can divide total by that later
        num_comparisons += 1

        num = 0
        den = 0
        num_nonint = 0
        den_nonint = 0
        num_coding = 0
        den_coding = 0
        num_coding_nonint = 0
        den_coding_nonint = 0
        for chrm in gp.chrms:

            ## do the same reading in of sequence for this strain,
            ## relative to reference, and also excluding introgressed
            ## sites
            coord_fn = gp.analysis_out_dir_absolute + 'coordinates/' + \
                       gp.master_ref + '_to_' + strain_j + \
                       '_chr' + chrm + '.txt.gz'
            f_coord = gzip.open(coord_fn, 'rb')
            ref_ind_to_strain_ind = [try_int(line[:-1]) for line in f_coord.readlines()]
            
            strain_fn = d_j + strain_j + '_chr' + chrm + gp.fasta_suffix
            strain_j_seq = referize(read_fasta.read_fasta(strain_fn)[1][0].lower(),\
                                    ref_ind_to_strain_ind)
            strain_j_seq_coding = mark_included(strain_j_seq, ref_genes[chrm])

            strain_j_seq_nonint = copy.deepcopy(strain_j_seq)
            strain_j_seq_coding_nonint = copy.deepcopy(strain_j_seq_coding)
            if regions_by_chrm_and_strain[chrm].has_key(strain_j):
                strain_j_seq_nonint = mark_excluded(strain_j_seq,\
                                        regions_by_chrm_and_strain[chrm][strain_j])
                strain_j_seq_coding_nonint = mark_excluded(strain_j_seq_coding,\
                                        regions_by_chrm_and_strain[chrm][strain_j])

            ## count sites that differ between the two strains
            ## (ignoring any sites where one of the strains has 'N')
            ## and add to appropriate running total

            ## all sites
            num_chrm, den_chrm = count_diffs(strain_i_seqs[chrm], strain_j_seq)
            num += num_chrm
            den += den_chrm
            total_fracs[chrm] += float(num_chrm)/den_chrm

            # nonintrogressed
            num_chrm_nonint, den_chrm_nonint = count_diffs(strain_i_seqs_nonint[chrm],\
                                                           strain_j_seq_nonint)
            num_nonint += num_chrm_nonint
            den_nonint += den_chrm_nonint
            total_fracs_nonint[chrm] += float(num_chrm_nonint)/den_chrm_nonint

            ## all coding sites
            num_chrm_coding, den_chrm_coding = \
                count_diffs(strain_i_seqs_coding[chrm], strain_j_seq_coding)
            num_coding += num_chrm_coding
            den_coding += den_chrm_coding
            total_fracs_coding[chrm] += float(num_chrm_coding)/den_chrm_coding

            # coding, nonintrogressed
            num_chrm_coding_nonint, den_chrm_coding_nonint = \
                count_diffs(strain_i_seqs_coding_nonint[chrm],\
                            strain_j_seq_coding_nonint)
            num_coding_nonint += num_chrm_coding_nonint
            den_coding_nonint += den_chrm_coding_nonint
            total_fracs_coding_nonint[chrm] += \
                float(num_chrm_coding_nonint)/den_chrm_coding_nonint

            print num_comparisons, chrm, \
                total_fracs[chrm], \
                total_fracs_nonint[chrm], \
                1 - total_fracs_nonint[chrm]/total_fracs[chrm], \
                total_fracs_coding[chrm], \
                total_fracs_coding_nonint[chrm], \
                1 - total_fracs_coding_nonint[chrm]/total_fracs_coding[chrm]

        # and keep track across all chromosomes
        total_frac += float(num)/den
        total_frac_nonint += float(num_nonint)/den_nonint
        total_frac_coding += float(num_coding)/den_coding
        total_frac_coding_nonint += float(num_coding_nonint)/den_coding_nonint

        print num_comparisons, total_frac, total_frac_nonint, \
            1 - total_frac_nonint/total_frac, total_frac_coding, \
            total_frac_coding_nonint, 1 - total_frac_coding_nonint/total_frac_coding
        sys.stdout.flush()

# nucleotide diversity is the running total of fractions of sites that
# differ pairwise, divided by the total number of pairs of strains
# we've looked at
nuc_div = total_frac/num_comparisons
nuc_div_nonint = total_frac_nonint/num_comparisons
nuc_div_coding = total_frac_coding/num_comparisons
nuc_div_coding_nonint = total_frac_coding_nonint/num_comparisons

print nuc_div
print nuc_div_nonint
print nuc_div_coding
print nuc_div_coding_nonint

########
## write overall results and results for individual chromosome to file
########

f = open(gp.analysis_out_dir_absolute + tag + '/polymorphism/' + \
         'nucleotide_diversity_c.txt', 'w')
f.write('chromosome\tpi\tpi_nonint\tpi_coding\tpi_coding_nonint\n')
f.write('all\t' + str(nuc_div) + '\t' + str(nuc_div_nonint) + \
        '\t' + str(nuc_div_coding) + '\t' + str(nuc_div_coding_nonint) + '\n')
for chrm in gp.chrms:
    f.write(chrm + '\t' + str(total_fracs[chrm]/num_comparisons) + '\t' + \
            str(total_fracs_nonint[chrm]/num_comparisons) + '\t' + \
            str(total_fracs_coding[chrm]/num_comparisons) + '\t' + \
            str(total_fracs_coding_nonint[chrm]/num_comparisons) + '\n')
f.close()
