# compare set of genes I've called to set called in Strope et al (100
# genomes paper)

import re
import sys
import os
import math
import Bio.SeqIO
import copy
import gene_predictions
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../align/')
import align_helpers
sys.path.insert(0, '../misc/')
import read_table
import read_fasta
import write_fasta
import mystats

tag = sys.argv[1]

fn_genes_strope = '../../data/introgressed_genes_strope.tsv'
f_strope = open(fn_genes_strope, 'r')
lines = [line.strip().split('\t') for line in f_strope.readlines()]
f_strope.close()
genes_strope = {}
sys_standard_strope = {}
strains = lines[0][7:-2]
for line in lines[1:]:
    n_int_other = lines[7:-1].count('s')
    n_del = lines[7:-1].count('0')
    strains_int_par = []
    for i in range(len(strains)):
        if line[7+i] == 'P':
            strains_int_par.append(strains[i])
    n_int_par = len(strains_int_par)
    genes_strope[line[2]] = (n_int_par, n_int_other, n_del, strains_int_par, \
                             line[1], line[4])
    sys_standard_strope[line[1]] = line[2]
    
fn_regions = gp.analysis_out_dir_absolute + tag + '/' + \
             'introgressed_blocks_filtered_par_' + tag + '_summary_plus.txt'
# dict keyed by region: {strain:, start:, end:, etc}
regions, l = read_table.read_table_rows(fn_regions, '\t')
region_to_genes = {}
for chrm in gp.chrms:
    fn_genes_regions = gp.analysis_out_dir_absolute + tag + '/' + \
                       'genes_for_each_region_chr' + chrm + '_' + tag + '.txt'
    # keyed by region: {'num_genes':,'gene_list'}
    region_to_genes_current = \
        gene_predictions.read_genes_for_each_region_summary(fn_genes_regions)
    region_to_genes.update(region_to_genes_current)
genes_by_strain = {}
for region in regions:
    if not genes_by_strain.has_key(regions[region]['strain']):
        genes_by_strain[regions[region]['strain']] = set([])
    [genes_by_strain[regions[region]['strain']].add(gene) \
     for gene in [x[0] for x in region_to_genes[region]['gene_list']]]

genes = {}
for chrm in gp.chrms:
    fn_genes = gp.analysis_out_dir_absolute + '/' + tag + '/' + \
               'strains_for_each_gene_chr' + chrm + '_' + tag + '.txt'
    f_genes = open(fn_genes, 'r')
    lines = [line.strip().split('\t') for line in f_genes.readlines()]
    f_genes.close()
    for line in lines:
        gene = line[0]
        current_strains = []
        current_fracs = []
        strains = line[2::2]
        fracs = line[3::2]
        for i in range(len(strains)):
            if gene in genes_by_strain[strains[i]]:
                current_strains.append(strains[i])
                current_fracs.append(float(fracs[i]))
        num_strains = len(current_strains)
        if num_strains > 0:
            genes[line[0]] = (num_strains, strains, fracs)


fn_all_genes = '../../data/S288c_verified_orfs.tsv'
f_all_genes = open(fn_all_genes, 'r')
lines = [line.strip().split('\t') for line in f_all_genes.readlines()]
f_all_genes.close()
all_genes = {}
sys_standard = {}
for line in lines:
    chrm = line[6][3:]
    start = int(line[4])
    end = int(line[5])
    strand = line[3]
    all_genes[line[1]] = (line[0], chrm, start, end, strand)
    sys_standard[line[0]] = line[1]

# TODO fix my gene list then get rid of this
all_genes = {}
for chrm in gp.chrms:
    fn_all_genes = gp.analysis_out_dir_absolute + 'S288c_chr' + chrm + '_genes.txt'
    f_all_genes = open(fn_all_genes, 'r')
    lines = [line.strip().split('\t') for line in f_all_genes.readlines()]
    f_all_genes.close()
    for line in lines:
        start = int(line[1])
        end = int(line[2])
        strand = 'NA'
        all_genes[line[0]] = ('NA', chrm, start, end, strand)

    
fn_paralogs = '../../data/S288c_paralogs.tsv'
f_paralogs = open(fn_paralogs, 'r')
lines = [line.strip().split('\t') for line in f_paralogs.readlines()]
f_paralogs.close()
paralogs = {}
for line in lines:
    if line[0] != "":
        paralogs[line[0]] = line[3]


f_s = open('compare_to_strope/genes_strope_only.txt', 'w')
f_m = open('compare_to_strope/genes_me_only.txt', 'w')
f_sm = open('compare_to_strope/genes_both.txt', 'w')
f_mp = open('compare_to_strope/genes_me_paralogs.txt', 'w')
f_sp = open('compare_to_strope/genes_strope_paralogs.txt', 'w')
c_s = 0
c_m = 0
c_sm = 0
c_s_p = 0
c_m_p = 0
c_sm_p = 0
for gene in genes:
    if gene in paralogs:
        f_mp.write(gene + '\n')
    if gene in genes_strope or all_genes[gene][0] in genes_strope:
        f_sm.write(gene + '\n')
        c_sm += 1
        if gene in paralogs:
            c_sm_p += 1
    else:
        f_m.write(gene + '\n')
        c_m += 1
        if gene in paralogs:
            c_m_p += 1
for gene in genes_strope:
    if gene in paralogs:
        f_sp.write(gene + '\n')
    if gene in genes or (gene in sys_standard and sys_standard[gene] in genes):
        continue
    elif not (gene in all_genes or (gene in sys_standard and sys_standard[gene] in all_genes)):
        continue
    elif genes_strope[gene][0] == 0:
        continue
    elif genes_strope[gene][-1] != 'Verified':
        continue
    else:
        f_s.write(gene + '\n')
        c_s += 1
        if gene in paralogs:
            c_s_p +=1
f_s.close()
f_m.close()
f_sm.close()
f_mp.close()
f_sp.close()

print 'number strope only:', c_s
print 'number me only:', c_m
print 'number strope and me:', c_sm
print 'number strope only paralogs', c_s_p
print 'number me only paralogs', c_m_p
print 'number strope and me paralogs', c_sm_p
print 'number paralogs', len(paralogs)

print paralogs.keys()
