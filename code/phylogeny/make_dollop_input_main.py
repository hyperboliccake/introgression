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

# read in genes for each region
fn_genes_for_each_region = gp.analysis_out_dir_absolute + tag + '/' + \
                           'genes_for_each_region_' + tag + '.txt'
introgressed_genes_strains = {}
for line in open(fn_genes_for_each_region, 'r').readlines():
    line = line.split('\t')
    region = line[0]
    if region in regions:
        strain = regions[region]['strain']
        g = line[2::2]
        for gene in g:
            if not introgressed_genes_strains.has_key(gene):
                introgressed_genes_strains[gene] = set([])
            introgressed_genes_strains[gene].add(strain)

f_gene_list = open('dollop_gene_list.txt', 'w')
genes = sorted(introgressed_genes_strains.keys())
for gene in genes:
    f_gene_list.write(gene + '\n')
f_gene_list.close()

f_dollop = open('infile_dollop', 'w')

f_dollop.write(str(len(s)+2) + ' ' + str(len(introgressed_genes_strains)) + '\n')
f_dollop.write(pad(gp.master_ref) + '0' * len(genes) + '\n')
f_dollop.write(pad(gp.alignment_ref_order[1]) + '0' * len(genes) + '\n')

for strain, d in s:
    f_dollop.write(pad(strain))
    for gene in genes:
        if strain in introgressed_genes_strains[gene]:
            f_dollop.write('1')
        else:
            f_dollop.write('0')
    f_dollop.write('\n')

f_dollop.close()

f_matrix = open('strain_gene_introgressed_matrix.tsv', 'w')
f_matrix.write('strain' + '\t' + '\t'.join(genes) +'\n')
f_matrix.write(gp.master_ref + '\t' + '\t'.join(['0' for g in genes]) + '\n')
f_matrix.write(gp.alignment_ref_order[1] + '\t' + '\t'.join(['0' for g in genes]) + '\n')
for strain, d in s:
    f_matrix.write( strain)
    for gene in genes:
        if strain in introgressed_genes_strains[gene]:
            f_matrix.write('\t1')
        else:
            f_matrix.write('\t0')
    f_matrix.write('\n')
f_matrix.close()
