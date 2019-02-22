import sys
import os
import gzip
import predict
from collections import defaultdict
from summarize_region_quality import *
import gene_predictions
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc/')
import read_fasta
import read_table
import seq_functions

tag = sys.argv[1]

fn = gp.analysis_out_dir_absolute + tag + \
     '/introgressed_blocks_filtered_par_' + tag + '_summary_plus.txt'
regions_filtered, l = read_table.read_table_rows(fn, "\t")

gene_strains = defaultdict(set)
strain_genes = defaultdict(lambda: defaultdict(set))

for chrm in gp.chrms:
    
    fn = gp.analysis_out_dir_absolute + tag + \
         '/genes_for_each_region_chr' + chrm + '_' + \
         tag + '.txt'
    f = open(fn, 'r')
    line = f.readline()
    while line != '':
        line = line.split('\t')
        region_id = line[0]
        if region_id in regions_filtered:
            strain = regions_filtered[region_id]["strain"]
            for gene in line[2::2]:
                gene_strains[gene].add(strain)
                strain_genes[chrm][strain].add(gene)
        line = f.readline()
    f.close()

gene_counts = {}
for gene in gene_strains:
    gene_counts[gene] = len(gene_strains[gene])

f_out = open(gp.analysis_out_dir_absolute + tag + \
             '/genes_for_each_strain_filtered_' + \
             tag + '.txt', 'w')
f_out.write('strain\tchromosome\tnum_genes\n')
for chrm in gp.chrms:
    for strain in strain_genes[chrm]:
        f_out.write(strain + '\t' + chrm + '\t' + \
                    str(len(strain_genes[chrm][strain])) + '\n')
f_out.close()

f_out = open(gp.analysis_out_dir_absolute + tag + \
             '/genes_strain_hist_' + \
             tag + '.txt', 'w')
f_out.write('gene\tnum_strains\n')
for gene in sorted(gene_counts.keys()):
    f_out.write(gene + '\t' + str(gene_counts[gene]) + '\n')
f_out.close()
