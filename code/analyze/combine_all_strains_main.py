# input a gene or start/end coordinates
# output a multiple alignment file
# - for gene, relies on annotations/orfs
# - for coordinates, relies on alignments

import re
import sys
import os
import math
import Bio.SeqIO
import copy
from combine_all_strains import *
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
gene = sys.argv[2]
chrm = sys.argv[3]

print 'getting gene sequence from reference strain'
ref = 'S288c'
ref_gene_coords_fn = gp.analysis_out_dir_absolute + ref + '_chr' + chrm + '_genes.txt'
ref_seq_fn = gp.ref_dir[ref] + gp.ref_fn_prefix[ref] + '_chr' + chrm + gp.fasta_suffix
ref_gene_seq = get_ref_gene_seq(gene, ref_gene_coords_fn, ref_seq_fn)
query_fn = gene + '.txt'
f = open(query_fn, 'w')
f.write(ref_gene_seq + '\n')
f.close()

print 'getting gene sequences from all strains'
s = align_helpers.get_strains(align_helpers.flatten(gp.non_ref_dirs.values()))
s.append((gp.ref_fn_prefix[gp.alignment_ref_order[1]], gp.ref_dir[gp.alignment_ref_order[1]]))
strain_gene_seqs = get_gene_seqs(query_fn, s, chrm)

print 'writing all gene sequences to file'
keys = sorted(strain_gene_seqs.keys())
headers = [key + ' ' + strain_gene_seqs[key][0] for key in keys]
seqs = [strain_gene_seqs[key][1] for key in keys]
strains = [ref] + keys
headers = [ref + ' ' + gene] + headers
seqs = [ref_gene_seq] + seqs
gene_seqs_fn = gp.analysis_out_dir_absolute + tag + '/genes/' + gene + '/' + \
               gene + gp.fasta_suffix
if not os.path.isdir(gp.analysis_out_dir_absolute + tag + '/genes/' + gene):
    os.makedirs(gp.analysis_out_dir_absolute + tag + '/genes/' + gene)
write_fasta.write_fasta(headers, seqs, gene_seqs_fn)

print 'finding (filtered) regions that overlap gene'
# read in filtered regions
fn_regions = gp.analysis_out_dir_absolute + tag + '/' + \
             'introgressed_blocks_filtered_par_' + tag + '_summary_plus.txt'
regions, l = read_table.read_table_rows(fn_regions, '\t')

# figure out which strains are introgressed/which regions overlap gene
print 'finding introgressed regions overlapping gene'
fn_genes_regions = gp.analysis_out_dir_absolute + tag + '/' + \
                   'genes_for_each_region_chr' + chrm + '_' + tag + '.txt'
region_to_genes = gene_predictions.read_genes_for_each_region_summary(fn_genes_regions)
#strains = [x[0] for x in s]
regions_overlapping = {}
# TODO does this actually ensure that regions are sorted appropriately
# in fasta headers below?
region_keys_ordered = sorted(regions.keys(), key=lambda x: int(x[1:]))
for region in region_keys_ordered:
    if regions[region]['chromosome'] == chrm and \
       gene in [x[0] for x in region_to_genes[region]['gene_list']]:
        strain = regions[region]['strain']
        if not regions_overlapping.has_key(strain):
            regions_overlapping[strain] = []
        regions_overlapping[strain].append(region)

print 'writing all gene sequences to file, with introgressed bases capitalized'
for i in range(len(seqs)):
    strain = strains[i]
    seq = seqs[i].lower()
    header = headers[i]
    seqs[i] = seq
    headers[i] = header
    if strain not in regions_overlapping:
        continue
    g = strain_gene_seqs[strain]
    t, labels = read_table.read_table_columns(gp.analysis_out_dir_absolute + \
                                              tag + '/' + \
                                              'site_summaries/predictions_' + strain + \
                                              '_chr' + chrm + \
                                              '_site_summary.txt.gz', '\t')
    ref_ind_to_strain_ind = dict(zip(t['ps_ref'], t['ps_strain']))
    for region in regions_overlapping[strain]:
        header += ' ' + region
        start_strain = math.ceil(float(ref_ind_to_strain_ind[regions[region]['start']]))
        end_strain = math.floor(float(ref_ind_to_strain_ind[regions[region]['end']]))
        start_relative = int(max(start_strain - int(g[2]), 0))
        end_relative = int(end_strain - int(g[2]))
        seq = seq[:start_relative] + \
              seq[start_relative:end_relative+1].upper() + \
              seq[end_relative+1:]
    seqs[i] = seq
    headers[i] = header

gene_seqs_fn = gp.analysis_out_dir_absolute + tag + '/genes/' + gene + '/' + gene + \
               '_introgressed' + gp.fasta_suffix
write_fasta.write_fasta(headers, seqs, gene_seqs_fn)


print 'aligning gene sequences'
gene_seqs_aligned_fn = gene_seqs_fn.replace(gp.fasta_suffix, gp.alignment_suffix)
cmd_string = gp.mafft_install_path + '/mafft ' + \
             ' --reorder --preservecase ' + \
             gene_seqs_fn + ' > ' + gene_seqs_aligned_fn
os.system(cmd_string)
