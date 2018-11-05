# get set of genes that are well-aligned (generally about the same
# length in all the strains, etc) and concatenate all the alignments
# in one big fasta file; also separate out all the genes that
# have/don't have any amount of intogression and put those alignments
# in separate files; also create files with all gaps removed

import sys
import os
from get_gene_set import *
import math
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../align/')
import align_helpers
sys.path.insert(0, '../analyze/')
import gene_predictions
import combine_all_strains
sys.path.insert(0, '../misc/')
import read_table
import read_fasta
import write_fasta
import mystats

tag = 'u3_i.001_tv_l1000_f.01'
suffix = '_filtered'

threshold = .1

# look up all cerevisiae reference genes
print 'getting reference gene coordinates'
all_ref_genes = {}
for chrm in gp.chrms:
    genes_fn = gp.analysis_out_dir_absolute + gp.master_ref + \
               '_chr' + chrm + '_genes.txt'
    all_ref_genes[chrm] = {}
    for line in open(genes_fn, 'r').readlines():
        gene, start, end = line.split('\t')
        all_ref_genes[chrm][gene] = (int(start), int(end))

print 'getting reference gene sequences'
# files for storing all gene sequences, unaligned and unaligned
genes_dir = gp.analysis_out_dir_absolute + tag + '/genes/'
for chrm in gp.chrms:
    for gene in all_ref_genes[chrm]:
        seq = read_fasta.read_fasta(gp.ref_dir[gp.master_ref] + \
                                    gp.ref_fn_prefix[gp.master_ref] + \
                                    '_chr' + chrm + gp.fasta_suffix)[1][0]
        start, end = all_ref_genes[chrm][gene]
        seq = seq[start:end + 1]
        if not os.path.isdir(genes_dir + gene):
            os.makedirs(genes_dir + gene)
        f = open(genes_dir + gene + '/' + gene + \
                 '_from_alignment' + gp.fasta_suffix, 'w')
        f.write('> ' + gp.master_ref + '\n')
        f.write(seq + '\n')
        f.close()

gp_dir = '../'
s = align_helpers.get_strains(align_helpers.flatten(gp.non_ref_dirs.values()))


print 'getting introgressed genes'
# read in filtered regions
fn_regions = gp.analysis_out_dir_absolute + tag + '/' + \
             'introgressed_blocks' + suffix + '_par_' + tag + '_summary_plus.txt'
regions, l = read_table.read_table_rows(fn_regions, '\t')

# read in genes for each region
fn_genes_for_each_region = gp.analysis_out_dir_absolute + tag + '/' + \
                           'genes_for_each_region_' + tag + '.txt'
introgressed_genes = set()
for line in open(fn_genes_for_each_region, 'r').readlines():
    line = line.split('\t')
    region = line[0]
    if region in regions:
        g = line[2::2]
        for gene in g:
            introgressed_genes.add(gene)


# loop through all genes by chromosome
print 'getting gene sequences for all other strains and aligning them'
genes = {}
genes_nonint = []
genes_int = []
genes_unflagged = []

flagged = set()

for chrm in gp.chrms:
    print chrm

    # get gene sequences for this strain and all genes 
    
    # for each strain, get the corresponding part of chromosome
    # sequence (don't worry about ORF); assess alignment quality, flag
    # if it doesn't meet a threshold
    
    # get indexing from reference genome to strain genomes based on
    # alignments

    # for par reference which doesn't have site summary file
    align_fn = gp_dir + gp.alignments_dir + \
               '_'.join(gp.alignment_ref_order) + '_chr' + chrm + \
               '_mafft' + gp.alignment_suffix
    t = combine_all_strains.get_inds_from_alignment(align_fn, False)
    other_ref_strain = gp.ref_fn_prefix[gp.alignment_ref_order[1]]
    ref_ind_to_strain_ind = dict(zip([int(x) for x in t['ps_ref']], \
                                     [float(x) for x in t['ps_strain']]))
    chrom_seq = read_fasta.read_fasta(gp.ref_dir[other_ref_strain] + \
                                      gp.ref_fn_prefix[other_ref_strain] + \
                                      '_chr' + chrm + gp.fasta_suffix)[1][0]
    for gene in all_ref_genes[chrm].keys():
        ref_start, ref_end = all_ref_genes[chrm][gene]

        start = int(ref_ind_to_strain_ind[ref_start])
        end = int(ref_ind_to_strain_ind[ref_end])
        seq = chrom_seq[start:end + 1]
        f = open(genes_dir + gene + '/' + gene + \
                 '_from_alignment' + gp.fasta_suffix, 'a')
        f.write('> ' + other_ref_strain + '\n')
        f.write(seq + '\n')
        f.close()

    for strain, d in s:
        print '*', strain
        sys.stdout.flush()

        align_fn = gp_dir + gp.alignments_dir + \
                   '_'.join(gp.alignment_ref_order) + '_' + \
                   strain + '_chr' + chrm + \
                   '_mafft' + gp.alignment_suffix
        t = combine_all_strains.get_inds_from_alignment(align_fn, False, 0, 2)
        ref_ind_to_strain_ind = dict(zip([int(x) for x in t['ps_ref']], \
                                         [float(x) for x in t['ps_strain']]))
        chrom_seq = read_fasta.read_fasta(d + '/' + strain + \
                                    '_chr' + chrm + gp.fasta_suffix)[1][0]

        for gene in all_ref_genes[chrm].keys():
            ref_start, ref_end = all_ref_genes[chrm][gene]
            start = int(ref_ind_to_strain_ind[float(ref_start)])
            end = int(ref_ind_to_strain_ind[float(ref_end)])
            seq = chrom_seq[start:end + 1]
            ref_length = ref_end - ref_start + 1
            if len(seq) > (1 + threshold) * ref_length or \
               len(seq) < (1 - threshold) * ref_length:
                flagged.add(gene)
            f = open(genes_dir + gene + '/' + gene + \
                     '_from_alignment' + gp.fasta_suffix, 'a')
            f.write('> ' + strain + '\n')
            f.write(seq + '\n')
            f.close()

    # if we haven't failed to meet alignment quality threshold in any
    # strains, append this gene to the current lists;
    # check whether gene overlaps introgressed region in any strains, keep
    # separate lists
    for gene in all_ref_genes[chrm].keys():
        
        gene_seqs_fn = genes_dir + gene + '/' + gene + \
                       '_from_alignment' + gp.fasta_suffix
        gene_seqs_aligned_fn = gene_seqs_fn.replace(gp.fasta_suffix, \
                                                    gp.alignment_suffix)
        # make alignments for _all_ the genes
        cmd_string = gp.mafft_install_path + '/mafft ' + \
                     ' --quiet --reorder --preservecase ' + \
                     gene_seqs_fn + ' > ' + gene_seqs_aligned_fn
        os.system(cmd_string)

        if gene not in flagged:
            genes_unflagged.append(gene)
            if gene in introgressed_genes:
                genes_int.append(gene)
            else:
                genes_nonint.append(gene)

        
# then just for unflagged, int, and nonint, append alignments so
# there's one master file with big entry for each strain
concatenate_fastas([genes_dir + gene + '/' + gene + '_from_alignment' + gp.alignment_suffix for gene in genes_unflagged], gp.analysis_out_dir_absolute + tag + '/genes_for_phylogeny_all' + gp.fasta_suffix, True)
f = open(gp.analysis_out_dir_absolute + tag + '/genes_for_phylogeny_all.txt', 'w')
for gene in genes_unflagged:
    f.write(gene + '\n')
f.close()

concatenate_fastas([genes_dir + gene + '/' + gene + '_from_alignment' + gp.alignment_suffix for gene in genes_int], gp.analysis_out_dir_absolute + tag + '/genes_for_phylogeny_int' + gp.fasta_suffix, True)
f = open(gp.analysis_out_dir_absolute + tag + '/genes_for_phylogeny_int.txt', 'w')
for gene in genes_int:
    f.write(gene + '\n')
f.close()

concatenate_fastas([genes_dir + gene + '/' + gene + '_from_alignment' + gp.alignment_suffix for gene in genes_nonint], gp.analysis_out_dir_absolute + tag + '/genes_for_phylogeny_nonint' + gp.fasta_suffix, True)
f = open(gp.analysis_out_dir_absolute + tag + '/genes_for_phylogeny_nonint.txt', 'w')
for gene in genes_nonint:
    f.write(gene + '\n')
f.close()


# then build phylogeny!




