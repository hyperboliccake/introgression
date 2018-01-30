# Given a list of introgresssed regions in this format:
# strain\tchromosome\tpredicted_species\tstart\tend\tnumber_non_gap
#
# generate a set of annotations in the ../../results/analyze/tag/ folder:
#
# for each region called introgressed, generate a file in
# results/analyze/tag/regions/ that is the alignment for the two references and
# the strain with the introgressed region,
# S288c_CBS432_strain_chrX_start-end.fasta. In this file, the aligned
# bases within coding sequence are upper case. In addition, there is a
# corresponding file S288c_CBS432_strain_chrX_start-end.genes.txt
# listing the genes that overlap this region, and the indices of
# the bases they overlap, in this format: 
# gene_name\t0-149\t25236-25385 
# gene_name\t200-600\t....
# 
# also generate a file in results/tag/gene_alignments/ for each
# introgressed gene, which contains one threeway alignment for each
# strain in which the gene was called introgressed...followed by all
# the genes that weren't
#

# todo in future:
# for each gene that is called introgressed in at least one strain,
# create folder gene/ containing the alignment
# of cerevisiae and paradoxus references to all the introgressed
# versions (gene_introgressed.fasta), and also to all of the versions
# (gene_all.fasta).

# TODO: 
## _annotated file should be .txt not .maf
## also modify so that 80 characters per line
## and extra row showing summary of which references match


import re
import sys
import os
import copy
from gene_predictions import *
import predict
import pickle
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc/')
import read_fasta


sys.exit()

##======
# read in analysis parameters
##======

refs, strains, args = predict.process_args(sys.argv)

resume = False
open_mode = 'w'
if resume:
    open_mode = 'a'

##======
# read in introgressed regions
##======

all_species_from = args['species'][1:] + ['unknown'] #TODO

gp_dir = '../'
all_blocks_fn = {}
for species_from in all_species_from:
    blocks_fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
                'introgressed_blocks_' + species_from + '_' + args['tag'] + '.txt'
    all_blocks_fn[species_from] = blocks_fn
# introgressed regions keyed by strain and then chromosome: 
# (start, end, number_non_gap, region_id)
# start and end indices are relative to (unaligned) master ref sequence
all_regions = {}
for species_from in all_species_from:
    regions = predict.read_blocks(all_blocks_fn[species_from], region_id=True)
    all_regions[species_from] = regions
regions = all_regions['par'] # TODO
species_from = 'par'

##======
# extract alignments and genes for introgressed regions
##======


# produce region summmary file with all the same info, but also with
# region ids (r1-rn), and with genes overlapping each region

# produce file with one line per introgressed gene, list of all
# strains with introgressed regions overlapping that gene

# one file for each region with alignment (fasta)

# above but with surrounding context and some extra annotation lines:
# part that's actually called introgressed, part that's within
# gene(s), part that matches reference 1, reference 2, etc; also
# incolude list of the genes that are shown in that region, in order
# they appear

# for alignment files (input)
fn_align_prefix = gp_dir + gp.alignments_dir
fn_align_prefix += '_'.join([refs[s][0] for s in args['species']]) + '_'

# for annotated region files (output)
fn_region_prefix = gp.analysis_out_dir_absolute + '/' + args['tag'] + '/regions/'
if not os.path.isdir(fn_region_prefix):
    os.makedirs(fn_region_prefix)

# summary files
fn_region_summary = gp.analysis_out_dir_absolute + '/' + args['tag'] + '/' + \
                   'introgressed_blocks_' + species_from + '_' + args['tag'] + \
                   '_summary.txt'
f_region_summary = open(fn_region_summary, open_mode)
# TODO fix this so that refs is ordered (instead of dic) and in correct order!
refs_ordered = args['species']
write_region_summary_header(refs_ordered, f_region_summary)

fn_genes_regions = gp.analysis_out_dir_absolute + '/' + args['tag'] + '/' + \
                  'genes_for_each_region_' + args['tag'] + '.txt'
f_genes_regions = open(fn_genes_regions, open_mode)

fn_regions_strains = gp.analysis_out_dir_absolute + '/' + args['tag'] + '/' + \
                  'regions_for_each_strain_' + args['tag'] + '.txt'
f_regions_strains = open(fn_regions_strains, open_mode)

fn_genes_strains = gp.analysis_out_dir_absolute + '/' + args['tag'] + '/' + \
                  'genes_for_each_strain_' + args['tag'] + '.txt'
f_genes_strains = open(fn_genes_strains, open_mode)

fn_strains_genes = gp.analysis_out_dir_absolute + '/' + args['tag'] + '/' + \
                  'strains_for_each_gene_' + args['tag'] + '.txt'
f_strains_genes = open(fn_strains_genes, open_mode)

# for keeping track of all genes introgressed in each strain, and the
# fraction introgressed
# keyed by strain, then gene, total fraction introgressed
# TODO maybe don't just pick par (and exclude unknown)
strain_genes_dic = dict(zip(regions.keys(), [{} for strain in regions.keys()]))
# the inverse of above
gene_strains_dic = {}

chrms_completed = []

if resume:
    try:
        f_strain_genes_dic = open(gp.analysis_out_dir_absolute + '/' + args['tag'] + \
                                  '/' + 'strain_genes_dic.pkl', 'rb')
        strain_genes_dic = pickle.load(f_strain_genes_dic)
        f_strain_genes_dic.close()

        f_gene_strains_dic = open(gp.analysis_out_dir_absolute + '/' + args['tag'] + \
                                  '/' + 'gene_strains_dic.pkl', 'rb')
        gene_strains_dic = pickle.load(f_gene_strains_dic)
        f_gene_strains_dic.close()

        f_chrms_completed = open(gp.analysis_out_dir_absolute + '/' + args['tag'] + \
                                  '/' + 'chrms_completed.pkl', 'rb')
        chrms_completed =  pickle.load(f_chrms_completed)
        f_chrms_completed.close()
    except:
        pass

# just read genes from master reference (species[0]),
# since that's how the introgressed regions are indexed
master_ref = refs[args['species'][0]][0]

ref_labels = [refs[s][0] for s in args['species']]

# deal with each chromosome separately because memory
for chrm in gp.chrms:

    if chrm in chrms_completed:
        continue

    # genbank file to read from
    fn = gp.ref_gb_dir[master_ref] + master_ref + '_chr' + chrm + '.gb'

    # for storing genes once they're read (or reading genes from if
    # already exists)
    fn_genes = gp.analysis_out_dir_absolute + '/' + \
        master_ref + '_chr' + chrm + '_genes.txt'

    print 'reading genes on chromosome', chrm
    # dictionary keyed by name: (start, end)
    genes = read_genes(fn, fn_genes)
    print 'done reading genes'

    # loop through all strains that we've called introgression in, and
    # associate genes with the regions they overlap
    for strain in regions.keys():
        
        print '***', strain, chrm
        sys.stdout.flush()
        # skip this strain x chromosome if there are no introgressed
        # regions for it
        if strain not in regions or chrm not in regions[strain]:
            continue

        # read alignment blocks for this strain and chromosome
        fn_align = fn_align_prefix + \
            strain + '_chr' + chrm + '_mafft' +  gp.alignment_suffix
        alignment_headers, alignment_seqs = read_fasta.read_fasta(fn_align)

        labels = ref_labels + [strain]
        
        # mark each site as matching each reference or not
        ref_match_by_site = get_ref_match_by_site(alignment_seqs, labels)
        # mark each site as in a gene or not
        genes_by_site = get_genes_by_site(genes, alignment_seqs[0])
        # mark each site as introgressed or not
        introgressed_by_site = \
            get_introgressed_by_site(regions[strain][chrm], alignment_seqs[0])

        #====
        # region maf and annotated alignment files
        #====

        # loop through all introgressed regions for this strain and
        # chromosome and pull out appropriate parts of multiple
        # alignment and annotate with useful info
        for entry in regions[strain][chrm]:
            # write appropriate part of alignment to a file, and annotate with
            # gene info and write to another file
            fn_region_current_prefix = fn_region_prefix + entry[3]
            fn_region = fn_region_current_prefix + gp.alignment_suffix
            if not os.path.exists(os.path.dirname(fn_region)):
                os.makedirs(os.path.dirname(fn_region))

            # regions are indexed by (unaligned) master ref sequence
            write_region_alignment(alignment_headers, alignment_seqs, fn_region, \
                                   entry[0], entry[1], 0)
            

            # write region to file in annotated/readable format
            fn_region_annotated = fn_region_current_prefix + '_annotated' + \
                '.txt'
            relative_start, relative_end = \
                write_region_alignment_annotated(labels, alignment_seqs, \
                                                 fn_region_annotated, \
                                                 entry[0], entry[1], 0, \
                                                 genes, ref_match_by_site, 
                                                 genes_by_site, \
                                                 introgressed_by_site, 100)

            #====
            # region summary file with extra info
            #====

            # strain chromosome predicted_species start end number_non_gap
            # number_match_ref1 number_match_ref2 number_match_only_ref1
            # number_match_ref2_not_ref1 number_mismatch_all_ref

            write_region_summary_line(entry, strain, chrm, species_from, \
                                      alignment_seqs, labels, \
                                      relative_start, relative_end, f_region_summary)

            #====
            # genes for each region summary file
            #====

            # region_id num_genes gene frac_intd gene frac_intd

            frac_intd = write_genes_for_each_region_summary_line(entry[3], \
                                                                 genes_by_site, \
                                                                 genes, \
                                                                 relative_start, \
                                                                 relative_end, \
                                                                 alignment_seqs[0], \
                                                                 f_genes_regions)
            for gene in frac_intd:
                # keep track of all genes for each strain...
                if not strain_genes_dic[strain].has_key(gene):
                    strain_genes_dic[strain][gene] = 0
                strain_genes_dic[strain][gene] += frac_intd[gene]
                # and all strains for each gene
                if not gene_strains_dic.has_key(gene):
                    gene_strains_dic[gene] = {}
                if not gene_strains_dic[gene].has_key(strain):
                    gene_strains_dic[gene][strain] = 0
                gene_strains_dic[gene][strain] += frac_intd[gene]

    
    f_strain_genes_dic = open(gp.analysis_out_dir_absolute + '/' + args['tag'] + \
                              '/' + 'strain_genes_dic.pkl', 'wb')
    pickle.dump(strain_genes_dic, f_strain_genes_dic)
    f_strain_genes_dic.close()

    f_gene_strains_dic = open(gp.analysis_out_dir_absolute + '/' + args['tag'] + \
                              '/' + 'gene_strains_dic.pkl', 'wb')
    pickle.dump(gene_strains_dic, f_gene_strains_dic)
    f_gene_strains_dic.close()

    chrms_completed.append(chrm)
    f_chrms_completed = open(gp.analysis_out_dir_absolute + '/' + args['tag'] + \
                             '/' + 'chrms_completed.pkl', 'wb')
    pickle.dump(chrms_completed, f_chrms_completed)
    f_chrms_completed.close()

    print 'saving intermediate results completed successfully'
        
#====
# strains for each gene summary file
#====

# (could do this for one chromsoome at a time if we wanted)    
# gene num_strains strain frac_intd strain frac_intd
    
write_strains_for_each_gene_lines(gene_strains_dic, f_strains_genes)

#====
# regions for each strain summary file
#====

# strain num_regions region length region length

write_regions_for_each_strain(regions, f_regions_strains)

#====
# genes for each strain summary file
#====

# strain num_genes gene frac_intd gene frac_intd

write_genes_for_each_strain(strain_genes_dic, f_genes_strains)


f_region_summary.close()
f_genes_regions.close()
f_regions_strains.close()
f_genes_strains.close()
f_strains_genes.close()
        

            
