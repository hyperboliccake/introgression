# Given a list of introgressed regions in this format:
# Sigma1278b.chrX, + strand, regionStart-regionEnd, blockStart, blockEnd
#
# generate a set of annotations in the ../../results/analyze/ folder:
#
# for each region called introgressed, generate a file in
# results/analyze/regions/ that is the alignment for the two references and
# the strain with the introgressed region,
# S288c_CBS432_strain_chr_start-end.fasta. In this file, the aligned
# bases within coding sequence are upper case. In addition, there is a
# corresponding file S288c_CBS432_strain_chr_start-end.genes.txt
# listing the genes that overlap with this region, and the indices of
# the bases they overlap, in this format: 
# gene_name, 0-149, 25236-25385 
# gene_name, 200-600, ....
# 
# also generate a file in results/gene_alignments/ for each
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
from process import *
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../sim/')
import sim_analyze_hmm_bw as sim
sys.path.insert(0, '../misc/')
import read_maf


#####
# set up stuff
#####

tag, topology, species_to, species_from1, species_from2, \
    num_samples_species_to, num_samples_species_from1, num_samples_species_from2, \
    N0_species_to, N0_species_from1, N0_species_from2, \
    migration_from1, migration_from2, \
    expected_tract_lengths, \
    expected_num_tracts, \
    has_ref_from1, has_ref_from2, \
    rho, outcross_rate, theta, num_sites, num_reps = \
    sim.process_args(sys.argv)



# reference names for actual species
states = []
i = -1
if species_from2 != None:
    states = [sys.argv[i]]
    i -= 1
states = [sys.argv[i]] + states
i -= 1
states = [sys.argv[i]] + states

# is one of the states actually unknown (and thus not a reference)?
unknown_state = False
if species_from2 != None:
    if not has_ref_from2:
        unknown_state = True
elif not has_ref_from1:
    unknown_state = True

refs = copy.deepcopy(states)
if unknown_state:
    refs = refs[:-1]
master_ref = refs[0]


# TODO don't hard code this obvs
# shortened names
#species_to_ref_name = []
#species_to_ref_name[species_to] = refs[0]
#species_to_ref_name[species_from1] = refs[1]
#if len(refs) == 3:
#    species_to_ref_name[species_from2] = refs[2]
species = [species_to, species_from1, species_from2]
species_to_ref = {species_to:refs[0]}
if has_ref_from1:
    species_to_ref[species_from1] = refs[1]
if has_ref_from2:
    species_to_ref[species_from2] = refs[2]
ref_codes = [species[i][0].upper() for i in range(len(refs))]
ref_to_code = dict(zip(refs, ref_codes))
print ref_to_code


#####
# read in introgressed regions
#####

gp_dir = '../'
fn_all_regions = gp.analysis_out_dir_absolute + 'introgressed_hmm_' + tag + '.txt'
# introgressed regions keyed by strain and then chromosome
regions = read_regions(fn_all_regions)
"""
s = regions.keys()[0]
t = regions.keys()[1]
c = 'IV'
regions_abbr = {s:{c:{}}}
regions_abbr[t] = {c:{}}
regions_abbr[s][c] = regions[s][c][:10]
regions_abbr[t][c] = regions[t][c][:10]
regions = regions_abbr
"""
#####
# extract alignments for introgressed regions
#####

"""
introgressed region coords
alignment blocks
gene coords

for each strain x chrom in gene file
    (maybe make sure we hit all strain x chroms? 
    still do other parts if we don't have gene file?)
    read alignment block file, store all those blocks for now
    loop through all introgressed regions for this strain x chrom
         pull out sequences from appropriate alignment block
         write these to file
         also annotate with gene info and write to separate file
"""

fn_align_prefix = gp_dir + gp.alignments_dir
fn_region_prefix = gp.analysis_out_dir_absolute + '/' + tag + '/regions/'
for r in gp.alignment_ref_order:
    if r in refs:
        fn_align_prefix += r + '_'
        #fn_region_prefix += r + '_'

# just read genes from master reference for now, since that's how the
# introgressed regions are indexed
for chrm in gp.chrms:
    fn = gp.gb_master_dir + master_ref + '/' + \
        master_ref + '_chr' + chrm + '.gb'
    # for storing genes once they're read (or reading genes from if
    # already exists)
    fn_genes = gp.analysis_out_dir_absolute + '/' + \
        master_ref + '_chr' + chrm + '_genes.txt'

    print 'reading genes on chromosome', chrm
    genes = read_genes(fn, fn_genes)
    print 'done reading genes'

    for strain in regions:
        
        print '***', strain, chrm
        sys.stdout.flush()
        # skip this strain x chromosome if there are no introgressed
        # regions for it
        if strain not in regions or chrm not in regions[strain]:
            continue

        # read alignment blocks for this strain and chromosome
        fn_align = fn_align_prefix + \
            strain + '_chr' + chrm + gp.alignment_suffix
        alignment_blocks = read_maf.read_mugsy(fn_align, len(refs))
        
        # loop through all introgressed regions for this strain and
        # chromosome
        for entry in regions[strain][chrm]:
            #print 'entry', entry

            # find corresponding alignment block, which we've already kept
            # track of
            current_alignment_block = alignment_blocks[entry['block_label']]

            # write appropriate part of block to a file, and annotate with
            # gene info and write to another file
            fn_region_current_prefix = fn_region_prefix + entry['region_id']
            fn_region = fn_region_current_prefix + gp.alignment_suffix

            if not os.path.exists(os.path.dirname(fn_region)):
                os.makedirs(os.path.dirname(fn_region))
            fn_region_annotated = fn_region_current_prefix + '_annotated' + \
                '.txt'
            write_region_alignment(current_alignment_block, entry, genes, \
                                       strain, master_ref, refs, \
                                       ref_to_code, \
                                       species_to_ref[species_to], entry['predicted_reference'], \
                                       fn_region, fn_region_annotated, \
                                       context = 100)


######
# keep track of all introgressed genes and the regions they're
# introgressed in
#####

introgressed_genes = {} 

for strain in regions:
    for chrm in regions[strain]:
        for entry in regions[strain][chrm]:
            for gene in entry['genes']:
                x = (entry['region_id'], strain, \
                         entry['genes_introgressed_fractions'][gene], \
                         entry['number_non_gap'], \
                         entry['ref_from_count'])
                if gene not in introgressed_genes:
                    introgressed_genes[gene] = []
                introgressed_genes[gene].append(x)

#####
# write file adding region id and gene info to list of regions
#####

fn = gp.analysis_out_dir_absolute + tag + '/introgressed_hmm_' + tag + '.txt'

fn_all_regions_genes = fn_all_regions[:-4] + '_genes.txt'
f = open(fn_all_regions_genes, 'w')
f.write('strain\tchromosome\talignment_block_label\tstrand\tpredicted_reference\tregion_start\tregion_end\tnumber_non_gap_sites\tgenes\tref_from_count\n')
for strain in regions:
    for chrm in regions[strain]:
        for entry in regions[strain][chrm]:
            f.write(entry['region_id'] + '\t')
            f.write(strain + '\t' + chrm + '\t')
            f.write(entry['block_label'] + '\t' + entry['strand'] + '\t')
            f.write(entry['predicted_reference'] + '\t')
            f.write(str(entry['region_start']) + '\t' + str(entry['region_end']) + '\t')
            f.write(str(entry['number_non_gap']) + '\t')
            f.write(' '.join(entry['genes']) + '\t')
            f.write(str(entry['ref_from_count']) + '\n')
f.close()


#####
# summarize introgressed gene information in a few different ways
#####

# overall file with list of all introgressed genes and summary info;
# one file for each introgressed gene with row for each strain
# introgressed in

fn_all = gp.analysis_out_dir_absolute + tag + '/introgressed_hmm_' + tag + \
    '_genes_summary.txt'

summarize_gene_info(fn_all, introgressed_genes, tag, threshold=0)

# same summaries, but filter out all the regions with not enough
# support (< x sites that match reference for species predicted to be
# introgressed from)
fn_all_filtered = gp.analysis_out_dir_absolute + tag + '/introgressed_hmm_' + tag + \
    '_genes_summary_filtered.txt'

summarize_gene_info(fn_all_filtered, introgressed_genes, tag, threshold=8)



