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

#####
# read in introgressed regions
#####

gp_dir = '../'
fn_all_regions = gp.analysis_out_dir_absolute + 'introgressed_hmm_' + tag + '.txt'
# introgressed regions keyed by strain and then chromosome
regions = read_regions(fn_all_regions)
#s = regions.keys()[0]
#t = regions.keys()[1]
#c = 'IV'
#regions_abbr = {s:{c:{}}}
#regions_abbr[t] = {c:{}}
#regions_abbr[s][c] = regions[s][c][:10]
#regions_abbr[t][c] = regions[t][c][:10]
#regions = regions_abbr

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
                         entry['number_non_gap'])
                if gene not in introgressed_genes:
                    introgressed_genes[gene] = []
                introgressed_genes[gene].append(x)

#####
# write file adding region id and gene info to list of regions
#####

fn = gp.analysis_out_dir_absolute + tag + '/introgressed_hmm_' + tag + '.txt'

fn_all_regions_genes = fn_all_regions[:-4] + '_genes.txt'
f = open(fn_all_regions_genes, 'w')
f.write('strain\tchromosome\talignment_block_label\tstrand\tpredicted_reference\tregion_start\tregion_end\tnumber_non_gap_sites\tgenes\n')
for strain in regions:
    for chrm in regions[strain]:
        for entry in regions[strain][chrm]:
            f.write(entry['region_id'] + '\t')
            f.write(strain + '\t' + chrm + '\t')
            f.write(entry['block_label'] + '\t' + entry['strand'] + '\t')
            f.write(entry['predicted_reference'] + '\t')
            f.write(str(entry['region_start']) + '\t' + str(entry['region_end']) + '\t')
            f.write(str(entry['number_non_gap']) + '\t')
            f.write(' '.join(entry['genes']) + '\n')
f.close()


#####
# summarize introgressed gene information in a few different ways
#####

# overall file with list of all introgressed genes and summary info;
# one file for each introgressed gene with row for each strain
# introgressed in

fn_all = gp.analysis_out_dir_absolute + tag + '/introgressed_hmm_' + tag + \
    '_genes_summary.txt'
f_all = open(fn_all, 'w')

for gene in introgressed_genes:
    # keyed by strain, because gene can be broken across multiple
    # alignment blocks/regions for the same strain
    sum_introgressed_fraction = {}
    sum_number_non_gap = {}
    fn_gene = gp.analysis_out_dir_absolute + tag + '/genes/' + gene + '.txt'
    if not os.path.exists(os.path.dirname(fn_gene)):
        os.makedirs(os.path.dirname(fn_gene))
    f_gene = open(fn_gene, 'w')
    for entry in introgressed_genes[gene]:
        region_id, strain, introgressed_fraction, number_non_gap = entry
        if strain not in sum_introgressed_fraction:
            sum_introgressed_fraction[strain] = 0
            sum_number_non_gap[strain] = 0
        sum_introgressed_fraction[strain] += introgressed_fraction
        sum_number_non_gap[strain] += number_non_gap
        f_gene.write(region_id + '\t' + strain + '\t' + \
                         str(introgressed_fraction) + '\t' + str(number_non_gap) + '\n')

    f_gene.close()

    # now do averaging over strains
    num_strains = len(sum_introgressed_fraction)
    avg_introgressed_fraction = sum(sum_introgressed_fraction.values()) / \
        float(num_strains)
    avg_number_non_gap = sum(sum_number_non_gap.values()) / \
        float(num_strains)

    avg_introgressed_fraction /= float(num_strains)
    avg_number_non_gap /= float(num_strains)
    f_all.write(gene + '\t' + str(num_strains) + '\t' + \
                    str(avg_introgressed_fraction) + '\t' + \
                            str(avg_number_non_gap) + '\n')
f_all.close()


"""
f = open(gp.gb_all, 'r')


line = f.readline()
result = read_one_strain_chrm(f, line)
c = 0
while result != None:
    # current strain, chromosome, and set of genes
    genes, strain, chrm, line = result
    print strain, chrm

    # get next set of genes for different strain and chromosome
    result = read_one_strain_chrm(f, line)

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
        # find corresponding alignment block, which we've already kept
        # track of
        current_alignment_block = alignment_blocks[entry['block_label']]

        # write appropriate part of block to a file, and annotate with
        # gene info and write to another file
        fn_region_current_prefix = fn_region_prefix + \
            strain + '_chr' + chrm + '_' + \
            str(entry['region_start']) + '-' + str(entry['region_end'])
        fn_region = fn_region_current_prefix + gp.alignment_suffix
        fn_region_annotated = fn_region_current_prefix + '_annotated' + \
            gp.alignment_suffix
        write_region_alignment(current_alignment_block, entry, genes, \
                                   strain, master_ref, refs, \
                                   fn_region, fn_region_annotated, \
                                   context = 100)
        print fn_region
    sys.exit()
"""

"""
##### 
# For each introgressed region, extract relevant part of alignment to a separate file
#####

fn_align_prefix = gp_dir + gp.alignments_dir
fn_region_prefix = gp_dir + gp.regions_out_dir
for r in gp.alignment_ref_order:
    if r in refs:
        fn_align_prefix += r + '_'
        fn_region_prefix += r + '_'
        

alignment_blocks = {}
for strain in regions:
    alignment_blocks[strain] = {}
    print strain
    for chrm in regions[strain]:
        print '   ', chrm
        
        fn_align = fn_align_prefix + \
            strain + '_chr' + chrm + gp.alignment_suffix

        alignment_blocks[strain][chrm] = read_maf.read_mugsy(fn_align, len(refs))

        for region in regions[strain][chrm]:
            fn_region = fn_align_prefix

            alignment_block = alignment_blocks[strain][chrm][region['block_label']]

            fn_region = fn_region_prefix + \
                strain + '_chr' + chrm + '_' + \
                str(region['region_start']) + '-' + str(region['region_end']) + \
                gp.alignment_suffix

            write_region_alignment(alignment_block, region, strain, \
                                       master_ref, refs, fn_region)


#####
# Identify genes (if any) that overlap each alignment block and
# annotate the alignment file for each region to note where the
# introgressed boundaries and genes are
#####

# you might think that we could just read in the genes from the master
# reference sequence, but really we want to do it for all the
# non-reference strains because gaps would be ambiguous otherwise

f = open(gp.gb_all, 'r')

line = f.readline()
result = read_one_strain_chrm(f, line)
while result != None:
    genes, strain, chrm, line = result
    
    result = read_one_strain_chrm(f, line)

    # now that we have all the genes, proceed with annotations (if
    # there's anything to annotate)

    # skip annotating if we aren't interested in this strain x
    # chromosome
    if strain not in regions or chrm not in regions[strain]:
        continue

    # otherwise, loop through all the introgressed regions for this
    # strain x chromosome and see whether any genes fall within them
    for i in range(len(regions[strain][chrm])):
        entry = regions[strain][chrm][i]
        fn_region_annotated = fn_region_prefix + \
            strain + '_chr' + chrm + '_' + \
            str(entry['start']) + '-' + str(entry['end']) + \
            '_annotated' + gp.alignment_suffix

        entry = write_annotated_region_alignment(\
            alignment_blocks[strain][chrm][entry['block_label']], \
                strain, entry, fn_region_annotated, genes, \
                introgressed_genes_strains, strains_introgressed_genes)
        regions[strain][chrm][i] = entry
"""
