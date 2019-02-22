import sys
import os
import gzip
from count_coding_changes import *
import annotate_positions
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc/')
import overlap
import read_table
import read_fasta

##======
# command line arguments
##======

tag = sys.argv[1]

##======
# read in introgressed regions
##======

# key region ids by chromosome and then strain 
fn_regions = gp.analysis_out_dir_absolute + tag + '/' + \
             'introgressed_blocks_filtered_par_' + tag + '_summary_plus.txt'
regions, l = read_table.read_table_rows(fn_regions, '\t')
region_ids_by_chrm_strain = {}
for r in regions.keys():
    strain = regions[r]['strain']
    chrm = regions[r]['chromosome']
    if not region_ids_by_chrm_strain.has_key(chrm):
        region_ids_by_chrm_strain[chrm] = {}
    if not region_ids_by_chrm_strain[chrm].has_key(strain):
        region_ids_by_chrm_strain[chrm][strain] = []
    region_ids_by_chrm_strain[chrm][strain].append(r)

# also get genes they overlap
fn = gp.analysis_out_dir_absolute + tag + '/' + \
     'genes_for_each_region_' + tag + '.txt'
region_genes = {}
f = open(fn, 'r')
line = f.readline()
while line != '':
    line = line[:-1].split('\t')
    if line[0] in regions:
        region_genes[line[0]] = line[2::2]
    line = f.readline()
f.close()


##======
# count sites within all regions that are coding/noncoding, plus some
# more details about coding changes
##======

other_ref = gp.alignment_ref_order[1]

region_totals = {}
gene_totals = {}
strain_totals = {}
totals = {'syn':0, 'non':0, 'syn_ref':0, 'non_ref':0, \
          'insert':0, 'delete':0, 'insert_ref':0, 'delete_ref':0, \
          'gene_delete':0, 'gene_delete_ref':0, \
          'ref_gene_only':0, 'strain_orf_only':0, \
          'coding':0, 'noncoding':0, 'frameshift':0}

for chrm in gp.chrms:

    print chrm

    # read in cer reference genes
    fn = gp.analysis_out_dir_absolute + gp.master_ref + '_chr' + chrm + \
         '_genes.txt'
    genes, l = read_table.read_table_rows(fn, '\t', header=False, key_ind=0)
    for gene in genes:
        genes[gene] = (int(genes[gene][0]), int(genes[gene][1]))

    # read in cer ref -> par ref position file
    fn = gp.analysis_out_dir_absolute + 'coordinates/' + gp.master_ref + \
         '_to_' + other_ref + '_chr' + chrm + '.txt.gz'
    master_to_other_ref_pos = [float(line[:-1]) \
                               for line in gzip.open(fn, 'rb').readlines()]

    # read in cer ref chromosome sequence
    fn = gp.ref_dir[gp.master_ref] + gp.ref_fn_prefix[gp.master_ref] + \
         '_chr' + chrm + gp.fasta_suffix
    master_seq = read_fasta.read_fasta(fn)[1][0]

    # read in par ref chromosome sequence
    fn = gp.ref_dir[other_ref] + gp.ref_fn_prefix[other_ref] + \
         '_chr' + chrm + gp.fasta_suffix
    other_ref_seq = read_fasta.read_fasta(fn)[1][0]

    # read in par ref ORFs
    fn = gp.ref_dir[other_ref] + 'orfs/' + other_ref + \
         '_chr' + chrm + '_orfs' + gp.fasta_suffix
    ref_orfs = annotate_positions.get_orfs(fn)

    for strain in region_ids_by_chrm_strain[chrm].keys():
        print '-', strain

        if not strain_totals.has_key(strain):
            strain_totals[strain] = {'syn':0, 'non':0, 'syn_ref':0, 'non_ref':0, \
                                     'ref_gene_only':0, 'strain_orf_only':0, \
                                     'coding':0, 'noncoding':0}
        
        # read in cer ref -> strain position file 
        fn = gp.analysis_out_dir_absolute + 'coordinates/' + gp.master_ref + \
             '_to_' + strain + '_chr' + chrm + '.txt.gz'
        master_to_strain_pos = [float(line[:-1]) \
                                for line in gzip.open(fn, 'rb').readlines()]

        # read in strain chromosome sequence
        fn = gp.non_ref_dirs[gp.master_ref][0] + strain + \
             '_chr' + chrm + gp.fasta_suffix
        strain_seq = read_fasta.read_fasta(fn)[1][0]

        # read in strain ORFs
        fn = gp.non_ref_dirs[gp.master_ref][0] + 'orfs/' + strain + \
             '_chr' + chrm + '_orfs' + gp.fasta_suffix
        orfs = annotate_positions.get_orfs(fn)

        for region in region_ids_by_chrm_strain[chrm][strain]:
            region_totals[region] = {'syn':0, 'non':0, 'syn_ref':0, 'non_ref':0, \
                                     'ref_gene_only':0, 'strain_orf_only':0, \
                                     'coding':0, 'noncoding':0}

            # is each site in region in a master ref gene and/or
            # strain ORF?
            # TODO the same thing but for variants (that match par ref)
            t_gene_orf = 0
            t_gene_not_orf = 0
            t_not_gene_orf = 0
            t_not_gene_not_orf = 0
            for site in range(int(regions[region]['start']), \
                              int(regions[region]['end'])):
                in_gene = overlap.contained_any(site, genes.values())
                in_orf = overlap.contained_any(master_to_strain_pos[site], orfs.keys())
                if in_gene:
                    if in_orf:
                        t_gene_orf += 1
                    else:
                        t_gene_not_orf += 1
                else:
                    if in_orf:
                        t_not_gene_orf += 1
                    else:
                        t_not_gene_not_orf += 1

            totals['ref_gene_only'] += t_gene_not_orf
            totals['strain_orf_only'] += t_not_gene_orf
            totals['coding'] += t_gene_orf
            totals['noncoding'] += t_not_gene_not_orf

            strain_totals[strain]['ref_gene_only'] += t_gene_not_orf
            strain_totals[strain]['strain_orf_only'] += t_not_gene_orf
            strain_totals[strain]['coding'] += t_gene_orf
            strain_totals[strain]['noncoding'] += t_not_gene_not_orf

            region_totals[region]['ref_gene_only'] += t_gene_not_orf
            region_totals[region]['strain_orf_only'] += t_not_gene_orf
            region_totals[region]['coding'] += t_gene_orf
            region_totals[region]['noncoding'] += t_not_gene_not_orf

            # now evalate individual genes
            for gene in region_genes[region]:

                gene_start = genes[gene][0]
                gene_end = genes[gene][1]

                # read multiple alignment for the gene, in which we've
                # previously selected the best orfs to match the gene
                fn = gp.analysis_out_dir_absolute + tag + '/genes/' + gene + '/' + \
                     gene + '_introgressed_filtered.maf'
                if not os.path.isfile(fn):
                    print 'do not have alignment for', gene
                    continue
                aligned_genes = get_aligned_genes(fn, \
                                                  [gp.master_ref, other_ref, strain])

                print gene, strain
                # for now, ignore cerevisiae reference genes that
                # don't map perfectly to an ORF in the strain and
                # paradoxus reference
                #if ambiguous(gene, gene_start, gene_end, master_to_strain_pos, orfs):
                #    continue
                #if ambiguous(gene, gene_start, gene_end, \
                #             master_to_other_ref_pos, ref_orfs):
                #    continue
                
                # extract gene sequence from references and strain
                g_master = master_seq[gene_start:gene_end+1]
                g_ref = other_ref_seq[int(master_to_other_ref_pos[gene_start]):\
                                      int(master_to_other_ref_pos[gene_end])+1]
                g_strain = strain_seq[int(master_to_strain_pos[gene_start]):\
                                      int(master_to_strain_pos[gene_end])+1]

                # get overlap between gene and introgressed region
                o_start, o_end = overlap.overlap_region(genes[gene][0], \
                                                        genes[gene][1], \
                                                        int(regions[region]['start']), \
                                                        int(regions[region]['end']))

                # count synonymous and non synonymous changes due to
                # paradoxus (deal with gene direction correctly)
                # t_syn, t_non = count_coding(g_master, g_ref, g_strain, \
                #                             o_start-gene_start, o_end-gene_start)

                # alternative method that deals with imperfect matches
                t_syn, t_non, t_syn_ref, t_non_ref, \
                    t_insert, t_delete, t_insert_ref, t_delete_ref, \
                    gene_delete, gene_delete_ref, frameshift = \
                    count_coding_with_gaps(aligned_genes[gp.master_ref], \
                                           aligned_genes[other_ref], \
                                           aligned_genes[strain], \
                                           o_start-gene_start, o_end-gene_start)

                # add to totals for region, gene, strain, and overall
                if not gene_totals.has_key(gene):
                    gene_totals[gene] = {'syn':0, 'non':0, 'syn_ref':0, 'non_ref':0, \
                                         'insert':0, 'delete':0, \
                                         'insert_ref':0, 'delete_ref':0, \
                                         'gene_delete':0, 'gene_delete_ref':0, \
                                         'frameshift':0}
                gene_totals[gene]['syn'] += t_syn
                gene_totals[gene]['non'] += t_non
                gene_totals[gene]['syn_ref'] += t_syn_ref
                gene_totals[gene]['non_ref'] += t_non_ref
                gene_totals[gene]['insert'] += t_insert
                gene_totals[gene]['delete'] += t_delete
                gene_totals[gene]['insert_ref'] += t_insert_ref
                gene_totals[gene]['delete_ref'] += t_delete_ref
                gene_totals[gene]['gene_delete'] += gene_delete
                gene_totals[gene]['gene_delete_ref'] += gene_delete_ref
                gene_totals[gene]['frameshift'] += frameshift

                region_totals[region]['syn'] += t_syn
                region_totals[region]['non'] += t_non
                region_totals[region]['syn_ref'] += t_syn_ref
                region_totals[region]['non_ref'] += t_non_ref

                strain_totals[strain]['syn'] += t_syn
                strain_totals[strain]['non'] += t_non
                strain_totals[strain]['syn_ref'] += t_syn_ref
                strain_totals[strain]['non_ref'] += t_non_ref

                totals['syn'] += t_syn
                totals['non'] += t_non
                totals['syn_ref'] += t_syn_ref
                totals['non_ref'] += t_non_ref
                totals['insert'] += t_insert
                totals['delete'] += t_delete
                totals['insert_ref'] += t_insert_ref
                totals['delete_ref'] += t_delete_ref
                totals['gene_delete'] += gene_delete
                totals['gene_delete_ref'] += gene_delete_ref
                totals['frameshift'] += frameshift

##======
# write output file
##======

fn = gp.analysis_out_dir_absolute + tag + '/' + 'coding_changes_summary_' + \
     tag + '.txt'
f = open(fn, 'w')
sep = '\t'
f.write('\t'.join(['label', 'label_type', 'total', 'total_type']) + '\n')

for key in totals.keys():
    f.write('all' + sep + 'all' + sep + str(totals[key]) + sep + key + '\n')

for strain in strain_totals:
    for key in strain_totals[strain].keys():
        f.write(strain + sep + 'strain' + sep + \
                str(strain_totals[strain][key]) + sep + key + '\n')

for gene in gene_totals:
    for key in gene_totals[gene].keys():
        f.write(gene + sep + 'gene' + sep + \
                str(gene_totals[gene][key]) + sep + key + '\n')

for region in region_totals:
    for key in region_totals[region].keys():
        f.write(region + sep + 'region' + sep + \
                str(region_totals[region][key]) + sep + key + '\n')

f.close()    




# new plan
# for each region

# for each site in region
# is it in ref gene and/or strain orf? (keep track of four totals)
# 

# for each gene
# get corresponding orfs in par and strain
# align
# get lengths
# count syn and non changes in strain due to par
# same as before but deal with gaps
# - categories:
#   multiples of 3
#   not multiples of 3 -> stop counting/ignore gene?

