import math
import numpy.random
import sys
import process_helpers
from summary_stats_helpers import *
sys.path.insert(0, '../misc/')
import mystats
sys.path.insert(0, '../')
import global_params as gp


# remove duplicates, aka why the fuck are chromsome XIV entries
# duplicated??
def remove_duplicates(regions):
    r = {}
    for strain in regions:
        r[strain] = {}
        for chrm in regions[strain]:
            entries = []
            for entry in regions[strain][chrm]:
                if entry not in entries:
                    entries.append(entry)
            r[strain][chrm] = entries
    return r

# distribution of gaps between adjacent introgressed sequeneces -->
# want to get an idea of whether we're not stitching nearby ones
# together
def get_gap_hist(regions):
    gap_hist = []
    for strain in regions:
        for chrm in regions[strain]:
            entries = sorted(regions[strain][chrm], \
                                 key = lambda x: x['region_start'])
            for i in range(1, len(entries)):
                # why does this happen? 2 reasons: usually because one
                # entry on forward strand and one on reverse but in
                # just a couple cases it seems like there's the same
                # sequence in multiple alignment blocks 
                if entries[i]['region_start'] - entries[i-1]['region_end'] < 0:
                    print '********'
                    print strain, chrm
                    print entries[i]
                    print entries[i-1]
                else:
                    gap_hist.append(entries[i]['region_start'] - \
                                        entries[i-1]['region_end'])

    fn = gp.analysis_out_dir_absolute + 'summary_gaps_' + tag + '.txt'
    f_out = open(fn, 'w')

    for g in gap_hist:
        f_out.write(str(g) + ' ')
    f_out.write('\n')

    f_out.write(str(mystats.mean(gap_hist)) + ' ' + \
                    str(mystats.std_err(gap_hist)) + ' ')
    bs = mystats.bootstrap(gap_hist)
    f_out.write(str(bs[0]) + ' ' + str(bs[1]) + '\n')

    f_out.close()

# lengths of introgressed regions
def get_length_hist(regions):

    length_hist = []
    for strain in regions:
        for chrm in regions[strain]:
            for entry in regions[strain][chrm]:
                length_hist.append(entry['region_end'] - entry['region_start'] + 1)
            
    fn = gp.analysis_out_dir_absolute + 'summary_lengths_' + tag + '.txt'
    f_out = open(fn, 'w')

    for l in length_hist:
        f_out.write(str(l) + ' ')
    f_out.write('\n')

    f_out.write(str(mystats.mean(length_hist)) + ' ' + \
                    str(mystats.std_err(length_hist)) + ' ')
    bs = mystats.bootstrap(length_hist)
    f_out.write(str(bs[0]) + ' ' + str(bs[1]) + '\n')

    f_out.close()

# lengths of introgressed regions only including ungapped columns
def get_non_gap_hist(regions):

    length_hist = []
    for strain in regions:
        for chrm in regions[strain]:
            for entry in regions[strain][chrm]:
                length_hist.append(entry['number_non_gap'])
            
    fn = gp.analysis_out_dir_absolute + 'summary_non_gap_' + tag + '.txt'
    f_out = open(fn, 'w')

    for l in length_hist:
        f_out.write(str(l) + ' ')
    f_out.write('\n')

    f_out.write(str(mystats.mean(length_hist)) + ' ' + \
                    str(mystats.std_err(length_hist)) + ' ')
    bs = mystats.bootstrap(length_hist)
    f_out.write(str(bs[0]) + ' ' + str(bs[1]) + '\n')

    f_out.close()


# average length and total amount of introgressed sequence in different strains
def get_lengths_by_strain(regions):

    fn = gp.analysis_out_dir_absolute + 'summary_lengths_by_strain_' + tag + '.txt'
    f_out = open(fn, 'w')
    for strain in regions:
        lengths = []
        for chrm in regions[strain]:
            lengths += [x['region_end'] - x['region_start'] + 1 \
                            for x in regions[strain][chrm]]
        f_out.write(strain + ' ')
        f_out.write(str(mystats.mean(lengths)) + ' ' + \
                        str(mystats.std_err(lengths)) + ' ')
        bs = mystats.bootstrap(lengths)
        f_out.write(str(bs[0]) + ' ' + str(bs[1]) + ' ')
        f_out.write(str(sum(lengths)))
        f_out.write('\n')

    f_out.close()

# get number of strains that each gene is introgressed in (given that
# it's introgressed in at least one strain)
def get_strains_by_gene(regions):
    gene_dic = {}
    for strain in regions:
        for chrm in regions[strain]:
            # can do it at this level because same gene will never be
            # on multiple chromosomes
            strain_genes = set()
            for entry in regions[strain][chrm]:
                for gene in entry['genes']:
                    strain_genes.add(gene)
            for gene in strain_genes:
                if gene not in gene_dic:
                    gene_dic[gene] = []
                gene_dic[gene].append(strain)

    genes_sorted = sorted(gene_dic.keys(), \
                              key = lambda x: len(gene_dic[x]))
                    
    # write strains for each gene
    fn = gp.analysis_out_dir_absolute + 'summary_strains_by_gene_' + tag + '.txt'
    f_out = open(fn, 'w')
    for gene in genes_sorted:
        f_out.write(gene + ' ')
        f_out.write(' '.join(gene_dic[gene]))
        f_out.write('\n')
    f_out.close()

    # write only number of strains
    fn = gp.analysis_out_dir_absolute + 'summary_strains_by_gene_count_' + tag + '.txt'
    f_out = open(fn, 'w')
    for gene in genes_sorted:
        f_out.write(gene + ' ')
        f_out.write(str(len(gene_dic[gene])))
        f_out.write('\n')
    f_out.close()
