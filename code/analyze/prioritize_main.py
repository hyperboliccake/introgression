# put introgressed genes in a file with relevant info to prioritize
# which ones look interesting

import sys
import os
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc/')
import read_table
import mystats

tag = sys.argv[1]
suffix = '_filtered'

# read in filtered regions
fn_regions = gp.analysis_out_dir_absolute + tag + '/' + \
             'introgressed_blocks' + suffix + '_par_' + tag + '_summary_plus.txt'
regions, l = read_table.read_table_rows(fn_regions, '\t')

# read in genes for each region
fn_genes_for_each_region = gp.analysis_out_dir_absolute + tag + '/' + \
                           'genes_for_each_region_' + tag + '.txt'
genes_for_each_region = {}
for line in open(fn_genes_for_each_region, 'r').readlines():
    line = line.split('\t')
    region = line[0]
    genes = line[2::2]
    fracs = [float(x) for x in line[3::2]]
    genes_for_each_region[line[0]] = dict(zip(genes, fracs))

# create dictionary keyed by gene to keep track of strains its
# introgressed in, and corresponding fraction introgressed, and
# regions
genes = {} #  gene:{'strain_fracs':{strain:frac, strain:frac}, 'regions':[region,...]}
for region in regions:
    for gene, frac in genes_for_each_region[region].iteritems():
        if not genes.has_key(gene):
            genes[gene] = {'strain_fracs':{}, 'regions':[]}
        strain = regions[region]['strain']
        if not genes[gene]['strain_fracs'].has_key(strain):
            genes[gene]['strain_fracs'][strain] = 0
        genes[gene]['strain_fracs'][strain] += frac
        genes[gene]['regions'].append(region)

# order genes by sum of fractions introgressed
genes_ordered = genes.keys()
genes_ordered.sort(key=lambda x: sum(genes[x]['strain_fracs'].values()), reverse=True)

# read in paralogs
fn_paralogs = '../../data/S288c_paralogs.tsv'
f_paralogs = open(fn_paralogs, 'r')
lines = [line.strip().split('\t') for line in f_paralogs.readlines()]
f_paralogs.close()
paralogs = {}
for line in lines:
    if line[0] != "":
        paralogs[line[0]] = line[3]


# write genes to file
fn = gp.analysis_out_dir_absolute + tag + \
     '/genes_prioritized.txt'
f = open(fn, 'w')
for gene in genes_ordered:

    # get reference identities
    #headers, seqs = read_fasta.read_fasta(fn_gene, 

    f.write(gene + '\t')
    if paralogs.has_key(gene):
        f.write(paralogs[gene])
    f.write('\t')
    strain_fracs = list(genes[gene]['strain_fracs'].iteritems())
    f.write(str(len(strain_fracs)) + '\t')
    f.write(str(mystats.mean([x[1] for x in strain_fracs])) + '\t')
    f.write(','.join([x[0] for x in strain_fracs]) + '\t')
    f.write(','.join([str(x[1]) for x in strain_fracs]))# + '\t')
    #f.write(','.join(genes[gene]['regions']) + '\t')
    f.write('\n')
f.close()


