# 
# columns::
# gene name
# average fraction introgressed
# number of strains introgressed
# average id with cer ref
# list of introgressed regions
# list of region lengths, correponding to regions
# list of strains
# list of introgressed fractions, corresponding to strains
# list of ids with cer ref (across whole gene) [from gene alignments]

import re
import sys
import os
import copy
from datetime import datetime, timedelta
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc/')
import read_fasta
import read_table
import overlap
import mystats
import seq_functions

two_days_ago = datetime.now() - timedelta(days=2)

##======
# read in analysis parameters
##======

tag = sys.argv[1]
suffix = ''
if len(sys.argv) == 3:
    suffix = sys.argv[2]

##======
# read in introgressed regions 
##======

gp_dir = '../'

blocks_fn = gp.analysis_out_dir_absolute + tag + '/' + \
            'introgressed_blocks' + suffix + '_par_' + tag + '_summary_plus.txt'
regions, l = read_table.read_table_rows(blocks_fn, '\t')

##======
# read in gene list
##======

genes_f = open('../../data/S288c_verified_orfs.tsv', 'r')
labels = ['systematic name', 'name', 'description', \
          'strand', 'start', 'end', 'chromosome']
gene_list = [dict(zip(labels, line[:-1].split('\t'))) for line in genes_f.readlines()]
for g in gene_list:
    g['start'] = int(g['start'])-1
    g['end'] = int(g['end'])-1
    g['chromosome'] = g['chromosome'][3:]
    if g['name'] == '""':
        g['name'] = g['systematic name']
genes_f.close()

paralogs, l = read_table.read_table_rows('../../data/S288c_paralogs.tsv', '\t', \
                                         header=False)

gene_summary = {}
for region_id in regions:
    region = regions[region_id]
    region_start = int(region['start'])
    region_end = int(region['end'])
    region_chrm = region['chromosome']
    strain = region['strain']
    for g in gene_list:
        if g['chromosome'] == region_chrm:
            o = overlap.overlap(region_start, region_end, g['start'], g['end'])
            if o > 0:
                if not gene_summary.has_key(g['name']):
                    gene_summary[g['name']] = {'regions':[], 'strains':[], \
                                               'fractions':[], 'lengths':[], \
                                               'strain_fractions':{}, \
                                               'chromosome':g['chromosome'], \
                                               'start':g['start'], 'end':g['end'], \
                                               'paralog':''}
                    if paralogs.has_key(g['name']):
                        gene_summary[g['name']]['paralog'] = paralogs[g['name']][3]

    
                gene_length = gene_summary[g['name']]['end'] - \
                              gene_summary[g['name']]['start'] + 1
                gene_summary[g['name']]['regions'].append(region_id)
                gene_summary[g['name']]['strains'].append(strain)
                gene_summary[g['name']]['lengths'].append(int(region['end']) - \
                                                          int(region['start']) + 1)
                if not gene_summary[g['name']]['strain_fractions'].has_key(strain):
                    gene_summary[g['name']]['strain_fractions'][strain] = 0
                gene_summary[g['name']]['strain_fractions'][strain] += \
                    o / float(gene_length)

f = open(gp.analysis_out_dir_absolute + tag + '/gene_summary' + \
         suffix + '_' + tag + '.txt', 'w')
f.write('\t'.join(['gene', \
                   'chromosome', \
                   'start', \
                   'end', \
                   'paralog', \
                   'num_strains', \
                   'avg_frac_intd', \
                   'average_cer_ref_id',\
                   'avg_region_length', \
                   'strains', \
                   'fractions', \
                   'cer_ref_ids', \
                   'regions', \
                   'lengths']) + '\n')

for gene in sorted(gene_summary.keys()):
    gene_summary[gene]['average fraction'] = \
        mystats.mean(gene_summary[gene]['strain_fractions'].values())
    gene_summary[gene]['average length'] = \
        mystats.mean(gene_summary[gene]['lengths'])

    strains = gene_summary[gene]['strain_fractions'].keys()
    strain_fractions = \
        [str(gene_summary[gene]['strain_fractions'][strain]) for strain in strains]
    
    r = gene_summary[gene]['regions']
    r_lengths = [str(x) for x in gene_summary[gene]['lengths']]

    fn_gene = gp.analysis_out_dir_absolute + tag + '/genes/' + \
              gene + '/' + gene + '_introgressed' + gp.alignment_suffix
    if not os.path.isfile(fn_gene):# or datetime.fromtimestamp(os.path.getmtime(fn_gene)) < two_days_ago:
        print 'skipping', gene
        continue
        #print 'getting gene alignment for', gene
        #sys.stdout.flush()
        #os.system('python combine_gene_all_strains_main.py ' + tag + ' ' +  \
        #          gene + ' ' + gene_summary[gene]['chromosome'])
    #print 'working on', gene
    headers, seqs = read_fasta.read_fasta(fn_gene)
    aligned_seqs = {}
    for i in range(len(headers)):
        s = headers[i].split()[0][1:]
        aligned_seqs[s] = seqs[i].lower()

    strain_ids = []
    try:
        for strain in strains:
            strain_ids.append(seq_functions.seq_id(aligned_seqs['S288c'], aligned_seqs[strain]))
    except:
        print gene, 'not finished'
        continue
    avg_strain_id = mystats.mean(strain_ids)
    strain_ids = [str(x) for x in strain_ids]
        
    
    f.write('\t'.join([gene, \
                       gene_summary[gene]['chromosome'], \
                       str(gene_summary[gene]['start']), \
                       str(gene_summary[gene]['end']), \
                       gene_summary[gene]['paralog'],\
                       str(len(strains)), \
                       str(gene_summary[gene]['average fraction']), \
                       str(avg_strain_id), \
                       str(gene_summary[gene]['average length']), \
                       ','.join(strains), \
                       ','.join(strain_fractions), \
                       ','.join(strain_ids), \
                       ','.join(r), \
                       ','.join(r_lengths)]) + '\n')



                #headers, seqs = read_fasta.read_fasta(gp.analysis_out_dir_absolute + \
                #                                      tag + '/regions/' + \
                #                                      region_id + \
                #                                      gp.alignment_suffix + '.gz')
                



f.close()
