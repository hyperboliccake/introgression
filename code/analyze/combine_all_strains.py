# input a gene or start/end coordinates
# output a multiple alignment file
# - for gene, relies on annotations
# - for coordinates, relies on alignments

import re
import sys
import os
import math
import Bio.SeqIO
import copy
import gene_predictions
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../sim/')
import sim_analyze_hmm_bw as sim
sys.path.insert(0, '../misc/')
import read_table
import mystats

def get_gene_seqs(fn, gene, chrm):
    # get gene sequence for each strain
    gb_records = Bio.SeqIO.parse(fn, 'genbank')
    strain_gene_seqs = {}
    strains = set([])
    for strain_chrm_record in gb_records:
        desc = strain_chrm_record.description
        m = re.search(' (?P<strain>[a-zA-Z0-9]+) chromosome (?P<chrm>[IVXM]+)', \
                      desc)
        chrm_current = m.group('chrm')
        strain = m.group('strain').lower()
        strains.add(strain)
        #if len(strain_gene_seqs) > 82:
        #    break
        print strain, chrm_current
        if chrm_current != chrm:
            continue
        for feature in strain_chrm_record.features:
            if feature.type == 'CDS' and feature.qualifiers.has_key('gene') and \
               feature.qualifiers['gene'][0] == gene:
                desc = strain_chrm_record.description
                m = re.search(\
                        ' (?P<strain>[a-zA-Z0-9]+) chromosome (?P<chrm>[IVXM]+)', \
                        desc)
                seq = str(feature.extract(strain_chrm_record.seq).lower())
                start = str(feature.location.start)
                end = str(feature.location.end)
                strand = str(feature.location.strand)
                locus_tag = feature.qualifiers['locus_tag'][0]
                strain_gene_seqs[strain] = {'seq':seq, \
                                            'chrm':chrm, \
                                            'start':start, \
                                            'end':end, \
                                            'strand':strand,\
                                            'locus_tag':locus_tag}
            
                print '- found gene in', strain
    return strain_gene_seqs, list(strains)

# because don't have gb file for paradoxus...
def get_gene_seqs_fsa(fn, gene, chrm):
    f = open(fn, 'r')
    line = f.readline()
    m = 'SCER:' + gene
    while line != '':
        if m in line:
            seq = ''
            line = f.readline()
            while not line.startswith('>'):
                seq += line[:-1]
                line = f.readline()
            f.close()

            seqfa = open(gp.ref_dir['CBS432'] + 'CBS432_chr' + chrm + '.fa', 'r').read()
            seqfa = seqfa.replace('\n', '')
            if seq in seqfa:
                print 'found paradoxus seq'
            else:
                print 'did not find paradoxus seq'
            fg = open('a.txt', 'w')
            fg.write(seq + '\n')
            fg.write(seqfa + '\n')
            fg.close()
            return seq.lower()

        line = f.readline()
        

tag = sys.argv[1]
gene = sys.argv[2]
chrm = sys.argv[3]

fx = open(gp.analysis_out_dir_absolute + tag + '/' + \
         'genes/' + gene + '.txt', 'w')

print 'getting gene sequences'
strain_gene_seqs, strains = get_gene_seqs(gp.gb_all, gene, chrm)
ref_seqs = {}
a = get_gene_seqs(gp.ref_gb_dir['S288c'] + 'S288c_chr' + chrm + '.gb', \
                                  gene, chrm)[0].values()[0]
ref_seqs['S288c'] = a['seq']
locus_tag = a['locus_tag']
ref_seqs['CBS432'] = get_gene_seqs_fsa(gp.ref_gb_dir['CBS432'] + 'CBS432.fsa', \
                                       locus_tag, chrm)

# read in filtered regions
fn_regions = gp.analysis_out_dir_absolute + tag + '/' + \
             'introgressed_blocks_filtered_par_' + tag + '_summary_plus.txt'
regions, l = read_table.read_table_rows(fn_regions, '\t')

# figure out which strains are introgressed/which regions overlap gene
print 'finding introgressed regions overlapping gene'
fn_genes_regions = gp.analysis_out_dir_absolute + tag + '/' + \
                   'genes_for_each_region_chr' + chrm + '_' + tag + '.txt'
region_to_genes = gene_predictions.read_genes_for_each_region_summary(fn_genes_regions)
regions_overlapping = dict(zip(strains, [[] for s in strains]))
print sorted(regions_overlapping.keys())
for region in regions:
    if regions[region]['chromosome'] == chrm and \
       gene in [x[0] for x in region_to_genes[region]['gene_list']]:
        strain = regions[region]['strain']
        print regions_overlapping
        regions_overlapping[strain].append(region)

for strain in strains:
    fx.write(strain + '\t' + '\t'.join(regions_overlapping[strain]))
    if strain not in strain_gene_seqs:
        fx.write('missing')
    fx.write('\n')
fx.close()
    

print '********', regions_overlapping
# put all gene sequences in one fasta file, with introgressed bases capitalized
f = open(gp.analysis_out_dir_absolute + tag + '/' + \
         'genes/' + gene + gp.fasta_suffix, 'w')
for ref in ref_seqs:
    f.write('> ' + ref + '\n' + ref_seqs[ref] + '\n')

for strain in strains:
    print strain
    fn_t = gp.analysis_out_dir_absolute + tag + '/' + \
           'site_summaries/' + \
           'predictions_' + strain + '_chr' + chrm +  '_site_summary.txt.gz'
    t = read_table.read_table_columns(fn_t, '\t')[0]
    ref_ind_to_strain_ind = dict(zip(t['ps_ref'], t['ps_strain']))
    g = strain_gene_seqs[strain]
    header = '> ' + strain + ' ' + g['chrm'] + ' ' + g['start'] + ' ' + \
             g['end'] + ' ' + g['strand']
    seq = g['seq']
    for region in regions_overlapping[strain]:
        header += ' ' + region
        start_strain = math.ceil(float(ref_ind_to_strain_ind[regions[region]['start']]))
        end_strain = math.floor(float(ref_ind_to_strain_ind[regions[region]['end']]))
        start_relative = max(start_strain - int(g['start']), 0)
        end_relative = end_strain - int(g['start'])
        seq = seq[:start_relative] + \
              seq[start_relative:end_relative+1].upper() + \
              seq[end_relative+1:]

    # write to file
    f.write(header + '\n'  + seq + '\n')
    
f.close()
