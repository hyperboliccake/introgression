# this is for adding a few columns to introgressed_blocks_par_summary file:
# - number of genes it overlaps
# - longest stretch of gaps

import sys
import os
import gzip
import gene_predictions
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc/')
import read_fasta

def write_region_summary_plus(fn, regions, fields):
    f = open(fn, 'w')
    f.write('region_id\t' + '\t'.join(fields) + '\n')
    keys = sorted(regions.keys(), key=lambda x: int(x[1:]))
    for region_id in keys:
        f.write(region_id + '\t')
        f.write('\t'.join([str(regions[region_id][field]) for field in fields]))
        f.write('\n')
    f.close()

def gap_columns(seqs):
    g = 0
    for i in range(len(seqs[0])):
        for seq in seqs:
            if seq[i] == gp.gap_symbol:
                g += 1
                break
    return g

def longest_consecutive(s, c):
    max_consecutive = 0
    current_consecutive = 0
    in_segment = False
    for i in range(len(s)):
        if s[i] == c:
            current_consecutive += 1
            in_segment = True
        else:
            if in_segment:
                max_consecutive = max(max_consecutive, current_consecutive)
                current_consecutive = 0
            in_segment = False
    return max_consecutive

tag = sys.argv[1]
suffix = ''
if len(sys.argv) == 3:
    suffix = sys.argv[2]

fn = gp.analysis_out_dir_absolute + tag + '/' + \
     'introgressed_blocks' + suffix + '_par_' + tag + '_summary.txt'
# copy pasta :(
fields = ['strain', 'chromosome', 'predicted_species', 'start', 'end', \
          'number_non_gap', 'number_match_ref1', 'number_match_ref2', \
          'number_match_only_ref1', 'number_match_ref2_not_ref1', \
          'number_mismatch_all_ref']
regions = gene_predictions.read_region_summary(fn)

region_genes = {}
for chrm in gp.chrms:
    fn_genes = gp.analysis_out_dir_absolute + tag + '/' + \
               'genes_for_each_region_chr' + chrm + '_' + tag + '.txt'
    d = gene_predictions.read_genes_for_each_region_summary(fn_genes)
    region_genes.update(d)

fn_out = gp.analysis_out_dir_absolute + tag + '/' + \
        'introgressed_blocks' + suffix + '_par_' + tag + '_summary_plus.txt'
fields.append('aligned_length')
fields.append('number_genes')
fields.append('number_gaps')
fields.append('longest_gap')

i = 0
for region_id in regions:
    if i % 100 == 0:
        sys.stdout.write(str(i) + '/' + str(len(regions)) + '\r')
        sys.stdout.flush()
    i += 1

    regions[region_id]['number_genes'] = region_genes[region_id]['num_genes']

    fn_align = gp.analysis_out_dir_absolute + tag + '/' + \
               'regions/'  + region_id + '.maf.gz'
    headers, seqs = read_fasta.read_fasta(fn_align, gz=True)

    regions[region_id]['aligned_length'] = len(seqs[0])

    number_gaps = gap_columns(seqs)
    regions[region_id]['number_gaps'] = number_gaps

    longest_gaps = [longest_consecutive(seq, gp.gap_symbol) for seq in seqs]
    regions[region_id]['longest_gap'] = max(longest_gaps)
    
write_region_summary_plus(fn_out, regions, fields)
    
