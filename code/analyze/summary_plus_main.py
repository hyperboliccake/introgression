import sys
import os
import gzip
from summary_plus import *
import gene_predictions
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc/')
import read_fasta

# this is for adding a few columns to introgressed_blocks_par_summary file:
# - number of genes it overlaps
# - longest stretch of gaps


tag = sys.argv[1]

# copy pasta :(
fields = ['strain', 'chromosome', 'predicted_species', 'start', 'end', \
          'number_non_gap', 'number_match_ref1', 'number_match_ref2', \
          'number_match_only_ref1', 'number_match_ref2_not_ref1', \
          'number_mismatch_all_ref']
regions = {}
for chrm in gp.chrms:
    fn = gp.analysis_out_dir_absolute + tag + '/' + \
         'introgressed_blocks_chr' + chrm + \
         '_par_' + tag + '_summary.txt'
    d = gene_predictions.read_region_summary(fn)
    regions.update(d)

region_genes = {}
for chrm in gp.chrms:
    fn_genes = gp.analysis_out_dir_absolute + tag + '/' + \
               'genes_for_each_region_chr' + chrm + '_' + tag + '.txt'
    d = gene_predictions.read_genes_for_each_region_summary(fn_genes)
    region_genes.update(d)

fn_out = gp.analysis_out_dir_absolute + tag + '/' + \
        'introgressed_blocks_par_' + tag + '_summary_plus.txt'
fields.append('aligned_length')
fields.append('number_genes')
fields.append('number_gaps')
fields.append('number_masked')
fields.append('number_masked_non_gap')
fields.append('longest_gap')
fields.append('longest_mask')
fields.append('distance_from_telomere')
fields.append('distance_from_centromere')


i = 0
for region_id in regions:
    if i % 100 == 0:
        sys.stdout.write(str(i) + '/' + str(len(regions)) + '\r')
        sys.stdout.flush()
    i += 1

    regions[region_id]['number_genes'] = region_genes[region_id]['num_genes']

    fn_align = gp.analysis_out_dir_absolute + tag + '/' + \
               'regions/'  + region_id + '_masked.maf.gz'
    if not os.path.isfile(fn_align):
        print 'no masking for', region_id
        fn_align = gp.analysis_out_dir_absolute + tag + '/' + \
               'regions/'  + region_id + '.maf.gz'

    headers, seqs = read_fasta.read_fasta(fn_align, gz=True)

    regions[region_id]['aligned_length'] = len(seqs[0])

    number_gaps = gap_columns(seqs)
    regions[region_id]['number_gaps'] = number_gaps

    number_masked, number_masked_non_gap = masked_columns(seqs)
    regions[region_id]['number_masked'] = number_masked
    regions[region_id]['number_masked_non_gap'] = number_masked_non_gap

    longest_gaps = [longest_consecutive(seq, gp.gap_symbol) for seq in seqs]
    regions[region_id]['longest_gap'] = max(longest_gaps)

    longest_masks = [longest_consecutive(seq, gp.masked_symbol) for seq in seqs]
    regions[region_id]['longest_mask'] = max(longest_masks)

    regions[region_id]['distance_from_telomere'] = \
        distance_from_telomere(int(regions[region_id]['start']), \
                               int(regions[region_id]['end']), \
                               regions[region_id]['chromosome'])

    regions[region_id]['distance_from_centromere'] = \
        distance_from_centromere(int(regions[region_id]['start']), \
                                 int(regions[region_id]['end']), \
                                 regions[region_id]['chromosome'])
    
write_region_summary_plus(fn_out, regions, fields)
