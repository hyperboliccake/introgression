import sys
import os
import gzip
import predict
from summary_plus import *
import gene_predictions
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc/')
import read_fasta
import read_table
import seq_functions

args = predict.process_predict_args(sys.argv[1:])
gp_dir = '../'

regions_all = {}
strains = set([])
for species_from in args['states']:

    fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
         'introgressed_blocks_' + species_from + \
         '_' + args['tag'] + '_labeled.txt'

    # region_id strain chromosome predicted_species start end number_non_gap
    d, labels = read_table.read_table_columns(fn, '\t')

    for s in args['known_states']:
        d['id_' + s] = [0 for i in range(len(d['region_id']))]

    regions_all[species_from] = d

    strains = strains.union(set(d['strain']))

# loop through chromosomes and strains, followed by species of
# introgression so that we only have to read each alignment in once
for chrm in gp.chrms[:1]:

    for strain in strains:

        fn = gp_dir + gp.alignments_dir + \
             '_'.join(gp.alignment_ref_order) + '_' + strain + \
             '_chr' + chrm + '_mafft' + gp.alignment_suffix
        headers, seqs = read_fasta.read_fasta(fn)

        ind_align = index_alignment_by_reference(seqs[0])

        for si in range(len(args['states'])):
            
            species_from = args['states'][si]
            regions = regions_all[species_from]

            for i in range(len(regions['region_id'])):

                if regions['chromosome'][i] != chrm or regions['strain'][i] != strain:
                    continue

                fn_region = gp.analysis_out_dir_absolute + args['tag'] + '/' \
                            'regions/' + regions['region_id'][i] + \
                            gp.fasta_suffix + '.gz'
                f_region = gzip.open(fn_region, 'wb')
                    
                # calculate:
                # - identity with each reference
                # - fraction of region that is gapped/masked

                slice_start = ind_align[int(regions['start'][i])]
                slice_end = ind_align[int(regions['end'][i])] + 1

                seqx = seqs[-1][slice_start:slice_end]

                for sj in range(len(args['known_states'])):
                    seqj = seqs[sj][slice_start:slice_end]
                    sid = seq_functions.seq_id(seqj, seqx)
                    regions_all[species_from]['id_' + args['known_states'][sj]][i] = sid
                    f_region.write('> ' + args['known_states'][sj] + '\n')
                    f_region.write(seqj + '\n')

                f_region.write('> ' + strain + '\n')
                f_region.write(seqx + '\n')

labels = labels + ['id_' + x for x in args['known_states']]

for species_from in args['states']:

    fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
         'introgressed_blocks_' + species_from + \
         '_' + args['tag'] + '_quality.txt'

    f = open(fn, 'w')
    f.write('\t'.join(labels) + '\n')

    for i in range(len(regions_all[species_from]['region_id'])):
        if regions_all[species_from]['chromosome'][i] != 'I':
            continue
        f.write('\t'.join([str(regions_all[species_from][label][i]) for label in labels]))
        f.write('\n')
    f.close()

"""
# copy pasta :(
fields = ['strain', 'chromosome', 'predicted_species', 'start', 'end', \
          'number_non_gap', 'number_match_ref1', 'number_match_ref2', \
          'number_match_only_ref1', 'number_match_ref2_not_ref1', \
          'number_mismatch_all_ref']
regions = {}
for chrm in gp.chrms:
    fn = gp.analysis_out_dir_absolute + tag + '/' + \
         'introgressed_blocks_chr' + chrm + \
         '_par_' + tag + '_quality.txt'
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
"""
