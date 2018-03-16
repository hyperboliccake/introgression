# ps_cer ps_strain cer_ref par_ref strain gene introgressed_region cer_masked par_masked strain_masked 

import re
import sys
import os
import copy
import gene_predictions 
from annotate_regions import *
import predict
import pickle
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc/')
import read_fasta

##======
# read in analysis parameters
##======

refs, strains, args = predict.process_args(sys.argv[1:])
chrm = sys.argv[1]

##======
# read in introgressed/unknown regions and alignments
##======

gp_dir = '../'

blocks_fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
            'introgressed_blocks_par_' + args['tag'] + '.txt'
regions = predict.read_blocks(blocks_fn, region_id='r')


blocks_fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
            'introgressed_blocks_unknown_' + args['tag'] + '.txt'
regions_unk = predict.read_blocks(blocks_fn, region_id='u')


# for alignment files (input)
fn_align_prefix = gp_dir + gp.alignments_dir
fn_align_prefix += '_'.join([refs[s][0] for s in args['species']]) + '_'

##======
# produce annotated files
##======

# for keeping track of all genes introgressed in each strain, and the
# fraction introgressed
# keyed by strain, then gene, total fraction introgressed
# TODO maybe don't just pick par (and exclude unknown)
strain_genes_dic = dict(zip(regions.keys(), [{} for strain in regions.keys()]))
# the inverse of above
gene_strains_dic = {}

# just read genes from master reference (species[0]),
# since that's how the introgressed regions are indexed
master_ref = refs[args['species'][0]][0]

ref_labels = [refs[s][0] for s in args['species']]

# deal with just the one chromosome

# genbank file to read from
fn = gp.ref_gb_dir[master_ref] + master_ref + '_chr' + chrm + '.gb'

# for storing genes once they're read (or reading genes from if
# already exists)
fn_genes = gp.analysis_out_dir_absolute + '/' + \
           master_ref + '_chr' + chrm + '_genes.txt'

print 'reading genes on chromosome', chrm
# dictionary keyed by name: (start, end)
genes = gene_predictions.read_genes(fn, fn_genes)
print 'done reading genes'

# loop through all strains that we've called introgression in, and
# associate genes with the regions they overlap
for strain in regions.keys():
        
    print '***', strain, chrm
    sys.stdout.flush()

    fn_out = gp.analysis_out_dir_absolute + args['tag'] + '/site_summaries/' + \
         'predictions_' + strain + '_chr' + chrm +  '_site_summary.txt.gz'
    if not os.path.exists(os.path.dirname(fn_out)):
        os.makedirs(os.path.dirname(fn_out))

    # skip this strain x chromosome if there are no introgressed
    # regions for it
    if strain not in regions or chrm not in regions[strain]:
        continue

    # read alignment blocks for this strain and chromosome
    fn_align = fn_align_prefix + \
               strain + '_chr' + chrm + '_mafft' +  gp.alignment_suffix
    alignment_headers, alignment_seqs = read_fasta.read_fasta(fn_align)

    # read masked (unaligned) sequences
    seq_masked_fns = [header.split()[-1] for header in alignment_headers]
    seq_masked_fns = [mfn[:-len(gp.fasta_suffix)] + '_masked' + gp.fasta_suffix \
                      for mfn in seq_masked_fns]
    seqs_masked = [read_fasta.read_fasta(mfn)[1][0] for mfn in seq_masked_fns]

    labels = ref_labels + [strain]
    
    # mark each site as matching each reference or not
    ref_match_by_site = gene_predictions.get_ref_match_by_site(alignment_seqs, labels)
    # mark each site as in a gene or not
    genes_by_site = gene_predictions.get_genes_by_site(genes, alignment_seqs[0])
    # mark each site as introgressed or not
    all_regions = [regions[strain][chrm]]
    if regions_unk.has_key(strain) and regions_unk[strain].has_key(chrm):
        all_regions.append(regions_unk[strain][chrm])
    block_by_site = get_block_by_site(all_regions, alignment_seqs[0])

    write_predictions_annotated(alignment_headers, alignment_seqs, 0, \
                                ref_labels + [strain], ref_match_by_site, \
                                genes_by_site, block_by_site, seqs_masked, fn_out)

    
