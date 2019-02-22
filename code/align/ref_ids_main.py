import re
import sys
import os
import copy
from collections import defaultdict
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc/')
import mystats
import seq_functions
import read_fasta

# get pairwise identities between all aligned references:
# - overall average
# - by chromosome
# - in windows across genome

window = 100

gp_dir = '../'

nrefs = len(gp.alignment_ref_order)

pair_chrm_ids = defaultdict(lambda: defaultdict(list))
for chrm in gp.chrms:
    print chrm
    fn = gp_dir + gp.alignments_dir + \
         '_'.join(gp.alignment_ref_order) + \
         '_chr' + chrm + '_mafft' + gp.alignment_suffix

    headers, seqs = read_fasta.read_fasta(fn)

    for i in range(nrefs):
        ref1 = gp.alignment_ref_order[i]
        for j in range(i+1, nrefs):
            print i, j
            ref2 = gp.alignment_ref_order[j]

            ids = seq_functions.seq_id_windowed(seqs[i], seqs[j], window)
        
            pair_chrm_ids[(ref1, ref2)][chrm] = ids

fs = open(gp.analysis_out_dir_absolute + 'ref_ids_summary_' + \
         '_'.join(gp.alignment_ref_order) + '.txt', 'w')
fs.write('pair\tchromosome\tmean\tmedian\n')

f = open(gp.analysis_out_dir_absolute + 'ref_ids_' + \
         '_'.join(gp.alignment_ref_order) + '.txt', 'w')
f.write('pair\tid\n')

for pair in pair_chrm_ids.keys():
    all_ids = []
    pair_string = ','.join(pair)
    for chrm in gp.chrms:
        ids = pair_chrm_ids[pair][chrm]
        fs.write(pair_string + '\t' + \
                 chrm + '\t' + \
                 str(mystats.mean(ids)) + '\t' + \
                 str(mystats.median(ids)) + '\n')
        all_ids += ids
    fs.write(pair_string + '\t' + \
             'all' + '\t' + \
             str(mystats.mean(all_ids)) + '\t' + \
             str(mystats.median(all_ids)) + '\n')

    for i in ids:
        f.write(pair_string + '\t' + str(i) + '\n')

f.close()
