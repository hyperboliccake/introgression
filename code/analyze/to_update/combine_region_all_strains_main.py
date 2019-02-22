# input a gene or start/end coordinates
# output a multiple alignment file
# - for gene, relies on annotations/orfs
# - for coordinates, relies on alignments

import re
import sys
import os
import math
import Bio.SeqIO
import copy
from combine_all_strains import *
import gene_predictions
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../align/')
import align_helpers
sys.path.insert(0, '../misc/')
import read_table
import read_fasta
import write_fasta
import mystats


tag = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])
chrm = sys.argv[4]

print 'getting range sequence from reference strain'
ref = 'S288c'
ref_seq_fn = gp.ref_dir[ref] + gp.ref_fn_prefix[ref] + '_chr' + chrm + gp.fasta_suffix
ref_range_seq = get_range_seq(start, end, ref_seq_fn)

print 'getting range sequences from all strains'
s = align_helpers.get_strains(align_helpers.flatten(gp.non_ref_dirs.values()))
s.append((gp.ref_fn_prefix[gp.alignment_ref_order[1]], gp.ref_dir[gp.alignment_ref_order[1]]))
# keyed by strain: (seq, start, end)
strain_range_seqs = get_range_seqs(s, chrm, start, end, tag)

print 'writing all range sequences to file'
keys = sorted(strain_range_seqs.keys())
headers = [key + ' ' + str(strain_range_seqs[key][1]) + ':' + \
           str(strain_range_seqs[key][2]) for key in keys]
seqs = [strain_range_seqs[key][0] for key in keys]
strains = [ref] + keys
headers = [ref + ' ' + str(start) + ':' + str(end)] + headers
seqs = [ref_range_seq] + seqs
label = 'chr' + chrm + '_' + str(start) + '-' + str(end)
range_seqs_fn = gp.analysis_out_dir_absolute + tag + '/ranges/' + \
                 label + '/' + label + gp.fasta_suffix
if not os.path.isdir(gp.analysis_out_dir_absolute + tag + '/ranges/' + label):
    os.makedirs(gp.analysis_out_dir_absolute + tag + '/ranges/' + label)
write_fasta.write_fasta(headers, seqs, range_seqs_fn)

suffixes = ['', '_filtered']
for suffix in suffixes:
    print ' '.join(['finding', suffix, 'regions that overlap range'])
    # read in filtered regions
    fn_regions = gp.analysis_out_dir_absolute + tag + '/' + \
                 'introgressed_blocks' + suffix + '_par_' + tag + '_summary_plus.txt'
    regions, l = read_table.read_table_rows(fn_regions, '\t')

    regions_overlapping = {}
    # TODO does this actually ensure that regions are sorted appropriately
    # in fasta headers below?
    region_keys_ordered = sorted(regions.keys(), key=lambda x: int(x[1:]))
    for region in region_keys_ordered:
        if regions[region]['chromosome'] == chrm and \
           ((int(regions[region]['start']) > start and \
             int(regions[region]['start']) < end) or \
            (int(regions[region]['end']) > start and \
             int(regions[region]['end']) < end)):
            strain = regions[region]['strain']
            if not regions_overlapping.has_key(strain):
                regions_overlapping[strain] = []
            regions_overlapping[strain].append(region)

    print ' '.join(['writing all range sequences to file, with', \
                   suffix, 'introgressed bases capitalized'])
    headers_current = copy.deepcopy(headers)
    seqs_current = copy.deepcopy(seqs)
    for i in range(len(seqs)):
        strain = strains[i]
        seq = seqs_current[i].lower()
        header = headers_current[i]
        seqs_current[i] = seq
        headers_current[i] = header
        if strain not in regions_overlapping:
            continue
        r = strain_range_seqs[strain]
        t, labels = read_table.read_table_columns(gp.analysis_out_dir_absolute + \
                                                  tag + '/' + \
                                                  'site_summaries/predictions_' + \
                                                  strain + \
                                                  '_chr' + chrm + \
                                                  '_site_summary.txt.gz', '\t')
        ref_ind_to_strain_ind = dict(zip(t['ps_ref'], t['ps_strain']))
        for region in regions_overlapping[strain]:
            header += ' ' + region
            start_strain = math.ceil(float(\
                                ref_ind_to_strain_ind[regions[region]['start']]))
            end_strain = math.floor(float(\
                                ref_ind_to_strain_ind[regions[region]['end']]))
            start_relative = int(max(start_strain - int(r[1]), 0))
            end_relative = int(end_strain - int(r[1]))
            seq = seq[:start_relative] + \
                  seq[start_relative:end_relative+1].upper() + \
                  seq[end_relative+1:]
        seqs_current[i] = seq
        headers_current[i] = header

    range_seqs_fn = gp.analysis_out_dir_absolute + tag + '/ranges/' + label + \
                    '/' + label + '_introgressed' + suffix + gp.fasta_suffix
    write_fasta.write_fasta(headers_current, seqs_current, range_seqs_fn)


    print 'aligning range sequences'
    range_seqs_aligned_fn = range_seqs_fn.replace(gp.fasta_suffix, gp.alignment_suffix)
    cmd_string = gp.mafft_install_path + '/mafft ' + \
                 ' --reorder --preservecase ' + \
                 range_seqs_fn + ' > ' + range_seqs_aligned_fn
    os.system(cmd_string)
