import re
import sys
import os
import copy
import gene_predictions
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../sim/')
import sim_analyze_hmm_bw as sim
sys.path.insert(0, '../misc/')
import mystats

# TODO make this more general table reading method
def read_region_summary_plus(fn):
    f = open(fn, 'r')
    labels = f.readline()[:-1].split('\t')
    regions = [line[:-1].split('\t') for line in f.readlines()]
    d = {}
    for region in regions:
        d[region[0]] = dict(zip(labels[1:], region[1:]))
    return labels[1:], d

def write_region_summary_plus_line(f, region_id, region, fields):
    f.write(region_id + '\t' + '\t'.join([region[field] for field in fields]))
    f.write('\n')

def passes_filters(region):
    
    # fraction gaps + masked filter
    fraction_gaps_masked_threshold = .5
    fraction_gaps_masked = \
        (float(region['number_gaps']) + float(region['number_masked_non_gap'])) / \
        (int(region['end']) - int(region['start']) + 1)
    if fraction_gaps_masked > fraction_gaps_masked_threshold:
        return False

    # number sites match only par filter
    number_match_only_threshold = 7
    number_match_only = int(region['number_match_ref2_not_ref1'])
    if number_match_only < number_match_only_threshold:
        return False

    # divergence from cer filter (idea is that poor alignments will
    # result in much larger divergence than we'd expect)
    id_ref1_threshold = .7
    id_ref1 = float(region['number_match_ref1']) / \
              (float(region['aligned_length']) - float(region['number_gaps']))
    if id_ref1 < id_ref1_threshold:
        return False
    
    return True

tag = sys.argv[1]

fn = gp.analysis_out_dir_absolute + tag + '/' + \
     'introgressed_blocks_par_' + tag + '_summary_plus.txt'
fields, region_summary = read_region_summary_plus(fn)

fn_out = gp.analysis_out_dir_absolute + tag + '/' + \
         'introgressed_blocks_filtered_par_' + tag + '_summary_plus.txt'
f_out = open(fn_out, 'w')
f_out.write('region_id\t' + '\t'.join(fields) + '\n')

for region_id in region_summary:
    region = region_summary[region_id]
    if passes_filters(region):
        write_region_summary_plus_line(f_out, region_id, region, fields)

f_out.close()
