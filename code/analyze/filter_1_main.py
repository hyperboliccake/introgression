# two levels of filtering:
# 1. remove regions that don't look confidently introgressed at all,
#    based on fraction gaps/masked, number of matches to S288c and not S288c
#    --> _filtered1 
# 2. remove regions that we can't confidently pin on a specific reference,
#    based on whether it matches similarly to other reference(s)
#    --> _filtered2

# just do the first level here, then run filter_2_thresholds_main.py
# to choose filtering thresholds for next level


import re
import sys
import os
import copy
import read_args
from filter_helpers import *
import summarize_region_quality
import global_params as gp
from misc import read_table
from misc import read_fasta

args = read_args.process_predict_args(sys.argv[1:])

for species_from in args['known_states'][1:]:

    print(species_from)

    fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
         'blocks_' + species_from + \
         '_' + args['tag'] + '_quality.txt'
    region_summary, fields = read_table.read_table_rows(fn, '\t')

    fields1i = fields + ['reason']
    fields1 = fields 

    fn_out1i = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
              'blocks_' + species_from + \
              '_' + args['tag'] + '_filtered1intermediate.txt'

    fn_out1 = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
              'blocks_' + species_from + \
              '_' + args['tag'] + '_filtered1.txt'

    f_out1i = open(fn_out1i, 'w')
    f_out1i.write('\t'.join(fields1i) + '\n')

    f_out1 = open(fn_out1, 'w')
    f_out1.write('\t'.join(fields1) + '\n')

    regions_fn = gp.analysis_out_dir_absolute + args['tag'] + '/regions/' + \
                 species_from + gp.fasta_suffix + '.gz'
    region_seqs = summarize_region_quality.read_region_file(regions_fn)

    for region_id in region_summary:

        region = region_summary[region_id]

        info_string = region_seqs[region_id]['info']['seq']
 
        # filtering stage 1: things that we're confident in calling not
        # S288c
        p, reason = passes_filters1(region, info_string, args['known_states'][0])
        region['reason'] = reason
        write_filtered_line(f_out1i, region_id, region, fields1i)

        if p:
            write_filtered_line(f_out1, region_id, region, fields1)

    f_out1i.close()
    f_out1.close()
