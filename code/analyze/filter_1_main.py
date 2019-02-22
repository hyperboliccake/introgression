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
import predict
from filter_helpers import *
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc/')
import read_table
import read_fasta

args = predict.process_predict_args(sys.argv[1:])

for species_from in args['known_states'][1:]:

    print species_from

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

    for region_id in region_summary:
        #print region_id, '****'
        region = region_summary[region_id]
        headers, seqs = read_fasta.read_fasta(gp.analysis_out_dir_absolute + \
                                              args['tag'] + \
                                              '/regions/' + region_id + '.fa.gz', \
                                              gz = True)
        info_string = seqs[-1]
        seqs = seqs[:-1]
 
        # filtering stage 1: things that we're confident in calling not
        # S288c
        p, reason = passes_filters1(region, info_string)
        region['reason'] = reason
        write_filtered_line(f_out1i, region_id, region, fields1i)

        if p:
            write_filtered_line(f_out1, region_id, region, fields1)

    f_out1i.close()
    f_out1.close()
