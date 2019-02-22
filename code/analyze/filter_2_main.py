# two levels of filtering:
# 1. remove regions that don't look confidently introgressed at all,
#    based on fraction gaps/masked, number of matches to S288c and not S288c
#    --> _filtered1 
# 2. remove regions that we can't confidently pin on a specific reference,
#    based on whether it matches similarly to other reference(s)
#    --> _filtered2

# do second level of filtering here, based on previously selected
# thresholds

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

args = predict.process_predict_args(sys.argv[2:])
threshold = float(sys.argv[1])

for species_from in args['known_states'][1:]:

    print species_from

    fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
         'blocks_' + species_from + \
         '_' + args['tag'] + '_filtered1.txt'
    region_summary, fields = read_table.read_table_rows(fn, '\t')

    fields2i = fields + ['alternative_states', 'alternative_ids', \
                         'alternative_P_counts']
    fields2 = fields

    fn_out2i = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
              'blocks_' + species_from + \
              '_' + args['tag'] + '_filtered2intermediate.txt'

    fn_out2 = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
              'blocks_' + species_from + \
              '_' + args['tag'] + '_filtered2.txt'

    f_out2i = open(fn_out2i, 'w')
    f_out2i.write('\t'.join(fields2i) + '\n')

    f_out2 = open(fn_out2, 'w')
    f_out2.write('\t'.join(fields2) + '\n')

    for region_id in region_summary:
        #print region_id, '****'
        region = region_summary[region_id]
        headers, seqs = read_fasta.read_fasta(gp.analysis_out_dir_absolute + \
                                              args['tag'] + \
                                              '/regions/' + region_id + '.fa.gz', \
                                              gz = True)
        info_string = seqs[-1]
        seqs = seqs[:-1]
 
        # filtering stage 2: things that we're confident in calling
        # introgressed from one species specifically
        p, alt_states, alt_ids, alt_P_counts = passes_filters2(region, seqs, threshold)
        region['alternative_states'] = ','.join(alt_states)
        region['alternative_ids'] = ','.join([str(x) for x in alt_ids])
        region['alternative_P_counts'] = ','.join([str(x) for x in alt_P_counts])
        write_filtered_line(f_out2i, region_id, region, fields2i)

        if p:
            write_filtered_line(f_out2, region_id, region, fields2)

    f_out2i.close()
    f_out2.close()
