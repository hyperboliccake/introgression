# explore different thresholds for calling introgressions for specific
# strains

# specifically, try a range of thresholds, and for each one, calculate
# fraction of introgressions we've classified as 1 strain or every
# possible combination of strains

# then we'll make some plots in R to see if there's a sort of obvious
# place to draw the line

import re
import sys
import os
import copy
from collections import defaultdict
import predict
from filter_helpers import *
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc/')
import read_table
import read_fasta

args = predict.process_predict_args(sys.argv[1:])

#thresholds = [.99, .98, .97, .96, .95, .94, .93, .92, .91, .9, .88, .85, .82, .8, .75, .7, .6, .5]
#thresholds = [.999, .995, .985, .975, .965, .955, .945, .935, .925, .915, .905, .89, .87, .86]
thresholds = [1]

open_mode = 'a'
f = open(gp.analysis_out_dir_absolute + args['tag'] + \
         '/filter_2_thresholds_' + args['tag'] + '.txt', open_mode)
if open_mode == 'w':
    f.write('threshold\tpredicted_state\talternative_states\tcount\n')
for threshold in thresholds:
    print threshold
    for species_from in args['known_states'][1:]:

        print '*', species_from

        fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
             'blocks_' + species_from + \
             '_' + args['tag'] + '_filtered1.txt'
        region_summary, fields = read_table.read_table_rows(fn, '\t')
        
        d = defaultdict(int)
        for region_id in region_summary:
            #print region_id, '****'
            region = region_summary[region_id]
            headers, seqs = read_fasta.read_fasta(gp.analysis_out_dir_absolute + \
                                                  args['tag'] + \
                                              '/regions/' + region_id + '.fa.gz', \
                                              gz = True)
            info_string = seqs[-1]
            seqs = seqs[:-1]
            
            p, alt_states, alt_ids, alt_P_counts = \
                passes_filters2(region, seqs, threshold)

            d[','.join(sorted(alt_states))] += 1
    
        for key in d:
            f.write(str(threshold) + '\t' + species_from + '\t' + \
                    key + '\t' + str(d[key]) + '\n')
f.close()
