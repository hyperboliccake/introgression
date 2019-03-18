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
import numpy as np
import read_args
import summarize_region_quality
from filter_helpers import *
import global_params as gp
from misc import read_table
from misc import read_fasta

args = read_args.process_predict_args(sys.argv[2:])
threshold = float(sys.argv[1])

for species_from in args['known_states'][1:]:

    print(species_from)

    fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
         'blocks_' + species_from + \
         '_' + args['tag'] + '_filtered1.txt'
    region_summary, fields = read_table.read_table_rows(fn, '\t')

    fields2i = fields + ['predicted_species_original', 'alternative_ids', \
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

    regions_fn = gp.analysis_out_dir_absolute + args['tag'] + '/regions/' + \
                 species_from + gp.fasta_suffix + '.gz'
    region_seqs = summarize_region_quality.read_region_file(regions_fn)

    for region_id in region_summary:

        region = region_summary[region_id]

        info_string = region_seqs[region_id]['info']['seq']
        seqs = np.asarray([list(region_seqs[region_id][ref]['seq']) \
                           for ref in args['known_states']])

        # filtering stage 2: things that we're confident in calling
        # introgressed from one species specifically
        p, alt_states, alt_ids, alt_P_counts = passes_filters2(region, seqs, \
                                                               threshold, \
                                                               args['known_states'])
        region['alternative_states'] = '/'.join(alt_states)
        region['alternative_ids'] = '/'.join([str(x) for x in alt_ids])
        region['alternative_P_counts'] = '/'.join([str(x) for x in alt_P_counts])

        region['predicted_species_original'] = region['predicted_species']
        region['predicted_species'] = region['alternative_states']
        write_filtered_line(f_out2i, region_id, region, fields2i)

        if p:
            write_filtered_line(f_out2, region_id, region, fields2)

    f_out2i.close()
    f_out2.close()
