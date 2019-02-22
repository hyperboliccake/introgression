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


tag = sys.argv[1]
suffix = ''
if len(sys.argv == 3):
    suffix = sys.argv[2]

fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
     'introgressed_blocks_par' + suffix + '_' + args['tag'] + '_summary_plus.txt'
region_summary = gene_predictions.read_region_summary(fn)

sep = '\t'

##======
# for plot: lengths of all introgressed regions
##======

# one table for each tag
# strain chrm region_length

lengths_all = []
for region in region_summary:
    length = int(region_summary[region]['end']) - \
             int(region_summary[region]['start']) + 1
