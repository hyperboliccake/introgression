# why the fuck are chromsome XIV entries duplicated??

import math
import numpy.random
import sys
import process_helpers
from summary_stats_helpers import *
sys.path.insert(0, '../misc/')
import mystats
sys.path.insert(0, '../')
import global_params as gp



tag = sys.argv[1]

fn = gp.analysis_out_dir_absolute + 'introgressed_hmm_' + tag + '_genes.txt'
regions = process_helpers.read_regions(fn)
regions = remove_duplicates(regions)

print 'gap hist'
get_gap_hist(regions)
print 'length hist'
get_length_hist(regions)
print 'non gap hist'
get_non_gap_hist(regions)
print 'lengths by strain'
get_lengths_by_strain(regions)
print 'strains by gene'
get_strains_by_gene(regions)
