import re
import sys
import os
import copy
from process import *
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../sim/')
import sim_analyze_hmm_bw as sim
sys.path.insert(0, '../misc/')
import read_maf

tag, topology, species_to, species_from1, species_from2, \
    num_samples_species_to, num_samples_species_from1, num_samples_species_from2, \
    N0_species_to, N0_species_from1, N0_species_from2, \
    migration_from1, migration_from2, \
    expected_tract_lengths, \
    expected_num_tracts, \
    has_ref_from1, has_ref_from2, \
    rho, outcross_rate, theta, num_sites, num_reps = \
    sim.process_args(sys.argv)


#####
# write gene summary file (gene chromosome start end)
#####

fn = gp.analysis_out_dir_absolute + tag +'introgressed_hmm_' + tag + \
    '_genes_summary_filtered.txt'
f = open(fn, 'r')

#####
# write strain list
#####

#####
# write file with regions for each introgressed gene
#####

gp_dir = '../'
fn_all_regions = gp.analysis_out_dir_absolute + 'introgressed_hmm_' + tag + '.txt'
# introgressed regions keyed by strain and then chromosome
regions = read_regions(fn_all_regions)
