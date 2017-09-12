import sys
import os
import process_args
import sim_process
from sim_actual import *

##======
# read in simulation parameters
##======
 
args = process_args.process_args()


##======
# loop all simulations and do several analyses
##======


gp_dir = '../'
# for reading output from ms
ms_f = open(gp_dir + gp.sim_out_dir + '/ms/' + gp.sim_out_prefix + tag + '.txt', 'r')

fill_symbol = '0'

for i in range(args['num_reps']):
    
    print i

    ##======
    # read in simulated sequences
    ##======
    
    # trees, recomb_sites, seg_sites, positions, seqs
    sim = sim_process.read_one_sim(ms_f, args['num_sites'], args['num_samples'])

    ##======
    # summarize properties of sequences
    ##======

    sim_stats(sim, args)

    ##======
    # calculate frequency of ILS
    ##======




    ##======
    # find introgressed/non-introgressed tracts
    ##======

    # fill in nonpolymorphic sites
    sim['seqs'] = sim_process.fill_seqs(sim['seqs'], sim['positions'], \
                                            args['num_sites'], fill_symbol)



##======
# output
##======


