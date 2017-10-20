import sys
import os
import process_args
import sim_process
from sim_actual import *
sys.path.append('..')
import global_params as gp

##======
# read in simulation parameters
##======
 
args, last_read = process_args.process_args(sys.argv)

##======
# loop through all simulations and do several analyses
##======

gp_dir = '../'
# for reading output from ms
ms_f = open(gp_dir + gp.sim_out_dir + '/ms/' + gp.sim_out_prefix + \
                args['tag'] + '.txt', 'r')
# summary output
out_f = open(gp_dir + gp.sim_out_dir + gp.sim_out_prefix + \
                args['tag'] + '_summary.txt', 'w')
# introgression output
introgression_f = open(gp_dir + gp.sim_out_dir + gp.sim_out_prefix + \
                           args['tag'] + '_introgressed_actual.txt', 'w')

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

    stats = sim_stats(sim, args)

    ##======
    # calculate frequency of ILS (or of possible ILS...)
    ##======

    concordance_info = calculate_ils(sim, args)

    ##======
    # find introgressed/non-introgressed tracts
    ##======

    introgression_stats, actual_state_seq = find_introgressed(sim, args)
    actual_state_seq_blocks = sim_process.convert_to_blocks(actual_state_seq, \
                                                            args['states'])

    ##======
    # output
    ##======

    # general summary statistics about simulated sequences
    write_output_line(stats, concordance_info, introgression_stats, out_f, i==0) 

    # specific locations of introgression (for comparing predictions
    # to)
    sim_process.write_introgression_blocks(actual_state_seq_blocks, introgression_f, \
                                           i, args['states'])

ms_f.close()
out_f.close()
introgression_f.close()


