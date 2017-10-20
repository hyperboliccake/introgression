import sys
import os
import process_args
import sim_process
from sim_predict import *
sys.path.append('..')
import global_params as gp

##======
# read in simulation parameters
##======

args, last_read = process_args.process_args(sys.argv)

##======
# loop through all simulations predict introgression
##======

gp_dir = '../'
# for reading output from ms
ms_f = open(gp_dir + gp.sim_out_dir + '/ms/' + gp.sim_out_prefix + \
                args['tag'] + '.txt', 'r')
# summary output
out_f = open(gp_dir + gp.sim_out_dir + gp.sim_out_prefix + \
                args['tag'] + '_hmm.txt', 'w')
# introgression output
introgression_f = open(gp_dir + gp.sim_out_dir + gp.sim_out_prefix + \
                           args['tag'] + '_introgressed_predicted.txt', 'w')

for i in range(args['num_reps']):
    
    print i

    ##======
    # read in simulated sequences
    ##======
    
    # trees, recomb_sites, seg_sites, positions, seqs
    sim = sim_process.read_one_sim(ms_f, args['num_sites'], args['num_samples'])

    ##======
    # predict introgressed/non-introgressed tracts
    ##======
    
    state_seq, hmm = predict_introgressed(sim, args, train=True)
    state_seq_blocks = sim_process.convert_to_blocks(state_seq, \
                                                     args['states'])

    ##======
    # output
    ##======

    # summary info about HMM
    write_hmm_line(hmm, out_f, i==0) 

    # specific locations of introgression (for comparing predictions
    # to)
    sim_process.write_introgression_blocks(state_seq_blocks, introgression_f, \
                                           i, args['states'])
    
ms_f.close()
out_f.close()
introgression_f.close()

