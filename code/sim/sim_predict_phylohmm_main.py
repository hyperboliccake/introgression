import sys
import os
import process_args
import sim_process
import sim_predict
from sim_predict_phylohmm import *
sys.path.append('..')
import global_params as gp

##======
# read in simulation parameters
##======

sim_tag = sys.argv[2]
sim_args = process_args.process_args_by_tag(sys.argv[1], sim_tag)
predict_args, last_read = process_args(sys.argv, sim_args, i=2)

##======
# loop through all simulations predict introgression
##======

gp_dir = '../'
# for reading output from ms
ms_f = open(gp_dir + gp.sim_out_dir + '/ms/' + gp.sim_out_prefix + \
            predict_args['tag'] + '.txt', 'r')
# summary output
out_f = open(gp_dir + gp.sim_out_dir + gp.sim_out_prefix + \
             sim_args['tag'] + '_phylohmm_' + \
             predict_args['predict_tag'] + '.txt', 'w')
# introgression output
introgression_f = open(gp_dir + gp.sim_out_dir + gp.sim_out_prefix + \
                       sim_args['tag'] + '_introgressed_predicted_phylohmm_' + \
                       predict_args['predict_tag'] + '.txt', 'w')

for i in range(sim_args['num_reps']):
    
    print i

    ##======
    # read in simulated sequences
    ##======
    
    # trees, recomb_sites, seg_sites, positions, seqs
    sim = sim_process.read_one_sim(ms_f, sim_args['num_sites'], sim_args['num_samples'])

    ##======
    # predict introgressed/non-introgressed tracts
    ##======

    state_seq, probs, init, emis, trans = \
        predict_introgressed(sim, sim_args, predict_args, i, gp_dir)

    state_seq_blocks = sim_process.convert_to_blocks(state_seq, \
                                                     sim_args['species'])

    ##======
    # output
    ##======

    # summary info about HMM
    # sim_predict.write_hmm_line(hmm, out_f, i==0) 

    # specific locations of introgression (for comparing predictions
    # to)
    sim_process.write_introgression_blocks(state_seq_blocks, introgression_f, \
                                           i, sim_args['species'])

ms_f.close()
out_f.close()
introgression_f.close()

