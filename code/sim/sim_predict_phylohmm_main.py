import sys
import os
import process_args
import sim_process
import sim_predict_phylohmm 
import gzip
sys.path.append('..')
import global_params as gp
sys.path.append('../analyze')
import predict

##======
# read in simulation parameters
##======

sim_tag = sys.argv[2]
sim_args = process_args.process_args_by_tag(sys.argv[1], sim_tag)
predict_args, last_read = sim_predict_phylohmm.process_args(sys.argv, sim_args, i=2)
write_ps = True

##======
# loop through all simulations predict introgression
##======

gp_dir = '../'
# for reading output from ms
ms_f = open(gp.sim_out_dir_absolute + '/ms/' + gp.sim_out_prefix + \
            sim_tag + '.txt', 'r')
# summary output
out_f = open(gp.sim_out_dir_absolute + gp.sim_out_prefix + \
             sim_tag + '_phylohmm_' + \
             predict_args['predict_tag'] + '.txt', 'w')
# introgression output
introgression_f = open(gp.sim_out_dir_absolute + gp.sim_out_prefix + \
                       sim_tag + '_introgressed_predicted_phylohmm_' + \
                       predict_args['predict_tag'] + '.txt', 'w')
# associated probabilities output
prob_f = open(gp.sim_out_dir_absolute + gp.sim_out_prefix + \
              sim_tag + '_introgressed_probs_predicted_phylohmm_' + \
              predict_args['predict_tag'] + '.txt', 'w')

# sites actually used in making predicitons (polymorphic sites)
ps_f = gzip.open(gp.sim_out_dir_absolute + gp.sim_out_prefix + \
                 sim_tag + '_positions_predicted_phylohmm_' + \
                 predict_args['predict_tag'] + '.txt.gz', 'wb')

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

    state_seq, probs, init, emis, trans, ps = \
        sim_predict_phylohmm.predict_introgressed(sim, sim_args, predict_args, i, gp_dir)

    state_seq_blocks = sim_process.convert_to_blocks(state_seq, \
                                                     sim_args['species'])

    ##======
    # output
    ##======

    if write_ps:
        predict.write_positions(ps, ps_f, str(i), 'I')

    # summary info about HMM
    # sim_predict.write_hmm_line(hmm, out_f, i==0) 

    # locations of introgression
    sim_process.write_introgression_blocks(state_seq_blocks, introgression_f, \
                                           i, sim_args['species'])

    # probabilities at each site
    sim_process.write_state_probs(probs, prob_f, i)

ms_f.close()
out_f.close()
introgression_f.close()
prob_f.close()
ps_f.close()
