# goal for plot:
# graph of probabilities with threshold
# predicted
# actual introgression
# reference
# SNPs

# do this by producing file:
# site coding rep predicted actual actual_ref prob_cer prob_par
#    0     ++  0        cer    cer        cer      .99       .1
#    .
#    .
#    .

import sys
import os
from plot_rep_setup import *
import process_args
import sim_process
import sim_predict
sys.path.append('..')
import global_params as gp

##======
# read in simulation parameters
##======

sim_tag = sys.argv[2]
sim_args = process_args.process_args_by_tag(sys.argv[1], sim_tag)
predict_args, last_read = sim_predict.process_args(sys.argv, sim_args, i=2)

predict_prob_block_types = ['predicted_' + predict_args['predict_tag']]
predict_path_block_types = ['predicted_viterbi_' + predict_args['predict_tag']]
actual_block_type = 'actual'
block_types = predict_prob_block_types + predict_path_block_types + [actual_block_type]

##======
# produce combined file
##======

gp_dir = '../'
# for reading output from ms
ms_f = open(gp_dir + gp.sim_out_dir + '/ms/' + gp.sim_out_prefix + \
               sim_args['tag'] + '.txt', 'r')

# all introgressed block files to read
introgression_fn_prefix = gp_dir + gp.sim_out_dir + gp.sim_out_prefix + \
                          sim_args['tag'] + '_introgressed_'
introgression_files = dict(zip(block_types, \
                               [open(introgression_fn_prefix + t + '.txt', 'r') \
                                for t in block_types]))
introgression_file_lines = dict(zip(block_types, \
                                    [introgression_files[k].readline() \
                                     for k in introgression_files]))

# all prob files to read
prob_fn_prefix = gp_dir + gp.sim_out_dir + gp.sim_out_prefix + \
                          sim_args['tag'] + '_introgressed_probs_'
prob_files = dict(zip(predict_prob_block_types, \
                      [open(prob_fn_prefix + t + '.txt', 'r') \
                       for t in predict_prob_block_types]))
prob_file_lines = dict(zip(predict_prob_block_types, \
                           [prob_files[k].readline() for k in prob_files]))


# indices of individuals in species we're predicting introgression in
inds = sim_args['species_to_indices'][sim_args['species_to']]
ref_ind = predict_args['ref_inds'][0]
inds.remove(ref_ind)

# combined output files, one per predicted strain and rep
combined_dir = gp_dir + gp.sim_out_dir + sim_args['tag'] + '/' + \
               predict_args['predict_tag'] + '/'
combined_fn_prefix = combined_dir + gp.sim_out_prefix + \
                     sim_args['tag'] + '_' + predict_args['predict_tag'] + \
                     '_combined_strain_'
if not os.path.exists(combined_dir):
    os.makedirs(combined_dir)

# loop through reps and then individuals
#for i in range(sim_args['num_reps']):
for i in range(100):

    print 'rep', i

    combined_files = dict(zip(inds, \
                              [open(combined_fn_prefix + str(ind) + '_rep' + str(i) + \
                                    '.txt', 'w') \
                               for ind in inds]))

    sim = sim_process.read_one_sim(ms_f, sim_args['num_sites'], sim_args['num_samples'])

    seqs_coded = sim_predict.set_up_seqs(sim, sim_args, predict_args)

    # read in blocks from all methods
    # keyed by individual, then block type, list of species at all sites
    blocks_dic = {} 

    # keyed by individual and then block_type and then species, list of probs
    probs = {}

    for t in block_types:
        # d is keyed by individual, then species
        d, rep, line = sim_process.read_introgression_blocks(\
                            introgression_files[t], \
                            introgression_file_lines[t], \
                            predict_args['states'])
        assert i == rep, str(i) + ' ' + str(rep)
        introgression_file_lines[t] = line
        d = sim_process.unblock(d, sim_args['num_sites'])
        # this is just converting the dictionary to have the block
        # type layer
        for ind in inds:
            if not blocks_dic.has_key(ind):
                blocks_dic[ind] = {}
            blocks_dic[ind][t] = d[ind]

        # want to keep track of (actual) introgression in reference also
        if t == actual_block_type:
            blocks_dic[ref_ind] = {}
            blocks_dic[ref_ind][t] = d[ref_ind]
        # probs only exist for predictions with probabilities
        elif t in predict_prob_block_types:
            probs_ind, rep, line = sim_process.read_state_probs(prob_files[t], \
                                                                prob_file_lines[t], \
                                                                predict_args['states'])
            assert i == rep, str(i) + ' ' + str(rep)
            prob_file_lines[t] = line
            # this is just converting the dictionary to have the block
            # type layer
            for ind in inds:
                if not probs.has_key(ind):
                    probs[ind] = {}
                probs[ind][t] = probs_ind[ind]

    write_combined_files(combined_files, inds, i, seqs_coded, blocks_dic, probs, \
                         actual_block_type, ref_ind, True)

for f in introgression_files.values() + prob_files.values() + combined_files.values():
    f.close()
