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
block_types = sys.argv[last_read+1:]

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
introgression_files = [open(introgression_fn_prefix + t + '.txt', 'r') \
                       for t in block_types]
introgression_file_lines = [f.readline() for f in introgression_files]

# all prob files to read
prob_fn_prefix = gp_dir + gp.sim_out_dir + gp.sim_out_prefix + \
                          sim_args['tag'] + '_introgressed_probs_'
prob_files = [open(prob_fn_prefix + t + '_' + \
                   predict_args['predict_tag'] + '.txt', 'r') \
              for t in block_types]
prob_file_lines = [f.readline() for f in prob_files]


# indices of individuals in species we're predicting introgression in
inds = sim_args['species_to_indices'][sim_args['species_to']]
for i in inds:
    if i in predict_args['ref_inds']:
        inds.remove(i)

# combined output files, one per predicted strain
combined_fn_prefix = gp_dir + gp.sim_out_dir + gp.sim_out_prefix + \
                     sim_args['tag'] + '_site_codings_strain_'
combined_files = dict(zip(inds, \
                          [open(coding_fn_prefix + str(i) + '.txt', 'w') \
                           for i in inds]))


# loop through reps and then individuals
for i in range(sim_args['num_reps']):

    sim = sim_process.read_one_sim(ms_f, sim_args['num_sites'], sim_args['num_samples'])

    seqs_coded = sim_predict.set_up_seqs(sim, sim_args, predict_args)

    # read in blocks from all methods
    # keyed by individual, then block type, list of species at all sites
    blocks_dic = {} 

    # keyed by individual and then block_type and then species, list of probs
    probs = {}

    for j in range(len(block_types)):
        # d is keyed by individual, then species
        d, rep, line = sim_process.read_introgression_blocks(introgression_files[j], \
                                                        introgression_file_lines[j], \
                                                        predict_args['states'])
        assert i == rep, str(i) + ' ' + str(rep)
        introgression_file_lines[j] = line
        d = sim_process.unblock(d, sim_args['num_sites'])

        probs_ind, rep, line = sim_process.read_state_probs(prob_files[j], \
                                                            prob_file_lines[j], \
                                                            predict_args['states'])
        assert i == rep, str(i) + ' ' + str(rep)
        prob_file_lines[j] = line

        # this is just converting the dictionaries to have the block
        # type layer
        for ind in d.keys():
            if not blocks_dic.has_key(ind):
                blocks_dic[ind] = {}
            blocks_dic[ind][block_types[j]] = d[ind]
            if not probs.has_key(ind):
                probs[ind] = {}
            probs[ind][block_types[j]] = probs[ind]

    write_combined_files(combined_files, i, seqs_coded, blocks_dic, probs, i==0)

for f in introgression_files + prob_files + combined_files.values():
    f.close()
