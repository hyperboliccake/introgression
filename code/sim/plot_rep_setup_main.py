# goal is to produce a plot showing actual introgressed sequence,
# where strain matches one/neither/both references, and multiple types
# of predictions

# this should be pretty similar to the code for plotting introgressed
# regions across all strains

# produce output here in a way easily parseable by R

# start and end coordinates (contained in below)
# at each site, code for matching some combination of reference sequences
# actual introgressed start and end coordinates, arbitrary number
# predicted (method #1) introgressed start and end coordinates, arbitrary number
# predicted (method #2) introgressed start and end coordinates, arbitrary number

# table 1 (coding): 
# site code rep
# 1    ++-  1
# 2    +++  1
# 3    --+
# 4    +++
# 5    +++

# table 2 (tracts):
# type                species  start  end   rep
# actual                       400    500   1
# actual (reference)           80     120   1
# predicted method 1           350    450   1
# predicted method 1           480    500   1
# predicted method 2           100    200   1
# predicted method 2           480    500   1

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

args, last_read = process_args.process_args(sys.argv)
# TODO fix
block_types = sys.argv[last_read+1+len(args['states']):]

##======
# loop through all simulations predict introgression
##======

gp_dir = '../'
# for reading output from ms
ms_f = open(gp_dir + gp.sim_out_dir + '/ms/' + gp.sim_out_prefix + \
               args['tag'] + '.txt', 'r')
# all introgressed block files to read
introgression_fn_prefix = gp_dir + gp.sim_out_dir + gp.sim_out_prefix + \
                          args['tag'] + '_introgressed_'
introgression_files = [open(introgression_fn_prefix + t + '.txt', 'r') \
                       for t in block_types]
introgression_file_lines = [f.readline() for f in introgression_files]

# indices of individuals in species with introgression (keep the
# reference in even though it's not interesting
inds = args['species_to_indices'][args['species_to']]
#inds.remove(args['ref_inds'][0])

# table 1 (coding of which references each strain matches at each
# site, 1 file per strain)
coding_fn_prefix = gp_dir + gp.sim_out_dir + gp.sim_out_prefix + \
               args['tag'] + '_site_codings_strain_'
coding_files = dict(zip(inds, \
                        [open(coding_fn_prefix + str(i) + '.txt', 'w') \
                         for i in inds]))


# table 2 (actual and predicted introgressed blocks, 1 file per
# strain)
blocks_fn_prefix = gp_dir + gp.sim_out_dir + gp.sim_out_prefix + \
               args['tag'] + '_introgressed_blocks_strain_'
blocks_files = dict(zip(inds, \
                        [open(blocks_fn_prefix + str(i) + '.txt', 'w') \
                         for i in inds]))

# loop through reps and then individuals
for i in range(args['num_reps']):

    sim = sim_process.read_one_sim(ms_f, args['num_sites'], args['num_samples'])

    seqs_coded = sim_predict.set_up_seqs(sim, args)
    write_coding_table(seqs_coded, coding_files, i, i == 0)

    # read in blocks from all methods
    # keyed by individual, then block type, then species
    blocks_dic = {} 
    for j in range(len(block_types)):
        # d is keyed by individual, then species
        d, rep, line = sim_process.read_introgression_blocks(introgression_files[j], \
                                                        introgression_file_lines[j], \
                                                        args['states'])
        assert i == rep, str(i) + ' ' + str(rep)
        introgression_file_lines[j] = line

        # this is just converting the dictionary to have the block
        # type layer
        for ind in d.keys():
            if not blocks_dic.has_key(ind):
                blocks_dic[ind] = {}
            if not blocks_dic[ind].has_key(block_types[j]):
                blocks_dic[ind][block_types[j]] = {}
            for species in d[ind]:
                if not blocks_dic[ind][block_types[j]].has_key(species):
                    blocks_dic[ind][block_types[j]][species] = []
                blocks_dic[ind][block_types[j]][species] += d[ind][species]

    write_blocks_table(blocks_dic, blocks_files, i, i == 0)

for f in introgression_files + coding_files.values() + blocks_files.values():
    f.close()
