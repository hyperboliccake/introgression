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

args = process_args.process_args(sys.argv)

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
    
        # use letters because phylo-hmm seems set up only for that
        sim[4] = convert_binary_to_nucleotides(seqs_filled)
        # also write to fasta
        seq_fn = gp_dir + gp.sim_out_dir + '/ms/' + gp.sim_out_prefix + \
            'sequence_' + tag + '_rep' + str(i) + '.fasta'
        # TODO unhardcode
        write_fasta(sim[4], ['C1', 'C2', 'P', 'OUTGROUP'], seq_fn)

    # create input file for phylo-hmm
    input_fn = gp_dir + gp.sim_out_dir + '/phylo-hmm/' + 'autoinput_' + \
        tag + '_rep' + str(i) + '.txt'
    working_dir = gen_input_file(seq_fn, input_fn, tag, i)

    # run phylo-hmm
    os.system('java -jar ~/software/phylo_hmm/phmm-0.1/dist/lib/phmm.jar < ' + input_fn)

    # write results in different format
    trees_to_states = {'p1':'cer', 'p2':'par'} # generalize this? worth it? nah
    init_new, emis_new, trans_new = process_phylo_output(sim, ref_inds, output_dic, \
                                                             f_out, \
                                                             trees_to_states, tag, i, \
                                                             num_samples_species_to, \
                                                             topology, species_to, \
                                                             index_to_species, states, \
                                                             f_tracts_predicted, \
                                                             working_dir + \
                                                             '/filtered_sites.txt')
#############

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

