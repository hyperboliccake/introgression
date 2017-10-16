import sys
import os
import copy
import itertools
import sim_predict
sys.path.append('..')
import global_params as gp

def convert_binary_to_nucleotides(seqs):
    n = ['A', 'T', 'G', 'C']
    seqs_n = [[] for s in range(len(seqs))]
    for i in range(len(seqs[0])):
        l0 = random.choice(n)
        l1 = random.choice(n)
        for s in range(len(seqs)):
            if seqs[s][i] == '0':
                seqs_n[s].append(l0)
            else:
                assert seqs[s][i] == '1'
                seqs_n[s].append(l1)
    return seqs_n

def predict_introgressed(sim, args):

    # fill in nonpolymorphic sites
    fill_symbol = '0'
    seqs_filled = sim_predict.fill_seqs(sim['seqs'], sim['recomb_sites'], \
                                        args['num_sites'], fill_symbol)

    # use letters because phylo-hmm seems set up only for that
    seqs_filled = convert_binary_to_nucleotides(seqs_filled)

    # code sequences by which references they match at each position
    ref_seqs = [seqs_filled[r] for r in args['ref_inds']]
    seqs_coded = sim_predict.code_seqs(seqs_filled, args['num_sites'], ref_seqs)

    #HERE

    seq_fn = gp_dir + gp.sim_out_dir + '/ms/' + gp.sim_out_prefix + \
             'sequence_' + tag + '_rep' + str(i) + '.fasta'
    # TODO unhardcode
    write_fasta(sim[4], ['C1', 'C2', 'P', 'OUTGROUP'], seq_fn)

    # create input file for phylo-hmm
    input_fn = gp_dir + gp.sim_out_dir + '/phylo-hmm/' + 'autoinput_' + \
        tag + '_rep' + str(i) + '.txt'
    working_dir = gen_input_file(seq_fn, input_fn, tag, i)

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


    # make predictions
    os.system('java -jar ~/software/phylo_hmm/phmm-0.1/dist/lib/phmm.jar < ' + input_fn)


    return predicted, hmm

