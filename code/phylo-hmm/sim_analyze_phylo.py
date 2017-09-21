import os
import sys
import copy
import itertools
import random
sys.path.append('../sim')
from concordance_functions import *
import sim_analyze_hmm_bw
sys.path.append('../hmm')
from hmm_bw import *
sys.path.append('..')
import global_params as gp

def read_predicted(fn, trees_to_states):
    
    # TODO make this deal with filtered sites? here or somewhere else?

    f = open(fn, 'r')
    line = f.readline()
    predicted_p = []
    predicted_g = []
    predicted = []
    while line != '':
        line = line.strip().split(',')
        # mebbe do something with these later
        predicted_p.append(line[0])
        predicted_g.append(line[1])
        predicted.append(trees_to_states[line[0]])
        # TODO figure out why 2 columns of genealogies (probably for diploids?)
        assert line[1] == line[2], line[1] + ' ' + line[2]
        line = f.readline()

    f.close()
    return predicted

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

def write_fasta(seqs, names, fn):
    f = open(fn, 'w')
    for i in range(len(names)):
        f.write('> ' + names[i] + '\n')
        f.write(''.join(seqs[i]) + '\n')
    f.close()

def get_actual(num_samples_species_to, trees, topology, species_to, index_to_species, recomb_sites):

    # sequence of states, one for each site and strain
    actual_state_seq = [[] for i in range(num_samples_species_to)]

    # and how many bases are introgressed in total in each strain
    num_introgressed = [0] * num_samples_species_to

    # loop through the trees for all blocks with no recombination
    # within them
    for ti in range(len(trees)):

        # note that species indices/labels are shifted to start at
        # 0 instead of 1
        t = trees[ti]

        # identify sequences that are introgressed from the one or
        # two other species, based on coalescent tree; could clean
        # this up a little

        # only supports two species 

        # introgressed is a list of species (one entry for
        # each individual in to species)
        introgressed, num_lineages_at_join_current = \
            sim_analyze_hmm_bw.find_introgressed_3(t, species_to, topology, index_to_species)

        # number of sites in the current block of sequence
        num_sites_t = recomb_sites[ti]

        # for all strains that have this block introgressed, add
        # the length of the block to the total number of
        # introgressed sites across all strains; also update the
        # state sequence
        for i in range(num_samples_species_to):
            if introgressed[i] != species_to:
                num_introgressed[i] += num_sites_t
            actual_state_seq[i] += [introgressed[i]] * num_sites_t
            
    return actual_state_seq, num_introgressed


def process_phylo_output(sim, ref_inds, output_dic, f_out, trees_to_states, tag, rep, \
                             num_samples_species_to, \
                             topology, species_to, index_to_species, \
                             states, f_tracts_predicted, filtered_sites_fn):

    # gah so much copy pasta here

    trees, recomb_sites, segsites, positions, seqs = sim

    # read predicted state sequence
    viterbi_fn = '../../results/sim/phylo-hmm/optimized.viterbi.sequence.'  + \
        tag + '.' + str(rep)
    predicted = read_predicted(viterbi_fn, trees_to_states)

    # actual state sequence
    actual_state_seq_unfiltered, num_introgressed = get_actual(num_samples_species_to, \
                                                                   trees, topology, \
                                                                   species_to, \
                                                                   index_to_species, recomb_sites)

    try:
        filtered_sites = [int(x) for x in \
                              open(filtered_sites_fn, 'r').readline().strip().split(' ')]
    except:
        print 'looks like none of the sites passed filtering'
        sys.exit()
    actual_state_seq = [[] for s in range(num_samples_species_to)]
    for i in filtered_sites:
        for s in range(num_samples_species_to):
            actual_state_seq[s].append(actual_state_seq_unfiltered[s][i])
    os.system('mv ' + filtered_sites_fn + ' ../' + gp.sim_out_dir + '/phylo-hmm/' + \
                  'filtered_sites_' + tag + '_rep' + str(rep) + '.txt')

    # TODO fix this stupidness; C1 is the reference, C2 what we're trying to predict
    predicted = [predicted]
    actual_state_seq = [actual_state_seq[1]]

    # see how well predictions did
    num_correct, num_introgressed_correct, actual_lens, predicted_lens, \
        num_predicted_tracts_actual, num_actual_tracts_predicted, \
        num_introgressed_tracts, num_not_introgressed_tracts, \
        num_predicted_introgressed_tracts, num_predicted_not_introgressed_tracts = \
        sim_analyze_hmm_bw.evaluate_predicted(predicted, actual_state_seq, species_to)

    groups = sim_analyze_hmm_bw.group_actual_predicted_bases(actual_state_seq, predicted, states)
    for group in groups:
        output_dic = sim_analyze_hmm_bw.update_value(output_dic, 'num_bases_actual_' + group[0] + \
                                                         '_predicted_' + group[1], \
                                                         groups[group])

    ### tracts
    inds_to_predict = range(num_samples_species_to)
    inds_to_predict.remove(ref_inds[0])
    blocks_predicted, blocks_actual = \
        sim_analyze_hmm_bw.evaluate_predicted_blocks(predicted, actual_state_seq, \
                                                         species_to, states, inds_to_predict)
    # remove blocks in reference individual
    blocks_actual = filter(lambda b: b[3] != ref_inds[0], blocks_actual)
    d_actual_predicted, d_predicted_actual, d_actual_counts, d_predicted_counts = \
        sim_analyze_hmm_bw.group_actual_predicted_blocks(blocks_actual, blocks_predicted, states, \
                                                             inds_to_predict)
    for group in d_actual_predicted:
        output_dic = sim_analyze_hmm_bw.update_value(output_dic, 'num_tracts_actual_' + group[0] + \
                                                         '_predicted_' + group[1], \
                                                         d_actual_predicted[group])
    for group in d_predicted_actual:
        output_dic = sim_analyze_hmm_bw.update_value(output_dic, 'num_tracts_predicted_' + group[0] + \
                                                         '_actual_' + group[1], \
                                                         d_predicted_actual[group])
    for state in states:
        output_dic = sim_analyze_hmm_bw.update_value(output_dic, 'num_tracts_actual_' + state, \
                                                         d_actual_counts[state])
        output_dic = sim_analyze_hmm_bw.update_value(output_dic, 'num_tracts_predicted_' + state, \
                                                         d_predicted_counts[state])

    # tract lengths, actual and predicted; list, not average
    for state in states:
        tract_lengths = [b[2] for b in filter(lambda x: x[0] == state, blocks_actual)]
        output_dic = sim_analyze_hmm_bw.update_value(output_dic, 'tract_lengths_actual_' + state, tract_lengths)

        tract_lengths = [b[2] for b in filter(lambda x: x[0] == state, blocks_predicted)]
        output_dic = sim_analyze_hmm_bw.update_value(output_dic, 'tract_lengths_predicted_' + state, tract_lengths)


    """
    # HMM parameters
    init, emis, trans = read_phylo_hmm_params()

    for i in range(len(states)):
        output_dic = sim_analyze_hmm_bw.update_value(output_dic, 'init_' + states[i], init[i])
            
    for i in range(len(states)):
        for j in range(len(states)):
            output_dic = sim_analyze_hmm_bw.update_value(output_dic, \
                                                             'trans_' + states[i] + '_' + states[j], \
                                                             trans[i][j])

    # so here emis_cer_par is going to be the probability that cer
    # state emits symbol that matches par, so *+; note that these
    # don't have to add to 1 (states includes unknown if
    # appropriate)
    for i in range(len(states)):
        to_states = states
        for j in range(len(states)):
            total = 0
            for s in emis[i].keys():
                # for unknown state, include --- (as well as '?')
                if states[j] == unknown_species:
                    if gp.match_symbol not in s:
                        total += emis[i][s]
                elif s[j] == gp.match_symbol:
                    total += emis[i][s]
            output_dic = update_value(output_dic, 'emis_'  + states[i] + '_' + \
                                          to_states[j], total)
                                          """
    #####
    # write results to file
    #####

    sim_analyze_hmm_bw.write_output_line(f_out, output_dic, False)
    f_out.flush()
    sys.stdout.flush()

    for b in blocks_predicted:
        f_tracts_predicted.write(str(rep) + ' ' +  ' '.join([str(x) for x in b]) + '\n')

    f_tracts_predicted.flush()

    #return init, emis, trans
    return None, None, None # TODO

def gen_input_file(sequence_fn, fn, tag, rep):
    """
    Initial mode:
    0) Build a new model.
    1) Load a pre-existing model.
    2) Exit.
    Choose an option: 
    0
    Input the basic file info path name:
    (note: see README for file format) 
    
    basic-info.txt
    Input the number of trees or states for this HMM:
    6
    
    Input the parental trees file path name:
    (note: see README for file format) 

    parental.trees
    
    Input the gene genealogies file path name:
    (note: see README for file format) 
    
    gene.trees
    
    Input outgroup taxon name, or empty string for no outgroup taxon: 
    
    OUTGROUP
    Empty working directory: 
    ./d
    Input non-zero substitution model rates in format <AG> <AC> <AT> <GC> <GT>: 
    0.964011 0.192217 0.270592 0.0806232 0.174283
    Input non-zero substitution model base frequencies in format <A> <G> <C> <T>: 
    0.1953 0.3077 0.3008 0.1962
    
    Across-row switching frequency gamma: 
    .1
    Hidden state switching frequency ratio term input file: 
    switching-frequency-ratio-terms

    Operate mode:
    0) Run Viterbi
    1) Learn Model using Baum Welch.
    3) Learn Model using a multivariate optimization heuristic that incorporates Brent's method
    4) Exit
    Choose an option:
    0

    Path to your output file: 
    initial.viterbi.sequence
    Which observation sequence would you like to use?
    0) Reuse previously read sequence.
    1) Load a new observation sequence.
    Choose an option:
    1
    Keep or discard parsimony-uninformative sites? true for keep, false for discard: 
    false
    Input observation sequence file path name : 
    sequence_nogaps.fasta
    Input alignment length: 86706
    Alignment length after optional filter step: 8194
    Building Psy Array
    100% Done!                             

    Total: 29 MiB; Free: 22 MiB; --> Used: 7 MiB.
    Max: 454 MiB.
    Begin computing Viterbi's algorithm : FORWARD PART 
    100% DONE!                                
    Enlarging... capacity = 10000
    Begin computing Viterbi's algorithm : BACKWARD PART 
    100% DONE!                                   
    Input HMM Viterbi log likelihood: |-21.866065387121736|
    Computing input HMM log likelihood for input sequences... 
    Computing input HMM log likelihood for input sequences DONE.
    Input HMM log likelihood: |-21.310316817591303|
    Operate mode:
    0) Run Viterbi
    1) Learn Model using Baum Welch.
    3) Learn Model using a multivariate optimization heuristic that incorporates Brent's method
    4) Exit
    Choose an option: 
    

    Which observation sequence would you like to use?
    0) Reuse previously read sequence.
    1) Load a new observation sequence.
    Choose an option: 
    Input parental-branch-length-parameter-to-edge map filename: 
    Input parental-branch-length-parameter strict inequalities filename: 
    Input length-parameter-set-constraints filename: 
    Input checkpoint file to restore from, or empty line for no restore: 
    Output posterior decoding probabilities file: 
    Output Viterbi-optimal hidden state sequence file: 
    Output model likelihoods file: 
    Output file with optimized model parameter values: 
    Initial search settings vector <setting 1> <setting 2> ... <setting s>, where <setting s> is one of CURRENT, DEFAULT, RANDOM
    
    Input parental-branch-length-parameter-to-edge map filename: 
    Input parental-branch-length-parameter strict inequalities filename: 
    Input length-parameter-set-constraints filename: 
    Input checkpoint file to restore from, or empty line for no restore: 
    Output posterior decoding probabilities file: 
    Output Viterbi-optimal hidden state sequence file: 
    Output model likelihoods file: 
    Output file with optimized model parameter values: 
    Initial search settings vector <setting 1> <setting 2> ... <setting s>, where <setting s> is one of CURRENT, DEFAULT, RANDOM
    Enable optimization flag vector <enable parental tree optimization flag> <enable gene genealogy optimization flag> <enable switching frequency optimization flag> <enable substitution model optimization flag>
    0) Run Viterbi
    1) Learn Model using Baum Welch.
    3) Learn Model using a multivariate optimization heuristic that incorporates Brent's method
    4) Exit
    Choose an option: 
    """

    initial_mode = '0' # build new model
    basic_info_fn = 'basic-info.txt'
    num_states = '6' # num parental trees * num genealogies
    parental_trees_fn = 'parental.trees'
    gene_trees_fn = 'gene.trees'
    outgroup_name = 'OUTGROUP'
    working_dir = '../../results/sim/phylo-hmm/working/' + tag # don't need one for every rep because we do one at a time
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)
    substitution_rates = '1 1 1 1 1' # <AG> <AC> <AT> <GC> <GT>
    base_frequencies = '0.25 0.25 0.25 0.25' # <A> <G> <C> <T>
    parental_tree_switching_freq = '.1'
    gene_tree_switching_fn = 'switching-frequency-ratio-terms'
    operate_mode = '0' # run viterbi
    # operate_mode = '1' # learn with Baum-Welch
    output_file_path = '../../results/sim/phylo-hmm/initial.viterbi.sequence.' + tag + '.' + str(rep)
    observation_sequence_option = '1' # read new sequence
    keep_uninformative_sites = 'false' # referred to as "optional filter step" later
    #sequence_fn = 'sequence.fasta'
    operate_mode_2 = '3' # learn with "a multivariate optimization heuristic that incorporates Brent's method"
    observation_sequence_option_2 = '0' # reuse previously read sequence
    length_params_fn = 'length-parameters'
    length_params_inequality_constraints_fn = 'length-parameter-inequality-constraints'
    length_params_constraint_sets_fn = 'length-parameter-constraint-sets'
    restore_fn = '' # blank for no file
    output_posterior_decoding_fn = '../../results/sim/phylo-hmm/optimized.posterior.decoding.probabilities.' + tag + '.' + str(rep)
    output_viterbi_optimized_fn = '../../results/sim/phylo-hmm/optimized.viterbi.sequence.'  + tag + '.' + str(rep)
    output_model_likelihoods_fn = '../../results/sim/phylo-hmm/optimized.model.likelihoods.'  + tag + '.' + str(rep)
    output_optimized_params_fn = '../../results/sim/phylo-hmm/optimized.model.parameters.' + tag + '.' + str(rep)
    initial_search_settings = 'CURRENT DEFAULT'
    enable_optimization = 'true true true true' # <enable parental tree optimization flag> <enable gene genealogy optimization flag> <enable switching frequency optimization flag> <enable substitution model optimization flag>
    operate_mode_3='4' # exit

    f = open(fn, 'w')
    f.write(initial_mode + '\n' + \
                basic_info_fn + '\n' + \
                num_states + '\n' + \
                parental_trees_fn + '\n' + \
                gene_trees_fn + '\n' + \
                outgroup_name + '\n' + \
                working_dir + '\n' + \
                substitution_rates + '\n' + \
                base_frequencies + '\n' + \
                parental_tree_switching_freq + '\n' + \
                gene_tree_switching_fn + '\n' + \
                operate_mode + '\n' + \
                output_file_path + '\n' + \
                observation_sequence_option + '\n' + \
                keep_uninformative_sites + '\n' + \
                sequence_fn + '\n' + \
                operate_mode_2 + '\n' + \
                observation_sequence_option_2 + '\n' + \
                length_params_fn + '\n' + \
                length_params_inequality_constraints_fn + '\n' + \
                length_params_constraint_sets_fn + '\n' + \
                restore_fn + '\n' + \
                output_posterior_decoding_fn + '\n' + \
                output_viterbi_optimized_fn + '\n' + \
                output_model_likelihoods_fn + '\n' + \
                output_optimized_params_fn + '\n' + \
                initial_search_settings + '\n' + \
                enable_optimization + '\n' + \
                operate_mode_3)

    return working_dir
