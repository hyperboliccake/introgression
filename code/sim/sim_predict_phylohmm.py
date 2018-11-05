import sys
import os
import copy
import itertools
import random
import sim_process
import sim_predict
sys.path.append('..')
import global_params as gp

def process_args(arg_list, all_sim_args, i=1):
    
    d = {}

    d['tag'] = arg_list[i]
    i += 1

    d['predict_tag'] = arg_list[i]
    i += 1

    d['threshold'] = float(arg_list[i])

    # TODO unhardcode
    d['states'] = ['cer', 'par']

    return d, i

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

def read_predicted_viterbi(fn, trees_to_states):
    
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

def process_phylo_viterbi_output(trees_to_states, tag, rep, filtered_sites_fn):

    # read predicted state sequence
    viterbi_fn = gp.sim_out_dir_absolute + '/phylo-hmm/optimized.viterbi.sequence.'  + \
        tag + '.' + str(rep)
    predicted = read_predicted_viterbi(viterbi_fn, trees_to_states)

    # move filtered sites file to appopriate output directory
    try:
        filtered_sites = [int(x) for x in \
                          open(filtered_sites_fn, 'r').readline().strip().split(' ')]
    except:
        print 'looks like none of the sites passed filtering'
        sys.exit()

    os.system('mv ' + filtered_sites_fn + ' ' + gp.sim_out_dir_absolute + '/phylo-hmm/' + \
                  'filtered_sites_' + tag + '_rep' + str(rep) + '.txt')

    # TODO implement getting hmm params
    return predicted, None, None, None

def read_predicted_posterior_decoding(fn, state_index_to_species, \
                                      predict_args,
                                      default_state):
    
    # TODO make this deal with filtered sites? here or somewhere else?

    f = open(fn, 'r')
    lines = [x[:-1].split(' ') for x in f.readlines()]
    f.close()

    num_sites = int(lines[-1][0]) + 1
    all_probs = [{} for i in range(num_sites)]
    for line in lines:
        site = int(line[0])
        state = state_index_to_species[int(line[1])]
        prob = float(line[2])
        if not all_probs[site].has_key(state):
            all_probs[site][state] = 0
        all_probs[site][state] += prob
        
    # not actually going to return probs, at least for now, since we want
    # to keep track of the probabilities for all states at each point
    predicted, probs = sim_process.get_max_path(all_probs)

    predicted = sim_process.threshold_predicted(predicted, probs, \
                                                predict_args['threshold'], \
                                                default_state)

    return predicted, all_probs

def process_phylo_posterior_decoding_output(state_index_to_species, tag, rep, \
                                            filtered_sites_fn, \
                                            predict_args, default_state):

    # read predicted state sequence
    posterior_decoding_fn = \
        gp.sim_out_dir_absolute + '/phylo-hmm/optimized.posterior.decoding.probabilities.' + \
        tag + '.' + str(rep)
    predicted, all_probs = \
        read_predicted_posterior_decoding(posterior_decoding_fn, \
                                          state_index_to_species, \
                                          predict_args, default_state)

    ps = [int(x) for x in open(filtered_sites_fn, 'r').readline().strip().split(' ')]

    # TODO implement getting hmm params
    return predicted, all_probs, None, None, None, ps


def predict_introgressed(sim, sim_args, predict_args, i, gp_dir):

    # fill in nonpolymorphic sites
    fill_symbol = '0'
    seqs_filled = sim_process.fill_seqs(sim['seqs'], sim['positions'], \
                                        sim_args['num_sites'], fill_symbol)

    # use letters because phylo-hmm seems set up only for that
    seqs_filled = convert_binary_to_nucleotides(seqs_filled)

    # and write to file
    seq_fn = gp.sim_out_dir_absolute + '/ms/' + gp.sim_out_prefix + \
             'sequence_' + sim_args['tag'] + '_rep' + str(i) + '.fasta'
    write_fasta(seqs_filled, ['C1', 'C2', 'P', 'OUTGROUP'], seq_fn) # TODO unhardcode

    # create input file for phylo-hmm
    input_fn = gp.sim_out_dir_absolute + '/phylo-hmm/' + 'autoinput_' + \
        sim_args['tag'] + '_rep' + str(i) + '.txt'
    working_dir = gen_input_file(seq_fn, input_fn, sim_args['tag'], i)

    # make predictions
    phylohmm_command = \
            'java -jar /tigress/anneec/software/phylo_hmm/phmm-0.1/dist/lib/phmm.jar < ' + input_fn
    print phylohmm_command
    os.system(phylohmm_command)

    # write results in different format
    """
    trees_to_states = {'p1':'cer', 'p2':'par'} # generalize this? worth it? nah
    state_seq, init, emis, trans = process_phylo_viterbi_output(trees_to_states, \
                                                        args['tag'], \
                                                        i, \
                                                        working_dir + \
                                                        '/filtered_sites.txt')
    """
    state_index_to_species = {0:'cer',1:'cer',2:'cer',3:'par',4:'par',5:'par'} 
    default_state = sim_args['species_to']
    state_seq, probs, init, emis, trans, ps = \
        process_phylo_posterior_decoding_output(state_index_to_species, \
                                                sim_args['tag'], i, \
                                                working_dir + '/filtered_sites.txt',
                                                predict_args, default_state)
    state_seq = sim_predict.fill_prediction(state_seq, ps, 0, \
                                            sim_args['num_sites'] - 1, \
                                            predict_args['states'])

    # TODO gah
    state_seq_dic = {'1': state_seq}
    probs_dic = {'1': probs}

    return state_seq_dic, probs_dic, init, emis, trans, ps

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
    basic_info_fn = 'phylohmm_inputs/basic-info.txt'
    num_states = '6' # num parental trees * num genealogies
    parental_trees_fn = 'phylohmm_inputs/parental.trees'
    gene_trees_fn = 'phylohmm_inputs/gene.trees'
    outgroup_name = 'OUTGROUP'
    working_dir = gp.sim_out_dir_absolute + '/phylo-hmm/working/' + tag # don't need one for every rep because we do one at a time
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    substitution_rates = '1 1 1 1 1' # <AG> <AC> <AT> <GC> <GT>
    base_frequencies = '0.25 0.25 0.25 0.25' # <A> <G> <C> <T>
    parental_tree_switching_freq = '.1'
    gene_tree_switching_fn = 'phylohmm_inputs/switching-frequency-ratio-terms'
    operate_mode = '0' # run viterbi
    # operate_mode = '1' # learn with Baum-Welch
    output_file_path = gp.sim_out_dir_absolute + '/phylo-hmm/initial.viterbi.sequence.' + tag + '.' + str(rep)
    observation_sequence_option = '1' # read new sequence
    keep_uninformative_sites = 'false' # referred to as "optional filter step" later
    #sequence_fn = 'sequence.fasta'
    operate_mode_2 = '3' # learn with "a multivariate optimization heuristic that incorporates Brent's method"
    observation_sequence_option_2 = '0' # reuse previously read sequence
    length_params_fn = 'phylohmm_inputs/length-parameters'
    length_params_inequality_constraints_fn = 'phylohmm_inputs/length-parameter-inequality-constraints'
    length_params_constraint_sets_fn = 'phylohmm_inputs/length-parameter-constraint-sets'
    restore_fn = '' # blank for no file
    output_posterior_decoding_fn = gp.sim_out_dir_absolute + '/phylo-hmm/optimized.posterior.decoding.probabilities.' + tag + '.' + str(rep)
    output_viterbi_optimized_fn = gp.sim_out_dir_absolute + '/phylo-hmm/optimized.viterbi.sequence.'  + tag + '.' + str(rep)
    output_model_likelihoods_fn = gp.sim_out_dir_absolute + '/phylo-hmm/optimized.model.likelihoods.'  + tag + '.' + str(rep)
    output_optimized_params_fn = gp.sim_out_dir_absolute + '/phylo-hmm/optimized.model.parameters.' + tag + '.' + str(rep)
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
