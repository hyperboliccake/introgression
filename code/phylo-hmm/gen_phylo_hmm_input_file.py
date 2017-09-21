import os
import sys
sys.path.insert(0, '../sim')
import sim_analyze_hmm_bw
sys.path.insert(0, '..')
import global_params as gp


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
    os.mkdir(working_dir)
    substitution_rates = '1 1 1 1 1' # <AG> <AC> <AT> <GC> <GT>
    base_frequencies = '0.25 0.25 0.25 0.25' # <A> <G> <C> <T>
    parental_tree_switching_freq = '.1'
    gene_tree_switching_fn = 'switching-frequency-ratio-terms'
    operate_mode = '0' # run viterbi
    # operate_mode = '1' # learn with Baum-Welch
    output_file_path = '../../results/sim/phylo-hmm/initial.viterbi.sequence.' + tag + '.' + rep
    observation_sequence_option = '1' # read new sequence
    keep_uninformative_sites = 'true' # referred to as "optional filter step" later
    #sequence_fn = 'sequence.fasta'
    operate_mode_2 = '3' # learn with "a multivariate optimization heuristic that incorporates Brent's method"
    observation_sequence_option_2 = '0' # reuse previously read sequence
    length_params_fn = 'length-parameters'
    length_params_inequality_constraints_fn = 'length-parameter-inequality-constraints'
    length_params_constraint_sets_fn = 'length-parameter-constraint-sets'
    restore_fn = '' # blank for no file
    output_posterior_decoding_fn = '../../results/sim/phylo-hmm/optimized.posterior.decoding.probabilities.' + tag + '.' + rep
    output_viterbi_optimized_fn = '../../results/sim/phylo-hmm/optimized.viterbi.sequence.'  + tag + '.' + rep
    output_model_likelihoods_fn = '../../results/sim/phylo-hmm/optimized.model.likelihoods.'  + tag + '.' + rep
    output_optimized_params_fn = '../../results/sim/phylo-hmm/optimized.model.parameters.' + tag + '.' + rep
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


