# TODO eliminate some of the copy pasta here from the code for
# analyzing my method

from sim_analyze_phylo import * 
sys.path.insert(0, '../sim')
from sim_analyze_hmm_bw import *
from concordance_functions import *
sys.path.insert(0, '..')
import global_params as gp

seq_gen = False

tag, topology, species_to, species_from1, species_from2, \
    num_samples_species_to, num_samples_species_from1, num_samples_species_from2, \
    N0_species_to, N0_species_from1, N0_species_from2, \
    migration_from1, migration_from2, \
    expected_length_introgressed, \
    expected_num_introgressed_tracts, \
    has_ref_from1, has_ref_from2, \
    rho, outcross_rate, theta, num_sites, num_reps = \
    process_args(sys.argv)

num_samples = num_samples_species_to + num_samples_species_from1 + num_samples_species_from2

# species_to always comes first
index_to_species = [species_to] * num_samples_species_to + \
    [species_from1] * num_samples_species_from1 + \
    [species_from2] * num_samples_species_from2

#####
# sequence and HMM symbols (mostly moved to general_params)
#####

fill_symbol = '0'

#####
# reference sequences for each species and states
#####


# take first index from each population to be reference sequence
ref_ind_species_to = 0
ref_ind_species_from1 = num_samples_species_to
ref_ind_species_from2 = num_samples_species_to + num_samples_species_from1
ref_inds = [ref_ind_species_to]
states = [species_to, species_from1]
unknown_species = None
if has_ref_from1:
    ref_inds.append(ref_ind_species_from1)
else:
    unknown_species = species_from1
if species_from2 != None:
    states.append(species_from2)
    if has_ref_from2:
        ref_inds.append(ref_ind_species_from2)
    else:
        unknown_species = species_from2

# if there are three species and the second is the species that has no
# reference, flip the order of the states so that the species with no
# reference always comes last; this will ensure that the indices of
# the species in states correspond to the indices of the references
# (and the sequence codings later); ACTUALLY just force the unknown
# species to come last
if species_from2 != None:
    assert has_ref_from1
#if species_from2 != None and not has_ref_from1:
#    states = states[0] + states[2] + states[1]

#####
# keep track of HMM parameters to average at end
#####

init_all = []
emis_all = []
trans_all = []

#####
# output files
#####

gp_dir = '../'
# for reading output from ms
ms_f = open(gp_dir + gp.sim_out_dir + '/ms/' +  gp.sim_out_prefix + tag + '.txt', 'r')

# for writing results of analysis
results_filename = gp_dir + gp.sim_out_dir + '/analyze/' + gp.sim_out_prefix + tag + '_summary_phylo.txt'
hmm_filename = gp_dir + gp.sim_out_dir + '/analyze/' + 'hmm_parameters_' + tag + '_phylo.txt'
avg_results_filename = gp_dir + gp.sim_out_dir + '/analyze/' + gp.sim_out_prefix + 'avg_' + tag + 'summary_phylo.txt'

# write results headers
# - for training on single sim
f_out = open(results_filename, 'w')
output_dic = make_output_dic(states, species_to)
write_output_line(f_out, output_dic, True)
# - and using averaged parameters
avg_f_out = open(results_filename, 'w')
avg_output_dic = make_output_dic(states, species_to)
write_output_line(avg_f_out, avg_output_dic, True)

# results files for tracts predicted to be and actually introgressed
# - training on single sim
f_tracts_predicted = open(gp_dir + gp.sim_out_dir + '/analyze/' + \
                              gp.sim_out_prefix + tag + '_introgressed_tracts_predicted_phylo.txt', 'w')
# - averaged params
avg_f_tracts_predicted = open(gp_dir + gp.sim_out_dir + '/analyze/' + \
                                  gp.sim_out_prefix + 'avg_' + tag + '_introgressed_tracts_predicted_phylo.txt', 'w')

#####
# for actual parameters, analyze simulation results and get hmm
# parameters
#####

# do this for individual simulations and for all together (?); that
# is, loop through all reps, and get results just from that one rep;
# then average all the HMM parameters and do analysis for the same
# reps again

# loop through all reps; for each, read in simulation (and store it),
# and predict introgressed tracts by training on those sequences
# using phylo-hmm
for i in range(num_reps):
    print i

    # trees, recomb sites, segsites, positions, seqs
    sim = read_one_sim(ms_f, num_sites, num_samples)
    assert sim != None, str(num_reps) + ' reps is not correct'

    seq_fn = ''
    if seq_gen:
        # simulated sequences already generated with ms and then seq-gen
        # DON'T USE THIS IN CURRENT FORM
        seq_fn = gp_dir + gp.sim_out_dir + '/seq-gen/' + gp.sim_out_prefix + 'seq_gen_' + tag + '_rep' + str(i) + '.fasta'
        sim[4] = read_fasta(seq_fn)
    else:
        # fill in the nonpolymorphic sites
        seqs_filled = fill_seqs(sim[4], sim[3], num_sites, fill_symbol)
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

    init_all.append(init_new)
    emis_all.append(emis_new)
    trans_all.append(trans_new)

f_out.close()
f_tracts_predicted.close()

"""
# get average of HMM parameters across all simulations and write them
# to file
init, emis, trans = average_hmm_params(init_all, emis_all, trans_all)
write_hmm_params(init, emis, trans, states, unknown_species, hmm_filename)

# loop through all reps again, this time using the averaged HMM
# parameters (no training)
for i in range(num_reps):
    # give averaged parameters this time
    analyze_one(sim, init, emis, trans, \
                    num_sites, fill_symbol, ref_inds, ref_seqs, \
                    num_samples_species_to, species_to, index_to_species, \
                    topology, states, unknown_species, \
                    avg_output_dic, avg_f_out, train=False)
"""
avg_f_out.close()
avg_f_tracts_predicted.close()


#####
# loop through all other parameter sets, and all reps for each set
#####

# laytah
"""
for params in param_sets:

    init, emis, trans = params

    f = open(outfilename, 'r')
    line = f.readline()

    for i in range(num_reps):
        print i
    
        analyze_one_rep(params, f)
"""


