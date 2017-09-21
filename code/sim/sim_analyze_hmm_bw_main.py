from sim_analyze_hmm_bw import *
from concordance_functions import *
sys.path.insert(0, '..')
import global_params as gp


# two options: just summary statistics for simulations, or also
# predictions
predict = False

#####
# get simulation arguments/options
#####

# have binary sequences been converted to nucleotides with seq-gen
# program?
seq_gen = False

# read in simulation arguments
tag, topology, species_to, species_from1, species_from2, \
    num_samples_species_to, num_samples_species_from1, num_samples_species_from2, \
    N0_species_to, N0_species_from1, N0_species_from2, \
    migration_from1, migration_from2, \
    expected_length_introgressed, \
    expected_num_introgressed_tracts, \
    has_ref_from1, has_ref_from2, \
    rho, outcross_rate, theta, num_sites, num_reps = \
    process_args(sys.argv)

num_samples = num_samples_species_to + \
    num_samples_species_from1 + num_samples_species_from2

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
results_filename = gp_dir + gp.sim_out_dir + '/analyze/' + gp.sim_out_prefix + tag + '_summary.txt'
hmm_filename = gp_dir + gp.sim_out_dir + '/analyze/' + 'hmm_parameters_' + tag + '.txt'

avg_results_filename = gp_dir + gp.sim_out_dir + '/analyze/' + gp.sim_out_prefix + 'avg_' + tag + '_summary.txt'

# write results headers
# - for training on single sim
f_out = open(results_filename, 'w')
output_dic = make_output_dic(states, species_to)
write_output_line(f_out, output_dic, True)
# - and using averaged parameters
avg_f_out = open(avg_results_filename, 'w')
avg_output_dic = make_output_dic(states, species_to)
write_output_line(avg_f_out, avg_output_dic, True)

# results files for tracts predicted to be and actually introgressed
# - for training on single sim
f_tracts_predicted = open(gp_dir + gp.sim_out_dir + '/analyze/' + \
                              gp.sim_out_prefix + tag + '_introgressed_tracts_predicted.txt', 'w')
f_tracts_actual = open(gp_dir + gp.sim_out_dir + '/analyze/' + \
                           gp.sim_out_prefix + tag + '_introgressed_tracts_actual.txt', 'w')
# - and using averaged parameters
avg_f_tracts_predicted = open(gp_dir + gp.sim_out_dir + '/analyze/' + \
                                  gp.sim_out_prefix + 'avg_' +  tag + '_introgressed_tracts_predicted.txt', 'w')
avg_f_tracts_actual = open(gp_dir + gp.sim_out_dir + '/analyze/' + \
                               gp.sim_out_prefix + 'avg_' + tag + '_introgressed_tracts_actual.txt', 'w')

#####
# theory (only need to do these calculations once)
#####

prob_topological_concordance = -1
prob_monophyletic_concordance = -1

"""
if theory:
    if theory_done:
        print 'reading theory results from file'
        f = open('out/' + results_filename, 'r')
        line = f.readline().split()
        assert(line[12] == 'prob_topological_concordance')
        assert(line[13] == 'prob_monophyletic_concordance')
        line = f.readline().split()
        prob_topological_concordance = float(line[12])
        prob_monophyletic_concordance =0 float(line[13])
        f.close()
"""

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
for i in range(num_reps):
    print i

    sim = read_one_sim(ms_f, num_sites, num_samples)
    assert sim != None, str(num_reps) + ' reps is not correct'
    
    if seq_gen:
        seq_fn = gp_dir + gp.sim_out_dir + '/seq-gen/' + gp.sim_out_prefix + 'seq_gen_' + tag + '_rep' + str(i) + '.fasta'
        sim[4] = read_fasta(seq_fn)
    else:
        # fill in the nonpolymorphic sites
        seqs_filled = fill_seqs(sim[4], sim[3], num_sites, fill_symbol)
        sim[4] = seqs_filled

    # fill in nonpolymorphic sites in sequences, and code for HMM
    # predictions
    seqs = sim[4]
    positions = sim[3]
    ref_seqs = [seqs[x] for x in ref_inds]
    # convert from binary to symbols indicating which reference
    # sequences each base matches
    seqs_coded = code_seqs(seqs, num_sites, ref_seqs, \
                               gp.match_symbol, gp.mismatch_symbol, \
                               gp.unknown_symbol, gp.unsequenced_symbol)

    init, emis, trans = initial_hmm_parameters(seqs_coded, species_to, \
                                                   index_to_species, states, \
                                                   unknown_species, \
                                                   gp.match_symbol, \
                                                   gp.mismatch_symbol, \
                                                   gp.unknown_symbol, \
                                                   expected_length_introgressed, \
                                                   expected_num_introgressed_tracts)

    # get parameters from training (but provide initial values)
    init_new, emis_new, trans_new = analyze_one(sim, seqs_coded, \
                                                    init, emis, trans, \
                                                    num_sites, fill_symbol, \
                                                    ref_inds, ref_seqs, \
                                                    num_samples_species_to, \
                                                    species_to, index_to_species, \
                                                    topology, states, unknown_species, \
                                                    len(states) - 1, \
                                                    output_dic, f_out, \
                                                    f_tracts_predicted, \
                                                    f_tracts_actual, i, \
                                                    train=True)
    init_all.append(init_new)
    emis_all.append(emis_new)
    trans_all.append(trans_new)

f_out.close()

ms_f.close()
ms_f = open(gp_dir + gp.sim_out_dir + '/ms/' +  gp.sim_out_prefix + tag + '.txt', 'r')

# get average of HMM parameters across all simulations and write them
# to file
init, emis, trans = average_hmm_params(init_all, emis_all, trans_all)
write_hmm_params(init, emis, trans, states, unknown_species, hmm_filename)

# loop through all reps again, this time using the averaged HMM
# parameters (no training)
for i in range(num_reps):

    print i

    sim = read_one_sim(ms_f, num_sites, num_samples)
    assert sim != None, str(num_reps) + ' reps is not correct'

    if seq_gen:
        seq_fn = gp_dir + gp.sim_out_dir + '/seq-gen/' + gp.sim_out_prefix + 'seq_gen_' + tag + '_rep' + str(i) + '.fasta'
        sim[4] = read_fasta(seq_fn)
    else:
        # fill in the nonpolymorphic sites
        seqs_filled = fill_seqs(seqs, positions, num_sites, fill_symbol)
        sim[4] = seqs_filled

    # fill in nonpolymorphic sites in sequences, and code for HMM
    # predictions
    seqs = sim[4]
    positions = sim[3]
    ref_seqs = [seqs[x] for x in ref_inds]
    # convert from binary to symbols indicating which reference
    # sequences each base matches
    seqs_coded = code_seqs(seqs, num_sites, ref_seqs, \
                               gp.match_symbol, gp.mismatch_symbol, \
                               gp.unknown_symbol, gp.unsequenced_symbol)


    # give averaged parameters this time
    init_new, emis_new, trans_new = analyze_one(sim, seqs_coded, \
                                                    init, emis, trans, \
                                                    num_sites, fill_symbol, ref_inds, \
                                                    ref_seqs, \
                                                    num_samples_species_to, \
                                                    species_to, index_to_species, \
                                                    topology, states, \
                                                    unknown_species, \
                                                    len(states) - 1, \
                                                    avg_output_dic, avg_f_out, \
                                                    avg_f_tracts_predicted, \
                                                    avg_f_tracts_actual, i, \
                                                    train=False)

avg_f_out.close()

#####
# loop through all other parameter sets, and all reps for each set
#####

"""
for params in param_sets:

    init, emis, trans = params

    f = open(outfilename, 'r')
    line = f.readline()

    for i in range(num_reps):
        print i
    
        analyze_one_rep(params, f)
"""


ms_f.close()
f_tracts_predicted.close()
f_tracts_actual.close()

