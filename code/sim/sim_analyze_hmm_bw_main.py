from sim_analyze_hmm_bw import *
from concordance_functions import *
sys.path.insert(0, '..')
import global_params as gp

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
outfilename = gp_dir + gp.sim_out_dir +  gp.sim_out_prefix + tag + '.txt'
results_filename = gp_dir + gp.sim_out_dir + gp.sim_out_prefix + tag + '_summary.txt'
hmm_filename = gp_dir + gp.sim_out_dir + 'hmm_parameters_' + tag + '.txt'

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
f_tracts_predicted = open(gp_dir + gp.sim_out_dir + gp.sim_out_prefix + tag + '_introgressed_tracts_predicted.txt', 'w')
f_tracts_actual = open(gp_dir + gp.sim_out_dir + gp.sim_out_prefix + tag + '_introgressed_tracts_actual.txt', 'w')

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

    sim = read_one_sim(f, num_sites, num_samples)
    assert sim != None, str(num_reps) + ' reps is not correct'
    
    # initial values of hmm parameters
    seqs = sim[4]
    init, emis, trans = initial_hmm_parameters(seqs, species_to, \
                                                   index_to_species, states, \
                                                   unknown_species, \
                                                   match_symbol, mismatch_symbol, \
                                                   unknown_symbol, \
                                                   expected_length_introgressed, \
                                                   expected_num_introgressed_tracts)

    # get parameters from training (but provide initial values)
    analyze_one(sim, init, emis, trans, \
                    num_sites, fill_symbol, ref_inds, ref_seqs, \
                    num_samples_species_to, species_to, index_to_species, \
                    topology, states, unknown_species, \
                    output_dic, f_out, train=True)

f_out.close()

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
avg_f_out.close()

#####
# loop through all other parameter sets, and all reps for each set
#####

for params in param_sets:

    init, emis, trans = params

    f = open(outfilename, 'r')
    line = f.readline()

    for i in range(num_reps):
        print i
    
        analyze_one_rep(params, f)



f.close()
f_tracts_predicted.close()
f_tracts_actual.close()

