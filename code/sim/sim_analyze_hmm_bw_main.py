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
# sequence and HMM symbols
#####

fill_symbol = '0'
unsequenced_symbol = 'N'
match_symbol = '+'
mismatch_symbol = '-'
unknown_symbol = '?'

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
# (and the sequence codings later)
if species_from2 != None and not has_ref_from1:
    states = states[0] + states[2] + states[1]

#####
# output files
#####

gp_dir = '../'
outfilename = gp_dir + gp.sim_out_dir +  gp.sim_out_prefix + tag + '.txt'
results_filename = gp_dir + gp.sim_out_dir + gp.sim_out_prefix + tag + '_summary.txt'

# write results headers
fout = open(results_filename, 'w')
output_dic = make_output_dic(states, species_to)
write_output_line(fout, output_dic, True)

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
# loop through all reps
#####

f = open(outfilename, 'r')
line = f.readline()
n = 0
while line != '' and n < num_reps:

    if line == '//\n':
        print n

        #####
        # read in results from next simulation
        #####

        # positions converted to integers within this function
        trees, recomb_sites, segsites, positions, seqs = read_sim(f, num_sites, num_samples)
        assert len(positions) == segsites

        #####
        # fill in nonpolymorphic sites in sequences, and code for HMM
        # predictions
        #####

        # fill in the nonpolymorphic sites
        seqs_filled = fill_seqs(seqs, positions, num_sites, fill_symbol)
        ref_seqs = [seqs_filled[x] for x in ref_inds]
        # convert from binary to symbols indicating which reference
        # sequences each base matches
        print positions
        seqs_coded = code_seqs(seqs_filled, num_sites, ref_seqs, \
                                   match_symbol, mismatch_symbol, \
                                   unknown_symbol, unsequenced_symbol)

        ########
        # figure out which sites are actually introgressed by
        # separately looking at the tree for each stretch without
        # recombination
        ########

        # sequence of states, one for each site and strain
        actual_state_seq = [[] for i in range(num_samples_species_to)]
        # keep track of whether each tree is concordant with the species tree
        concordant = []
        # and how many lineages of the to species are left when it
        # first joins another species
        num_lineages_at_join = []
        # and how many bases are introgressed in total in each strain
        num_introgressed = [0] * num_samples_species_to
        # loop through the trees for all blocks with no recombination
        # within them
        for ti in range(len(trees)):

            # note that species indices/labels are shifted to start at
            # 0 instead of 1
            t = trees[ti]
            
            # is this tree concordant with the species tree? (only
            # checks whether the to species is monophyletic, which
            # indicates that ILS not possible)
            if is_concordant(t, index_to_species, species_to):
                concordant.append(True)
            else:
                concordant.append(False)

            # identify sequences that are introgressed from the one or
            # two other species, based on coalescent tree; could clean
            # this up a little
            introgressed = None
            num_lineages_at_join_current = None
            # three species
            if species_from2 != None:
                introgressed, num_lineages_at_join_current = \
                    find_introgressed_3(t, species_to, topology, index_to_species)
            # two species
            else:
                # introgressed is a list of species (one entry for
                # each individual in to species)
                introgressed, num_lineages_at_join_current = \
                    find_introgressed_2(t, topology[2], species_to, index_to_species)
            print introgressed[0]
            # number of lineages that were present when all
            # populations joined
            num_lineages_at_join.append(num_lineages_at_join_current)

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

        ########
        # predict whether each site in each cer strain is
        # introgressed with a hidden markov model
        ########

        print 'HMM'

        # predicted is a list with one entry for each individual of
        # species_to; each entry gives predicted species for each
        # position
        predicted, hmm = predict_introgressed_hmm(seqs_coded, species_to, \
                                                      index_to_species, \
                                                      states, \
                                                      unknown_species, \
                                                      match_symbol, mismatch_symbol, \
                                                      unknown_symbol, \
                                                      expected_length_introgressed, \
                                                      expected_num_introgressed_tracts)

        num_correct, num_introgressed_correct, actual_lens, predicted_lens, \
            num_predicted_tracts_actual, num_actual_tracts_predicted, \
            num_introgressed_tracts, num_not_introgressed_tracts, \
            num_predicted_introgressed_tracts, num_predicted_not_introgressed_tracts = \
            evaluate_predicted(predicted, actual_state_seq, species_to)
       

        groups = group_actual_predicted_bases(actual_state_seq, predicted, states)
        for group in groups:
            output_dic = update_value(output_dic, 'num_bases_actual_' + group[0] + \
                                          '_predicted_' + group[1], \
                                          groups[group])

        ### tracts
        blocks_predicted, blocks_actual = \
            evaluate_predicted_blocks(predicted, actual_state_seq, \
                                          species_to, states)

        d_actual_predicted, d_predicted_actual, d_actual_counts, d_predicted_counts = \
            group_actual_predicted_blocks(blocks_actual, blocks_predicted, states, num_samples_species_to)
        for group in d_actual_predicted:
            output_dic = update_value(output_dic, 'num_tracts_actual_' + group[0] + \
                                          '_predicted_' + group[1], \
                                          d_actual_predicted[group])
        for group in d_predicted_actual:
            output_dic = update_value(output_dic, 'num_tracts_predicted_' + group[0] + \
                                          '_actual_' + group[1], \
                                          d_predicted_actual[group])
        for state in states:
            output_dic = update_value(output_dic, 'num_tracts_actual_' + state, \
                                          d_actual_counts[state])
            output_dic = update_value(output_dic, 'num_tracts_predicted_' + state, \
                                          d_predicted_counts[state])

        # tract lengths, actual and predicted; list, not average
        for state in states:
            tract_lengths = [b[2] for b in filter(lambda x: x[0] == state, blocks_actual)]
            output_dic = update_value(output_dic, 'tract_lengths_actual_' + state, tract_lengths)

            tract_lengths = [b[2] for b in filter(lambda x: x[0] == state, blocks_predicted)]
            output_dic = update_value(output_dic, 'tract_lengths_predicted_' + state, tract_lengths)


        # HMM parameters
        for i in range(len(states)):
            output_dic = update_value(output_dic, 'init_' + states[i], hmm.init[i])
            
        for i in range(len(states)):
            for j in range(len(states)):
                output_dic = update_value(output_dic, \
                                              'trans_' + states[i] + '_' + states[j], \
                                              hmm.trans[i][j])

        # so here emis_cer_par is going to be the probability that cer
        # state emits symbol that matches par, so *+; note that these
        # don't have to add to 1 (states includes unknown if
        # appropriate)
        for i in range(len(states)):
            to_states = states
            for j in range(len(states)):
                total = 0
                for s in hmm.emis[i].keys():
                    # for unknown state, include --- (as well as '?')
                    if states[j] == unknown_species:
                        if match_symbol not in s:
                            total += hmm.emis[i][s]
                    elif s[j] == match_symbol:
                        total += hmm.emis[i][s]
                output_dic = update_value(output_dic, 'emis_'  + states[i] + '_' + \
                                              to_states[j], total)

        ########
        # predict whether each site in each cer strain is
        # introgressed by windowing (each window is in each strain
        # is either completely introgressed or not introgressed)
        ########

        # meh

        ########
        # sequence identities
        ########

        s = seq_id(seqs_filled, index_to_species, states)
        for i in range(len(states)):
            for j in range(i, len(states)):
                output_dic = update_value(output_dic, \
                                              'avg_identity_' + states[i] + '_' + states[j], \
                                              s[i][j])

        #####
        # write results to file
        #####

        write_output_line(fout, output_dic, False)

        fout.flush()

        sys.stdout.flush()

        for b in blocks_predicted:
            f_tracts_predicted.write(str(n) + ' ' +  ' '.join([str(x) for x in b]) + '\n')
        for b in blocks_actual:
            f_tracts_actual.write(str(n) + ' ' + ' '.join([str(x) for x in b]) + '\n')

        f_tracts_predicted.flush()
        f_tracts_actual.flush()

        print 'DONE'

        n += 1

    line = f.readline()


f.close()
fout.close()
f_tracts_predicted.close()
f_tracts_actual.close()

