from sim_analyze_hmm_bw import *
from concordance_functions import *
sys.path.insert(0, '..')
import global_params as gp

tag, model, N0, include_bay, include_unk,\
    num_samples_cer, num_samples_par, num_samples_bay,\
    par_cer_migration, bay_cer_migration,\
    t_par_cer, t_bay_par_cer,\
    num_sites, rho, theta, outcross_rate, num_reps = \
    process_args(sys.argv)

num_samples = num_samples_par + num_samples_cer + num_samples_bay
mu = 1.84 * 10 ** -10
theta = mu * 2 * num_sites * N0

index_to_species = [gp.Species.cer] * num_samples_cer + \
    [gp.Species.par] * num_samples_par + \
    [gp.Species.bay] * num_samples_bay

# take first index from each population to be reference sequence
ref_ind_cer = 0
ref_ind_par = num_samples_cer
ref_ind_bay = None
if include_bay:
    ref_ind_bay = num_samples_cer + num_samples_par

species_to_predict = gp.Species.cer


#####
# stats to keep track of
#####

# for each sim, number of introgressed bases in each strain
num_introgressed_cer = [[0 for i in range(num_samples_cer)] for n in range(num_reps)]

# success of HMM predictions
fraction_correct_all = [[] for i in range(num_reps)]

fraction_introgressed_correct_all = [[] for i in range(num_reps)]

predicted_lens_all = [[] for i in range(num_reps)]

# actual lengths don't depend on prediction strategy
actual_lens_all = [[] for i in range(num_reps)]

fpr_cer = [0] * num_reps
fpr_cer_only_par = [0] * num_reps

# id between all cer and cer ref
avg_id_cer = [0] * num_reps

# id between all cer and par ref 
avg_id_cer_par = [0] * num_reps

num_lineages_at_join = [[] for n in range(num_reps)]

concordant = [[] for n in range(num_reps)]

prob_topological_concordance = 0

# symbol codings
# 0    1       2       3          4    5       6    7
# cer, cerpar, cerbay, cerparbay, par, parbay, bay, none
#hmm_symbol = {'cer':0, 'cerpar':1, 'cerbay':2, 'cerparbay':3, 'par':4, 'parbay':5, 'bay':6, 'none':7}

#hmm_symbols = []
#individual_symbols = ['+', '-', '?']
#for c in individual_symbols:
#    for p in individual_symbols:
#        for b in individual_symbols:
#            hmm_symbols.append(c + p + b)

fill_symbol = '0'
unsequenced_symbol = 'N'
match_symbol = '+'
mismatch_symbol = '-'
unknown_symbol = '?'


introgressed_states = [gp.Species.cer, gp.Species.par]
ref_inds = [ref_ind_cer, ref_ind_par]
if include_bay:
    introgressed_states.append(gp.Species.bay)
    ref_inds.append(ref_ind_bay)
elif include_unk:
    introgressed_states.append(gp.Species.unk)

#####
# loop through all reps
#####

gp_dir = '../'
outfilename = gp.sim_out_prefix + tag + '.txt'
results_filename = gp.sim_out_prefix + tag + '_summary.txt'

# write results headers
fout = open(gp_dir + gp.sim_out_dir + results_filename, 'w')
output_dic = make_output_dic(introgressed_states, species_to_predict)
write_output_line(fout, output_dic, True)

# results files for tracts predicted to be and actually introgressed
f_tracts_predicted = open(g;outdir + 'sim_out_' + tag + '_introgressed_tracts_predicted.txt', 'w')
f_tracts_actual = open(outdir + 'sim_out_' + tag + '_introgressed_tracts_actual.txt', 'w')

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
        prob_monophyletic_concordance = float(line[13])
        f.close()
"""

f = open(outdir + outfilename, 'r')
line = f.readline()
n = 0
while line != '' and n < num_reps:

    if line == '//\n':
        print n

        #####
        # read in results from next simulation
        #####

        trees, recomb_sites, segsites, positions, seqs = read_sim(f, num_sites, num_samples)
        assert len(positions) == segsites

        #####
        # fill in nonpolymorphic sites in sequences, and code for HMM
        # predictions
        #####

        seqs_filled = fill_seqs(seqs, positions, num_sites, fill_symbol)
        ref_seqs = [seqs_filled[x] for x in ref_inds]
        seqs_coded = code_seqs(seqs_filled, positions, num_sites, ref_seqs, \
                                   match_symbol, mismatch_symbol, unknown_symbol, unsequenced_symbol)
        
        ########
        # figure out which sites are actually introgressed by
        # separately looking at tree for each without recombination
        ########

        # sequence of states, one for each site and strain
        actual_state_seq = [[] for i in range(num_samples)]

        # stuff to keep track of
        concordant = []
        num_lineages_at_join = []
        num_introgressed = [0] * num_samples
        # loop through the trees for all blocks with no recombination
        # within them
        for ti in range(len(trees)):

            # note that species indices/labels are shifted to start at
            # 0 instead of 1
            t = trees[ti]
            
            # keep track of how often gene tree is concordant with
            # species tree, by whatever definition we're interested in
            if is_concordant(t, index_to_species, species_to_predict):
                concordant.append(True)
            else:
                concordant.append(False)
            # identify cerevisiae sequences that are introgressed from
            # paradoxus or bayanus, based on coalescent tree; return
            # result for all individuals (not just cerevisiae) to
            # simplify indexing and allow flexibility later on
            introgressed = None
            if include_bay:
                # TODO finish this function 
                introgressed = find_introgressed_3(t, t_par_cer, t_bay_par_cer, \
                                                       gp.Species.cer, gp.Species.par, \
                                                       gp.Species.bay, index_to_species)
            else:
                # introgressed is a list of species
                t_state_seqs, t_num_lineages = find_introgressed_2(t, t_par_cer, \
                                                                       gp.Species.cer, gp.Species.par, \
                                                                       index_to_species)

            # number of lineages that were present when all
            # populations joined
            num_lineages_at_join.append(t_num_lineages)

            # number of sites in the current block of sequence
            num_sites_t = recomb_sites[ti]

            # for all strains that have this block introgressed, add
            # the length of the block to the total number of
            # introgressed sites across all strains; also update the
            # state sequence
            for i in range(num_samples):
                if t_state_seqs[i] != index_to_species[i]:
                    num_introgressed[i] += num_sites_t
                actual_state_seq[i] += [t_state_seqs[i]] * num_sites_t

        num_introgressed_species_to_predict = []
        for i in range(num_samples):
            if index_to_species[i] == species_to_predict:
                num_introgressed_species_to_predict.append(num_introgressed[i])
        output_dic = update_value(output_dic, 'num_introgressed_bases_actual_' + species_to_predict, 
                                  num_introgressed_species_to_predict)
                
        ########
        # predict whether each site in each cer strain is
        # introgressed with a hidden markov model
        ########

        print 'HMM'

        # predicted is a list; for each index that wasn't the species
        # we wanted to predict introgression, the entry is None; for
        # other indices, the entry is a list of predicted species at
        # each position

        predicted, hmm = predict_introgressed_hmm(seqs_coded, species_to_predict, \
                                                      index_to_species, \
                                                      introgressed_states, \
                                                      match_symbol, mismatch_symbol, \
                                                      unknown_symbol, \
                                                      gp.expected_length_introgressed, \
                                                      gp.expected_num_introgressed_tracts)

        num_predicted_introgressed = []
        for i in range(num_samples):
            if index_to_species[i] == species_to_predict:
                x = 0
                for j in range(len(predicted[i])):
                    # this adds one if introgressed (true is 1, false is 0)
                    x += (predicted[i][j] != species_to_predict)
                num_predicted_introgressed.append(x)
        output_dic = update_value(output_dic, 'num_introgressed_bases_predicted_' + species_to_predict, num_predicted_introgressed)

        num_correct, num_introgressed_correct, actual_lens, predicted_lens, \
            num_predicted_tracts_actual, num_actual_tracts_predicted, \
            num_introgressed_tracts, num_not_introgressed_tracts, \
            num_predicted_introgressed_tracts, num_predicted_not_introgressed_tracts = \
            evaluate_predicted(predicted, actual_state_seq, \
                                   index_to_species, species_to_predict)
        output_dic = update_value(output_dic, 'num_bases_correct_' + species_to_predict, num_correct)

        # true positives
        output_dic = update_value(output_dic, 'num_actual_introgressed_bases_predicted_' + species_to_predict, \
                                      num_introgressed_correct)


        ### tracts
        blocks_predicted, blocks_actual = \
            evaluate_predicted_blocks(predicted, actual_state_seq, \
                                          index_to_species, species_to_predict)
        
        # positives
        output_dic = update_value(output_dic, 'num_introgressed_tracts_actual_' + species_to_predict, \
                                      count_blocks(blocks_actual, species_to_predict, ['positive']))

        # negatives
        output_dic = update_value(output_dic, 'num_not_introgressed_tracts_actual_' + species_to_predict, \
                                      count_blocks(blocks_actual, species_to_predict, ['negative']))

        # predicted positives
        output_dic = update_value(output_dic, 'num_introgressed_tracts_predicted_' + species_to_predict, \
                                      count_blocks(blocks_predicted, species_to_predict, ['positive']))

        # predicted negatives
        output_dic = update_value(output_dic, 'num_not_introgressed_tracts_predicted_' + species_to_predict, \
                                      count_blocks(blocks_predicted, species_to_predict, ['negative']))

        # true positives (number of actual that tracts overlap with
        # predicted introgressed block)
        output_dic = update_value(output_dic, 'num_actual_introgressed_tracts_predicted_' + species_to_predict, \
                                      count_blocks(blocks_actual, species_to_predict, ['true', 'positive']))

        # false positives (number of predicted tracts that don't
        # overlap with actually introgressed block)
        output_dic = update_value(output_dic, 'num_predicted_introgressed_tracts_not_actual_' + species_to_predict, \
                                      count_blocks(blocks_predicted, species_to_predict, ['false', 'positive']))

        # tract lengths, actual and predicted
        output_dic = update_value(output_dic, 'introgressed_tract_lengths_actual_' + species_to_predict, \
                                      [b[2] for b in filter(lambda x: x != species_to_predict, blocks_actual)])
        output_dic = update_value(output_dic, 'introgressed_tract_lengths_predicted_' + species_to_predict, \
                                      [b[2] for b in filter(lambda x: x != species_to_predict, blocks_predicted)])

        # HMM parameters
        for i in range(len(introgressed_states)):
            output_dic = update_value(output_dic, 'init_' + introgressed_states[i], hmm.emis[i])
            
        for i in range(len(introgressed_states)):
            for j in range(len(introgressed_states)):
                output_dic = update_value(output_dic, 'trans_' + introgressed_states[i] + '_' + introgressed_states[j], \
                                              hmm.trans[i][j])

        # so here emis_cer_par is going to be the probability that cer
        # state emits symbol that matches par, so *+; note that these
        # don't have to add to 1 (introgressed_states includes unknown
        # if appropriate)
        for i in range(len(introgressed_states)):
            to_states = introgressed_states
            # TODO deal with unknown state? not here really but not to
            # treat it like the other states in the hmm
            for j in range(len(introgressed_states)):
                total = 0
                for s in hmm.emis[i].keys():
                    if s[j] == match_symbol:
                        total += hmm.emis[i][s]
                output_dic = update_value(output_dic, 'emis_'  + introgressed_states[i] + '_' + \
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

        s = seq_id(seqs_coded, index_to_species, introgressed_states)
        for i in range(len(introgressed_states)):
            for j in range(i, len(introgressed_states)):
                output_dic = update_value(output_dic, \
                                              'avg_identity_' + introgressed_states[i] + '_' + introgressed_states[j], \
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


