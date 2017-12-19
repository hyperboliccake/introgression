import sys
import analyze
sys.path.append('..')
import global_params as gp
sys.path.append('../sim')
import sim_predict

##======
# read in analysis parameters
##======

refs, strains, args = predict.process_args(sys.argv)

# refs = {'cer':('S288c', '../../data/', 'S288C-SGD_R64'), ...]
# strains = {'cer':[('strain1', '../../data/'), ...], ...}

##======
# loop through all sequences and predict introgression
##======

gp_dir = '../'

out_f = open(gp.analysis_out_dir_absolute, 'w')
introgression_f
probs_f

for chrm in gp.chrms:

    for strain in strains:

        fn = fp_dir + gp.alignment_dir + '_'.join(args['species']) + '_' + strain + \
             '_mafft' + gp.alignment_suffix
        ref_seqs, predict_seq = read_aligned_seqs(fn, refs, strain)

        ##======
        # predict introgressed/non-introgressed tracts
        ##======

        state_seq, probs, hmm, hmm_init, ps = \
            predict.predict_introgressed(predict_seq, ref_seqs, args,\
                                         train = True, method='posterior')

        state_seq_blocks = sim_process.convert_to_blocks(state_seq, \
                                                         args['species'])
        ##======
        # output
        ##======
        
        # 
        predict.write_positions(ps)

        predict.write_blocks(state_seq_blocks, ps, , strain, chrm)
        

        # summary info about HMM (before training)
        sim_predict.write_hmm_line(hmm_init, out_init_f, i==0) 
        
        # summary info about HMM (after training)
        sim_predict.write_hmm_line(hmm, out_f, i==0) 

        # locations of introgression
        sim_process.write_introgression_blocks(state_seq_blocks, introgression_f, \
                                               strain, args['species'])

        # probabilities at each site
        sim_process.write_state_probs(probs, prob_f, strain)

out_f.close()
introgression_f.close()
prob_f.close()


"""


tag, topology, species_to, species_from1, species_from2, \
    num_samples_species_to, num_samples_species_from1, num_samples_species_from2, \
    N0_species_to, N0_species_from1, N0_species_from2, \
    migration_from1, migration_from2, \
    expected_tract_lengths, \
    expected_num_tracts, \
    has_ref_from1, has_ref_from2, \
    rho, outcross_rate, theta, num_sites, num_reps = \
    sim.process_args(sys.argv)

# reference names for actual species
states = []
i = -1
if species_from2 != None:
    states = [sys.argv[i]]
    i -= 1
states = [sys.argv[i]] + states
i -= 1
states = [sys.argv[i]] + states

# reference names for sim species (in same order)
sim_states = [species_to, species_from1]
if species_from2 != None:
    sim_states.append(species_from2)

# is one of the states actually unknown (and thus not a reference)?
unknown_state = False
if species_from2 != None:
    if not has_ref_from2:
        unknown_state = True
elif not has_ref_from1:
    unknown_state = True

refs = copy.deepcopy(states)
if unknown_state:
    refs = refs[:-1]
master_ref = refs[0]

gp_dir = '../'

# get the filenames of all the strains we're going to predict; this is
# probably not the simplest way to set this up but whatever
strain_names = []
for d in gp.non_ref_dirs[states[0]]:
    fns = os.listdir(d)
    fns = filter(lambda x: x[-3:] == '.fa', fns)
    strain_names += [x[:x.find('_')] for x in fns]
strain_names = list(set(strain_names))
f_strain_list = open(gp.analysis_out_dir_absolute + '/' + tag + '/strain_list.txt', 'w')
f_strain_list.write('\n'.join(strain_names) + '\n')
f_strain_list.close()

fn_hmm = gp_dir + gp.sim_out_dir + 'hmm_parameters_' + tag + '.txt'
init, emis, trans = read_hmm_params(fn_hmm, states, sim_states, unknown_state)

fn_out = gp.analysis_out_dir_absolute + 'introgressed_hmm_' + tag + '.txt'
f_out = open(fn_out, 'w')
f_out.write('strain\tchromosome\talignment_block_label\tstrand\tpredicted_reference\tregion_start\tregion_end\tnumber_non_gap_sites\n')


# one at a time because memory
for x in strain_names:
    print x    
    seqs, seq_inds, psx = get_seqs(x, refs, gp.chrms, \
                                       gp.match_symbol, gp.mismatch_symbol, \
                                       gp.unknown_symbol, gp.unsequenced_symbol, \
                                       master_ref)
    print len(seqs), 'seqs'
    sys.stdout.flush()
    print states
    print init
    print emis
    print trans
    predicted, hmm = predict_introgressed_hmm(seqs, states, \
                                                  sim_states, unknown_state, \
                                                  init, emis, trans)

    blocks = get_predicted_tracts(predicted, psx, master_ref, seq_inds)
    write_predicted_tracts(blocks, f_out)
    
    sys.stdout.flush()

f_out.close()
"""
