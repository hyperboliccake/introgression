import sys
import os
import predict
import gzip
sys.path.append('..')
import global_params as gp
sys.path.append('../sim')
import sim_process
sys.path.append('../align')
import align_helpers
sys.path.append('../misc')
import read_fasta

##======
# read in analysis parameters
##======

args = predict.process_predict_args(sys.argv[1:])

# refs: {'cer':('S288c', '../../data/', 'S288C-SGD_R64'), ...]
# strains: {'cer':[('strain1', '../../data/'), ...], ...}

strain_dirs = align_helpers.get_strains(align_helpers.flatten(gp.non_ref_dirs.values()))

##======
# output files and if and where to resume
##======

resume = False
open_mode = 'a'
if not resume:
    open_mode = 'w'

gp_dir = ''  # '../'

if not os.path.isdir(gp.analysis_out_dir_absolute + args['tag']):
    os.makedirs(gp.analysis_out_dir_absolute + args['tag'])

# positions
# TODO move this to more general location and make separate files for
# each strain x chrm
ps_fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + 'positions_' + \
        args['tag'] + '.txt.gz'

# introgressed blocks
blocks_f = {}
for s in args['states']:
    blocks_f[s] = open(gp.analysis_out_dir_absolute + args['tag'] + '/' + \
                       'blocks_' + s + '_' + args['tag'] + '.txt', \
                       open_mode)
    if not resume:
        predict.write_blocks_header(blocks_f[s])

# HMM parameters
hmm_init_fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + 'hmm_init_' + \
              args['tag'] + '.txt'
hmm_init_f = open(hmm_init_fn, open_mode)
hmm_fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + 'hmm_' + \
         args['tag'] + '.txt'
hmm_f = open(hmm_fn, open_mode)
emis_symbols = predict.get_emis_symbols(args['known_states'])
if not resume:
    predict.write_hmm_header(args['known_states'], args['unknown_states'], \
                             emis_symbols, hmm_init_f)
    predict.write_hmm_header(args['known_states'], args['unknown_states'], \
                             emis_symbols, hmm_f)

# posterior probabilities
probs_fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + 'probs_' + \
           args['tag'] + '.txt.gz'

# figure out which species x chromosomes we've completed already
# and also make sure file is gzipped properly
completed = {} # keyed by chrm, list of strains
if resume:
    probs_f = gzip.open(probs_fn, 'rb')
    probs_temp_fn = probs_fn[:-7] + '.temp.txt.gz'
    probs_temp_f = gzip.open(probs_temp_fn, 'wb')
    line = probs_f.readline()
    probs_temp_f.write(line)
    print 'reading already completed:'
    while line != '':
        x1 = line.find('\t')
        strain = line[:x1]
        x2 = line[x1+1:].find('\t') + x1 + 1
        chrm = line[x1+1:x2]
        if not completed.has_key(chrm):
            completed[chrm] = []
        completed[chrm].append(strain)
        print strain, chrm
        try:
            line = probs_f.readline()
            probs_temp_f.write(line)
        except Exception as e:
            print e
            break
    probs_f.close()
    probs_temp_f.close()
    os.system('mv ' + probs_temp_fn + ' ' + probs_fn)

    ps_f = gzip.open(ps_fn, 'rb')
    ps_temp_fn = ps_fn[:-7] + '.temp.txt.gz'
    ps_temp_f = gzip.open(ps_temp_fn, 'wb')
    line = ps_f.readline()
    ps_temp_f.write(line)
    while line != '':
        try:
            line = ps_f.readline()
            ps_temp_f.write(line)
        except Exception as e:
            print e
            break
    ps_f.close()
    ps_temp_f.close()
    os.system('mv ' + ps_temp_fn + ' ' + ps_fn)

    

write_ps = True
if write_ps:
    ps_f = gzip.open(ps_fn, open_mode + 'b')

probs_f = gzip.open(probs_fn, open_mode + 'b')

##======
# loop through all sequences and predict introgression
##======


for chrm in gp.chrms:

    for strain, strain_dir in strain_dirs:

        if resume and completed.has_key(chrm) and strain in completed[chrm]:
            print 'already finished:', strain, chrm
            continue

        print 'working on:', strain, chrm
        
        ref_prefix = '_'.join(gp.alignment_ref_order)
        fn = gp_dir + gp.alignments_dir + ref_prefix + '_' + strain + \
             '_chr' + chrm + '_mafft' + gp.alignment_suffix
        try:
            headers, seqs = read_fasta.read_fasta(fn)
        except Exception as e:
            print 'no alignment for', strain, chrm
            continue

        ref_seqs = seqs[:-1]
        predict_seq = seqs[-1]
        
        # TODO
        #ps_fn = gp.analysis_out_dir_absolute + '/positions/positions_' + \
        #        ref_prefix + '_' + strain + '_chr' + chrm + '.txt.gz'

        ##======
        # predict introgressed/non-introgressed tracts
        ##======

        state_seq, probs, hmm, hmm_init, ps = \
            predict.predict_introgressed(ref_seqs, predict_seq, args, \
                                         train = True)

        # hack
        state_seq_blocks = sim_process.convert_to_blocks({1:state_seq}, \
                                                         args['states'])[1]
        ##======
        # output
        ##======
        
        # the positions actually used in predictions (alignment
        # columns with no gaps)
        if write_ps:
            predict.write_positions(ps, ps_f, strain, chrm)

        # blocks predicted to be introgressed, separate files for each species
        for s in state_seq_blocks:
            predict.write_blocks(state_seq_blocks[s], ps, blocks_f[s], strain, chrm, s)        

        # summary info about HMM (before training)
        predict.write_hmm(hmm_init, hmm_init_f, strain, chrm, emis_symbols)
        
        # summary info about HMM (after training)
        predict.write_hmm(hmm, hmm_f, strain, chrm, emis_symbols) 

        # probabilities at each site
        predict.write_state_probs(probs, probs_f, strain, chrm)

for k in blocks_f:
    blocks_f[k].close()
ps_f.close()
hmm_init_f.close()
hmm_f.close()
probs_f.close()


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
