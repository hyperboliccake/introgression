import sys
import os
import predict
import gzip
import global_params as gp
from sim import sim_process
from align import align_helpers
from misc import read_fasta

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
    probs_f = gzip.open(probs_fn, 'rt')
    probs_temp_fn = probs_fn[:-7] + '.temp.txt.gz'
    probs_temp_f = gzip.open(probs_temp_fn, 'wt')
    line = probs_f.readline()
    probs_temp_f.write(line)
    print('reading already completed:')
    while line != '':
        x1 = line.find('\t')
        strain = line[:x1]
        x2 = line[x1+1:].find('\t') + x1 + 1
        chrm = line[x1+1:x2]
        if not completed.has_key(chrm):
            completed[chrm] = []
        completed[chrm].append(strain)
        print(strain, chrm)
        try:
            line = probs_f.readline()
            probs_temp_f.write(line)
        except Exception as e:
            print(e)
            break
    probs_f.close()
    probs_temp_f.close()
    os.system('mv ' + probs_temp_fn + ' ' + probs_fn)

    ps_f = gzip.open(ps_fn, 'rt')
    ps_temp_fn = ps_fn[:-7] + '.temp.txt.gz'
    ps_temp_f = gzip.open(ps_temp_fn, 'wt')
    line = ps_f.readline()
    ps_temp_f.write(line)
    while line != '':
        try:
            line = ps_f.readline()
            ps_temp_f.write(line)
        except Exception as e:
            print(e)
            break
    ps_f.close()
    ps_temp_f.close()
    os.system('mv ' + ps_temp_fn + ' ' + ps_fn)

    

write_ps = True
if write_ps:
    ps_f = gzip.open(ps_fn, open_mode + 't')

probs_f = gzip.open(probs_fn, open_mode + 't')

##======
# loop through all sequences and predict introgression
##======


for chrm in gp.chrms:

    for strain, strain_dir in strain_dirs:

        if resume and completed.has_key(chrm) and strain in completed[chrm]:
            print('already finished:', strain, chrm)
            continue

        print('working on:', strain, chrm)
        
        ref_prefix = '_'.join(gp.alignment_ref_order)
        fn = gp_dir + gp.alignments_dir + ref_prefix + '_' + strain + \
             '_chr' + chrm + '_mafft' + gp.alignment_suffix
        try:
            headers, seqs = read_fasta.read_fasta(fn)
        except Exception as e:
            print('no alignment for', strain, chrm)
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
        predict.write_state_probs(probs, probs_f, strain, chrm, hmm.hidden_states)

for k in blocks_f:
    blocks_f[k].close()
ps_f.close()
hmm_init_f.close()
hmm_f.close()
probs_f.close()
