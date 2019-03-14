import sys
import os
import gzip
import predict
import global_params as gp
from sim import sim_process
from align import align_helpers
from misc import read_fasta

# read in analysis parameters

args = predict.process_predict_args(sys.argv[1:])

strain_dirs = align_helpers.get_strains(
    align_helpers.flatten(gp.non_ref_dirs.values()))

# output files and if and where to resume

if not os.path.isdir(gp.analysis_out_dir_absolute + args['tag']):
    os.makedirs(gp.analysis_out_dir_absolute + args['tag'])

# positions
# TODO move this to more general location and make separate files for
# each strain x chrm
base_dir = f'{gp.analysis_out_dir_absolute}{args["tag"]}'

# introgressed blocks
blocks_f = {}
for s in args['states']:
    blocks_f[s] = open(f'{base_dir}/blocks_{s}_{args["tag"]}.txt', 'w')
    predict.write_blocks_header(blocks_f[s])

# HMM parameters
emis_symbols = predict.get_emis_symbols(args['known_states'])

hmm_init_f = open(f'{base_dir}/hmm_init_{args["tag"]}.txt', 'w')
predict.write_hmm_header(args['known_states'], args['unknown_states'],
                         emis_symbols, hmm_init_f)

hmm_f = open(f'{base_dir}/hmm_{args["tag"]}.txt', 'w')
predict.write_hmm_header(args['known_states'], args['unknown_states'],
                         emis_symbols, hmm_f)

# posterior probabilities

write_ps = True
if write_ps:
    ps_f = gzip.open(f'{base_dir}/positions_{args["tag"]}.txt', 'wt')

probs_f = gzip.open(f'{base_dir}/probs_{args["tag"]}.txt', 'wt')

# loop through all sequences and predict introgression


for chrm in gp.chrms:

    for strain, strain_dir in strain_dirs:

        print(f'working on: {strain} {chrm}')

        ref_prefix = '_'.join(gp.alignment_ref_order)
        fn = (f'{gp.alignments_dir}{ref_prefix}_{strain}'
              f'_chr{chrm}_mafft{gp.alignment_suffix}')

        if not os.path.exists(fn):
            print(fn)
            print(f'no alignment for {strain} {chrm}')
            continue

        headers, seqs = read_fasta.read_fasta(fn)

        ref_seqs = seqs[:-1]
        predict_seq = seqs[-1]

        # predict introgressed/non-introgressed tracts

        state_seq, probs, hmm, hmm_init, ps = \
            predict.predict_introgressed(ref_seqs, predict_seq,
                                         args, train=True)

        state_seq_blocks = sim_process.convert_to_blocks_one(state_seq,
                                                             args['states'])
        # output

        # the positions actually used in predictions
        # (alignment columns with no gaps)
        if write_ps:
            predict.write_positions(ps, ps_f, strain, chrm)

        # blocks predicted to be introgressed, separate files for each species
        for s in state_seq_blocks:
            predict.write_blocks(state_seq_blocks[s], ps, blocks_f[s],
                                 strain, chrm, s)

        # summary info about HMM (before training)
        predict.write_hmm(hmm_init, hmm_init_f, strain, chrm, emis_symbols)

        # summary info about HMM (after training)
        predict.write_hmm(hmm, hmm_f, strain, chrm, emis_symbols)

        # probabilities at each site
        predict.write_state_probs(probs, probs_f, strain,
                                  chrm, hmm.hidden_states)

for k in blocks_f:
    blocks_f[k].close()

ps_f.close()
hmm_init_f.close()
hmm_f.close()
probs_f.close()
