import os
import sys
import copy
import re
import numpy.random
sys.path.insert(0, '../hmm/')
from hmm_bw import *
sys.path.insert(0, '../sim/')
import sim_analyze_hmm_bw as sim
sys.path.insert(0, '../')
import global_params as gp

resume = False

def get_seqs_chrm(x, refs, chrm, match_symbol, mismatch_symbol, \
                      unknown_symbol, unsequenced_symbol, master_ref):
    # for now, just make each alignment block a separate sequence, but
    # maybe check at some point whether any of them should be
    # concatenated?

    # what to do about columns with gaps? would need a completely
    # separate model for predicting indels...so just ignore all
    # columns with gaps (this would include things where not all of
    # the references aligned)...but we don't want to throw off the
    # locations of sites

    mult_expected = len(refs) + 1
    seqs = []
    seq_inds = []
    scores = []
    psx = []
    fn = ''
    for r in gp.alignment_ref_order:
        if r in refs:
            fn += r + '_'
    fn += x + '_chr' + chrm + '.maf'
    f = open(gp.alignments_dir + fn, 'r')
    mult = -1
    label = None
    line = f.readline()
    while line != '':
        if line[0] == '#' or line[0] == '\n':
            line = f.readline()
            continue

        if line[0] == 'a':
            m = re.search('score=(?P<score>[0-9\.]+)', line)
            scores.append(float(m.group('score')))
            m = re.search('mult=(?P<mult>[0-9\.]+)', line)
            mult = int(m.group('mult'))
            m = re.search('label=(?P<label>[A-Za-z0-9_]+)', line)
            label = m.group('label')
            line = f.readline()
            
        if line[0] == 's':
            # only want sections aligned to _all_ references
            if mult < mult_expected:
                line = f.readline()
                continue

            seqx = ''
            startx = -1
            sizex = -1
            strandx = ''
            strand_master_ref = ''
            ref_seqs = {}
            current_master_ref_ind = -1
            for i in range(mult_expected):
                s, src, start, size, strand, src_size, text = line.split()
                src = src[:src.index('_')]
                if src == x:
                    seqx = text
                    # leave everything indexed from 1 instead of 0
                    startx = int(start)
                    sizex = int(size)
                    strandx = strand
                else:
                    if src == master_ref:
                        current_master_ref_ind = int(start)
                        strand_master_ref = strand
                    ref_seqs[src] = text
                line = f.readline()

            seq = []
            inds = []
            for i in range(len(seqx)):
                xi = seqx[i]
                # only keep this position in sequence if no gaps in
                # the alignment column
                keep = True
                if xi == gp.gap_symbol:
                    keep = False
                # determine which references the sequence matches at
                # this position
                symbol = ''
                for ref in refs:
                    ri = ref_seqs[ref][i]
                    if ri == gp.gap_symbol:
                        keep = False
                        break
                    elif xi == ri:
                        symbol += match_symbol
                    else:
                        symbol += mismatch_symbol
                if keep:
                    seq.append(symbol)
                    inds.append(current_master_ref_ind)
                if ref_seqs[master_ref][i] != gp.gap_symbol:
                    current_master_ref_ind += 1

            if len(seq) > 0:
                seqs.append(seq)
                seq_inds.append(inds)
                psx.append((x, chrm, label, strand_master_ref))

        else:
            line = f.readline()

    f.close()

    return seqs, seq_inds, psx

def get_seqs(x, refs, chrms, match_symbol, mismatch_symbol, \
                 unknown_symbol, unsequenced_symbol, master_ref):
    all_seqs = []
    all_seq_inds = []
    all_psx = []
    for chrm in chrms:
        seqs, seq_inds, psx = get_seqs_chrm(x, refs, chrm, \
                                                match_symbol, mismatch_symbol, \
                                                unknown_symbol, unsequenced_symbol, \
                                                master_ref)
        all_seqs += seqs
        all_seq_inds += seq_inds
        all_psx += psx

    return all_seqs, all_seq_inds, all_psx

def convert_predictions(path, states):
    new_path = []
    for p in path:
        new_path.append(states[p])
    return new_path

def read_hmm_params(fn, states, sim_states, unknown_state):

    # read parameters into dictionary first
    sim_states_to_states = dict(zip(sim_states, states))
    lines = [line.strip().split('\t') for line in open(fn, 'r').readlines()]
    d = {'init':{}, 'emis':{}, 'trans':{}}
    for line in lines:
        if line[0] == 'init':
            state = sim_states_to_states[line[1]]
            d['init'][state] = float(line[2])
        # this line will look like e.g. emis cer cer + par + bay - .2
        elif line[0] == 'emis':
            state = sim_states_to_states[line[1]]
            if state not in d['emis']:
                d['emis'][state] = {}
            symbol_list = []
            len_symbol_list = len(sim_states)
            if unknown_state:
                len_symbol_list -= 1
            for i in range(len_symbol_list):
                symbol_list.append(sim_states_to_states[line[2*i+2]])
                symbol_list.append(line[2*i+3])
            d['emis'][state][tuple(symbol_list)] = float(line[-1])
        else:
            assert line[0] == 'trans'
            from_state = sim_states_to_states[line[1]]
            to_state = sim_states_to_states[line[2]]
            if from_state not in d['trans']:
                d['trans'][from_state] = {}
            d['trans'][from_state][to_state] = float(line[3])

    # init
    init = []
    for state in states:
        init.append(d['init'][state])

    # emis
    refs = copy.deepcopy(states)
    if unknown_state:
        refs = refs[:-1]
    emis = []
    for state in states:
        emis_state = {}
        for symbol_list in d['emis'][state]:
            symbol = ''
            for ref in refs:
                x = symbol_list.index(ref) + 1
                symbol += symbol_list[x]
            emis_state[symbol] = d['emis'][state][symbol_list]
        emis.append(emis_state)

    # trans
    trans = []
    for state_from in states:
        row = []
        for state_to in states:
            row.append(d['trans'][state_from][state_to])
        trans.append(row)

    return init, emis, trans

def predict_introgressed_hmm(seqs, states, sim_states, unknown_state, init, emis, trans):

    # create a hidden markov model and determine which reference genome we
    # are most likely to be in at each variant site
    hmm = HMM()


    hmm.set_states(states)

    hmm.set_init(init)
    hmm.set_emis(emis)
    hmm.set_trans(trans)

    predicted = []

    sys.stdout.flush()

    for i in range(len(seqs)):
        if i % 100 == 0:
            print 'viterbi on seq', i
            sys.stdout.flush()
        hmm.set_obs(seqs[i])
        predicted.append(convert_predictions(hmm.viterbi(), states))

    return predicted, hmm

'''
def predict_introgressed_id(seqs):
    # ok so the aligned sequences can have gaps
    window_size = 1000
    window_shift = 500
    thresholds = [.7, .95, .96]
    predicted = []
    for seq in seqs:
        p = [0 for b in xrange(len(seq))]
        for i in range(0, len(seq) - window_size, window_shift):
            region = seq[i:i+window_size]
            count0 = region.count('0')
            count1 = region.count('1')
            count2 = region.count('2')
            #count3 = window_size - count0 - count1 - count2
            cer_match = float(count0 + count2) / len(region)
            par_match = float(count1 + count2) / len(region) 
            if cer_match > thresholds[0]:
                if cer_match < thresholds[1] and par_match > thresholds[2]:
                    # this automatically appends to the list if we've
                    # gone over the length, but we'll deal with that
                    # later
                    p[i:i+window_size] = [1] * window_size 
        predicted.append(p[:len(seq)])

    return predicted
'''

def get_predicted_tracts(predicted, all_ps, master_ref, seq_inds):

    blocks = []
    # loop through all sequences
    for s in xrange(len(predicted)):
        print 'alignment block: ', all_ps[s]
        alignment_block = all_ps[s]
        # loop through all sites in that sequence, keeping track of
        # introgressed blocks
        prev_species = None
        start_int = 0
        end_int = 0
        predicted[s] = predicted[s] + ['END']
        for i in xrange(len(predicted[s])):
            # staying in same species
            if predicted[s][i] == prev_species:
                end_int = i
            # ending previous block
            else:
                if prev_species != None and prev_species != master_ref:
                    blocks.append(list(alignment_block) + \
                                      [prev_species, \
                                           seq_inds[s][start_int], seq_inds[s][end_int], \
                                           end_int - start_int + 1])
                start_int = i
                end_int = i
                prev_species = predicted[s][i]
                
    return blocks

def write_predicted_tracts(blocks, f):
    for block in blocks:
        block = [str(x) for x in block]
        f.write('\t'.join(block) + '\n')

