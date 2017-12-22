import os
import sys
import copy
import re
import numpy.random
import itertools
sys.path.insert(0, '../hmm/')
import hmm_bw
sys.path.insert(0, '../sim/')
import sim_predict
import sim_process
sys.path.insert(0, '../')
import global_params as gp
sys.path.insert(0, '../misc')
import read_fasta

def process_predict_args(arg_list):

    d = {}
    i = 0

    d['tag'] = arg_list[i]
    i += 1

    d['improvement_frac'] = float(arg_list[i])
    i += 1

    d['threshold'] = float(arg_list[i])
    i += 1

    # expected length of introgressed tracts and fraction of sequence
    # introgressed
    expected_tract_lengths = {}
    expected_frac = {}

    species = []
    while i < len(arg_list):
        species.append(arg_list[i])
        i += 1
        expected_tract_lengths[species[-1]] = float(arg_list[i])
        i += 1
        expected_frac[species[-1]] = float(arg_list[i])
        i += 1
    d['species'] = species

    # TODO deal with noref

    expected_frac[species[0]] = 0
    expected_frac[species[0]] = 1 - sum(expected_frac.values())
    d['expected_frac'] = expected_frac

    # calculate these based on remaining bases, but after we know
    # which chromosome we're looking at
    expected_tract_lengths[species[0]] = 0    
    d['expected_tract_lengths'] = expected_tract_lengths
    d['expected_num_tracts'] = {}
    d['expected_bases'] = {}

    return d


def process_ref_args(fn):

    f = open(fn, 'r')
    line = f.readline()
    refs = {}
    while line != '':
        ref_species, ref_name, ref_seq_name, ref_seq_location = line.strip().split(' ')
        refs[ref_species] = (ref_name, ref_seq_name, ref_seq_location)
        line = f.readline()
    f.close()
    return refs

def process_strain_args(fn):
    f = open(fn, 'r')
    line = f.readline()
    strains = {}
    while line != '':
        line = line.strip().split(' ')
        species = line[0]
        strain_seq_location = line[1]
        strain_names = line[2:]
        if not strains.has_key(species):
            strains[species] = []
        for strain in strain_names:
            strains[species].append((strain, strain_seq_location))
        line = f.readline()
    f.close()
    return strains

def process_args(arg_list):

    refs = process_ref_args(arg_list[1])
    strains = process_strain_args(arg_list[2])
    args = process_predict_args(arg_list[3:])
    return refs, strains, args

'''
def read_seq(name, path, chrm):

    name, path = strain
    fn = path + '/' + name + '/' + name + 'chr' + chrm + '.fa'
    headers, seqs = read_fasta.read_fasta(fn)
    return seqs[0]

def read_ref_seqs(refs, chrm):

    ref_seqs = []
    for species in refs:
        seq = read_seq(refs[species][1], refs[species][2], chrm)
        ref_seqs.append(seq)
    return ref_seqs
'''

def read_aligned_seqs(fn, refs, strain, species_order):
    headers, seqs = read_fasta.read_fasta(fn)
    d = {}
    for i in range(len(seqs)):
        name = headers[i][1:].split(' ')[0]
        d[name] = seqs[i]

    ref_seqs = []
    for s in species_order:
        ref_seqs.append(d[refs[s][0]])
    predict_seq = d[strain]

    return ref_seqs, predict_seq

def set_expectations(args, n):

    species_to = args['species'][0]
    species_from = copy.deepcopy(args['species'])
    species_from.remove(species_to)

    for s in species_from:
        args['expected_num_tracts'][s] = \
            args['expected_frac'][s] * n / args['expected_tract_lengths'][s]
        args['expected_bases'][s] = args['expected_num_tracts'][s] * \
                                 args['expected_tract_lengths'][s]

    args['expected_bases'][species_to] = \
        n - sum([args['expected_bases'][s] for s in species_from])

    args['expected_num_tracts'][species_to] = \
        sum([args['expected_num_tracts'][s] for s in species_from]) + 1

    args['expected_tract_lengths'][species_to] = \
        args['expected_bases'][species_to] / args['expected_num_tracts'][species_to]

def ungap_and_code_helper(predict_seq, ref_seqs, index_ref):

    # assume the first reference is what we want to index from
    # and ref seqs are in general in correct order

    ps = []
    seq = []
    ind = 0
    for i in range(len(predict_seq)):
        xi = predict_seq[i]
        # only keep this position in sequence if no gaps in
        # the alignment column
        keep = True
        if xi == gp.gap_symbol:
            keep = False
        # determine which references the sequence matches at
        # this position
        symbol = ''
        for r in range(len(ref_seqs)):
            ri = ref_seqs[r][i]
            if ri == gp.gap_symbol:
                keep = False
                break
            elif xi == ri:
                symbol += gp.match_symbol
            else:
                symbol += gp.mismatch_symbol
        if keep:
            seq.append(symbol)
            ps.append(ind)
        if ref_seqs[index_ref][i] != gp.gap_symbol:
            ind += 1
    return seq, ps
        
def ungap_and_code(predict_seq, ref_seqs, index_ref = 0):
    
    seq, ps = ungap_and_code_helper(predict_seq, ref_seqs, index_ref)
    ref_seqs_coded = []
    for r in ref_seqs:
        ref_seq_coded, ps_r = ungap_and_code_helper(r, ref_seqs, index_ref)
        ref_seqs_coded.append(ref_seq_coded)
    return seq, ref_seqs_coded, ps

def predict_introgressed(ref_seqs, predict_seq, predict_args, \
                         train=True, method='posterior'):

    # code sequence by which reference it matches at each site
    seq_coded, ref_seqs_coded, ps = ungap_and_code(predict_seq, ref_seqs)

    set_expectations(predict_args, len(predict_seq))

    # initial values for initial, emission, and transition
    # probabilities
    predict_args['states'] = predict_args['species'] # hackity
    predict_args['unknown_species'] = None # hack
    predict_args['ref_inds'] = [1,2] # hack
    seqs_for_hmm_params = [ref_seqs_coded[0]] + [seq_coded] + ref_seqs_coded[1:]
    species_to_indices = {predict_args['species'][0]:[0,1]}
    for i in range(1, len(predict_args['species'])):
        species_to_indices[predict_args['species'][i]] = [1+i]
    init, emis, trans = sim_predict.initial_hmm_parameters(seqs_for_hmm_params, \
                                                           species_to_indices, \
                                                           predict_args['species'][0], \
                                                           len(seq_coded), \
                                                           predict_args)

    ######
    # make predictions
    ######

    default_state = predict_args['species'][0]

    # new Hidden Markov Model
    hmm = hmm_bw.HMM()

    # set obs
    hmm.set_obs([seq_coded])

    # set states and initial probabilties
    hmm.set_states(predict_args['species'])
    hmm.set_init(init)
    hmm.set_emis(emis)
    hmm.set_trans(trans)

    hmm_init = copy.deepcopy(hmm)

    # optional Baum-Welch parameter estimation
    if train:
        hmm.go(predict_args['improvement_frac'])

    if method == "posterior":
        predicted = {}
        all_probs = {}
        # for all obs sequences, each site is a dic with one prob for each
        # state
        p = hmm.posterior_decoding()
        path, path_probs = sim_process.get_max_path(p[0])
        path_t = sim_process.threshold_predicted(path, path_probs, \
                                                 predict_args['threshold'], \
                                                 default_state)

        return path_t, p[0], hmm, hmm_init, ps
        
    if method == 'viterbi':
        hmm.set_obs(predict_seq)
        predicted = convert_predictions(hmm.viterbi(), predict_args['states'])
        return predicted, hmm, hmm_init, ps

    else:
        print 'invalid method'

def write_positions(ps, f, strain, chrm):
    sep = '\t'
    f.write(strain + sep + chrm + sep + sep.join([str(x) for x in ps]) + '\n')
    f.flush()

def read_positions(fn):
    # dictionary keyed by strain and then chromsome
    f = open(fn, 'r')
    line = f.readline()
    d = {}
    while line != '':
        line = line.strip().split('\t')
        strain = line[0]
        chrm = line[1]
        ps = [int(x) for x in line[2:]]
        if not d.has_key(strain):
            d[strain] = {}
        d[strain][chrm] = ps
        line = f.readline()
    f.close()
    return d

def write_blocks_header(f):
    sep = '\t'
    f.write('strain' + sep + 'chromosome' + sep + 'predicted_species' + sep + \
            'start' + sep + 'end' + sep + 'number_non_gap' + '\n')
    f.flush()

def write_blocks(state_seq_blocks, ps, f, strain, chrm, species_pred):
    # one file for each species
    # file format is:
    # strain chrm predicted_species start end number_non_gap
    sep = '\t'
    for block in state_seq_blocks:
        start, end = block
        f.write(strain + sep + chrm + sep + species_pred + sep + \
                str(ps[start]) + sep + str(ps[end]) + sep + \
                str(end - start + 1) + '\n')
    f.flush()

def read_blocks(fn):
    # return dictionary of (start, end, number_non_gap), keyed by strain
    # and then chromosome
    f = open(fn, 'r')
    line = f.readline()
    d = {}
    while line != '':
        strain, chrm, species, start, end, number_non_gap = line.strip().split('\t')
        if not d.has_key(strain):
            d[strain] = {}
        if not d[strain].has_key(chrm):
            d[strain][chrm] = []
        d[strain][chrm].append((int(start), int(end), int(number_non_gap)))
        line = f.readline()
    f.close()
    return d

def write_hmm_header(states, f):

    sep = '\t'
    header_string = 'strain' + sep + 'chromosome' + sep

    symbols = [gp.match_symbol, gp.mismatch_symbol, gp.unknown_symbol]
    emis_symbols = [''.join(x) for x in \
                           list(itertools.product(symbols, repeat=len(states)))]
    emis_symbols.sort()

    # initial
    for s in states:
        header_string += 'init_' + s + sep

    # emission
    for s in states:
        for symbol in emis_symbols:
            header_string += 'emis_' + s + '_' + symbol + sep

    # transition
    for s1 in states:
        for s2 in states:
            header_string += 'trans_' + s1 + '_' + s2 + sep
    
    f.write(header_string[:-(len(sep))] + '\n')
    f.flush()

def write_hmm(hmm, f, strain, chrm):
    
    sep = '\t'
    emis_symbols = sorted(hmm.emis[0].keys())

    line_string = strain + sep + chrm + sep

    # initial
    for i in range(len(hmm.states)):
        line_string += str(hmm.init[i]) + sep

    # emission
    for i in range(len(hmm.states)):
        for symbol in emis_symbols:
            line_string += str(hmm.emis[i][symbol]) + sep

    # transition
    for i in range(len(hmm.states)):
        for j in range(len(hmm.states)):
            line_string += str(hmm.trans[i][j]) + sep

    f.write(line_string[:-(len(sep))] + '\n')
    f.flush()

def write_state_probs(probs, f, strain, chrm):

    # probs is list of sites, each site dic keyed by state
    
    # file format is:
    # strain\tchrm\tcer:.1,.2,.3\tpar:.9,.8,.7
    # strain\tchrm\tcer:.4,.2,.3\tpar:.6,.8,.7

    sep  ='\t'
    f.write(chrm + sep + chrm)

    for state in probs[0].keys():
        f.write(sep + state + ':')
        probs_string = ','.join(["{0:.5f}".format(site[state]) \
                                     for site in probs])
        f.write(probs_string)
    f.write('\n')
    f.flush()


"""


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

"""
