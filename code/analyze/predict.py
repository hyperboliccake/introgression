import copy
import gzip
import itertools
from collections import defaultdict, Counter
from hmm import hmm_bw
from sim import sim_predict
from sim import sim_process
import global_params as gp
from misc import read_fasta
import numpy as np

def read_aligned_seqs(fn, strain):
    headers, seqs = read_fasta.read_fasta(fn)
    d = {}
    for i in range(len(seqs)):
        name = headers[i][1:].split(' ')[0]
        d[name] = seqs[i]

    ref_seqs = []
    for ref in gp.alignment_ref_order:
        ref_seqs.append(d[ref])
    predict_seq = d[strain]

    return ref_seqs, predict_seq


def set_expectations(args, n):

    species_to = gp.alignment_ref_order[0]
    species_from = gp.alignment_ref_order[1:]

    args['expected_num_tracts'] = {}
    args['expected_bases'] = {}
    for s in species_from:
        args['expected_num_tracts'][s] = \
            args['expected_frac'][s] * n / args['expected_length'][s]
        args['expected_bases'][s] = args['expected_num_tracts'][s] * \
                                    args['expected_length'][s]

    args['expected_bases'][species_to] = \
        n - sum([args['expected_bases'][s] for s in species_from])

    args['expected_num_tracts'][species_to] = \
        sum([args['expected_num_tracts'][s] for s in species_from]) + 1

    args['expected_length'][species_to] = \
        args['expected_bases'][species_to] / args['expected_num_tracts'][species_to]

def ungap_and_code_helper(predict_seq, ref_seqs, index_ref):

    # index_ref is index of reference strain to index relative to 

    ps = [] # positions that we're keeping in analysis
    seq = [] # sequence of emitted symbols, e.g. '++-++'
    ind = 0 # current position in reference strain
    
    # TODO: move positions file to more general location (maybe? but
    # it can change based on references so unless there's a deeper
    # organizational structure it needs to be for the specific run),
    # and here first see if we have this info stored in a file already

    for i in range(len(predict_seq)):
        xi = predict_seq[i]
        # only keep this position in sequence if no gaps in
        # the alignment column
        keep = True
        if xi == gp.gap_symbol or xi == gp.unsequenced_symbol:
            keep = False
        # determine which references the sequence matches at
        # this position
        symbol = ''
        for r in range(len(ref_seqs)):
            ri = ref_seqs[r][i]
            if ri == gp.gap_symbol or ri == gp.unsequenced_symbol:
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


def poly_sites(seq, ref_seqs, ps):
    ps_poly = []
    seq_poly = []
    ref_seqs_poly = [[] for r in ref_seqs]

    for i in range(len(ps)):
        if set(seq[i]) != set([gp.match_symbol]):
            ps_poly.append(ps[i])
            seq_poly.append(seq[i])
            for j in range(len(ref_seqs)):
                ref_seqs_poly[j].append(ref_seqs[j][i])
    
    return seq_poly, ref_seqs_poly, ps_poly

def get_symbol_freqs(seq):

    num_states = len(seq[0])
    num_sites = len(seq)

    individual_symbol_freqs = []
    for s in range(num_states):
        d = defaultdict(int)
        for i in range(num_sites):
            d[seq[i][s]] += 1
        for sym in d:
            d[sym] /= float(num_sites)
        individual_symbol_freqs.append(d)

    symbol_freqs = defaultdict(int)
    for i in range(num_sites):
        symbol_freqs[seq[i]] += 1
    for sym in symbol_freqs:
        symbol_freqs[sym] /= float(num_sites)

    # for each state, how often seq matches that state relative to
    # others
    weighted_match_freqs = []
    for s in range(num_states):
        weighted_match_freqs.append(individual_symbol_freqs[s][gp.match_symbol])
    weighted_match_freqs = norm_list(weighted_match_freqs)

    return individual_symbol_freqs, symbol_freqs, weighted_match_freqs

def norm_list(l):
    scale = float(sum(l))
    for i in range(len(l)):
        l[i] /= scale
    return l
            

def norm_dict(d):
    scale = float(sum(d.values()))
    for k in d:
        d[k] /= scale
    return d
    args['expected_length'][species_to] = \
        args['expected_bases'][species_to] /\
        args['expected_num_tracts'][species_to]


def ungap_and_code(predict_seq, ref_seqs, index_ref=0):
    # index_ref is index of reference strain to index relative to
    # build character array
    sequences = np.array([list(predict_seq)] +
                         [list(r) for r in ref_seqs])

    # make boolean for valid characters
    isvalid = np.logical_and(sequences != gp.gap_symbol,
                             sequences != gp.unsequenced_symbol)

    # positions are where everything is valid, index where the reference is
    # valid.  The +1 removes the predict sequence at index 0
    positions = np.where(
        np.all(isvalid[:, isvalid[index_ref+1, :]], axis=0))[0]

    matches = np.where(sequences[0] == sequences[1:],
                       gp.match_symbol,
                       gp.mismatch_symbol)

    # 1: indexing removes currently examined sequence
    matches = [''.join(row)
               for row in np.transpose(matches[:, np.all(isvalid, axis=0)])]

    # NOTE list is for unit test comparisons
    return matches, positions


def poly_sites(sequences, positions):
    seq_len = len(sequences[0])
    # check if seq only contains match_symbol
    retain = np.vectorize(
        lambda x: x.count(gp.match_symbol) != seq_len)(sequences)
    indices = np.where(retain)[0]
    ps_poly = [positions[i] for i in indices]
    seq_poly = [sequences[i] for i in indices]
    return seq_poly, ps_poly


def get_symbol_freqs(sequence):

    individual = []
    weighted = []

    symbols = defaultdict(int, Counter(sequence))
    total = len(sequence)
    for k in symbols:
        symbols[k] /= total

    sequence = np.array([list(s) for s in sequence])

    # look along species
    for s in np.transpose(sequence):
        s = ''.join(s)
        counts = Counter(s)
        weighted.append(counts[gp.match_symbol])
        total = sum(counts.values())
        for k in counts:
            counts[k] /= total
        individual.append(defaultdict(int, counts))

    total = sum(weighted)
    weighted = [w / total for w in weighted]
    return individual, symbols, weighted


def initial_probabilities(known_states, unknown_states,
                          expected_frac, weighted_match_freqs):

    init = []
    expectation_weight = .9
    for s, state in enumerate(known_states):
        expected = expected_frac[state]
        estimated = weighted_match_freqs[s]
        init.append(expected * expectation_weight +
                    estimated * (1 - expectation_weight))

    for state in unknown_states:
        expected_frac = expected_frac[state]
        init.append(expected_frac)

    return init / np.sum(init)


def emission_probabilities(known_states, unknown_states, symbols):

    probabilities = {
        gp.mismatch_symbol + gp.match_symbol: 0.9,
        gp.match_symbol + gp.match_symbol: 0.09,
        gp.mismatch_symbol + gp.mismatch_symbol: 0.009,
        gp.match_symbol + gp.mismatch_symbol: 0.001,
    }

    mismatch_bias = .99

    for s in range(len(unknown_states)):
        state = unknown_states[s]
        emis.append(defaultdict(float))
        for symbol in symbols:
            match_count = symbol.count(gp.match_symbol)
            mismatch_count = symbol.count(gp.mismatch_symbol)
            emis[s + len(known_states)][symbol] = (match_count * (1 - mismatch_bias) + \
                                                    mismatch_count * mismatch_bias)
        emis[s + len(known_states)] = norm_dict(emis[s + len(known_states)])

    return emis

def emission_probabilities_old(known_states, unknown_states, symbol_freqs):

    # doesn't use prior expectations, but probably should if we ever
    # want to include multiple unknown states

    own_bias = .99
    
    emis = []
    for s in range(len(known_states)):
        state = known_states[s]
        emis.append(defaultdict(float))
        for symbol in symbol_freqs:
            f = symbol_freqs[symbol]
            match = symbol[s] == gp.match_symbol
            if match:
                emis[s][symbol] = f * own_bias
            else:
                emis[s][symbol] = f * (1 - own_bias)
        emis[s] = norm_dict(emis[s])
    for s in range(len(unknown_states)):
        state = unknown_states[s]
        emis.append(defaultdict(float))
        for symbol in symbol_freqs:
            f = symbol_freqs[symbol]
            match_count = symbol.count(gp.match_symbol)
            mismatch_count = symbol.count(gp.mismatch_symbol)
            emis[s + len(known_states)][symbol] = f
            emis[s + len(known_states)][symbol] *= (match_count * (1 - own_bias) + \
                                                    mismatch_count * own_bias)
        emis[s + len(known_states)] = norm_dict(emis[s + len(known_states)])

    return emis

def transition_probabilities(known_states, unknown_states, \
                             expected_frac, expected_length):

    num_per_category = 2 ** (len(known_states) - 2)
    for key in probabilities:
        probabilities[key] *= num_per_category

    # for known states
    symbol_array = np.array([list(s) for s in symbols], dtype='<U1')
    # for unknown states
    symbol_length = symbol_array.shape[1]
    number_matches = (symbol_array == gp.match_symbol).sum(axis=1)
    # combine first column with rest to generate probabilities
    first_column = np.tile(symbol_array[:, 0:1], (1, len(known_states)))
    symbol_array = np.core.defchararray.add(
        first_column, symbol_array[:, 0:len(known_states)])
    # index into probabilities and normalize
    emissions = np.vectorize(probabilities.__getitem__)(symbol_array)
    emissions /= sum(emissions)

    # convert to match * (1-bias) + mismatch * bias, simplified
    number_matches = (number_matches + mismatch_bias *
                      (symbol_length - 2 * number_matches))
    number_matches /= sum(number_matches)
    # repeat for each unknown state
    number_matches = np.transpose(
        np.tile(number_matches, (len(unknown_states), 1)))

    result = [defaultdict(float,
                          {k: v for k, v in
                           zip(symbols, emissions[:, i])})
              for i in range(emissions.shape[1])]
    result.extend([defaultdict(float,
                               {k: v for k, v in
                                zip(symbols, number_matches[:, i])})
                   for i in range(number_matches.shape[1])])

    return result


def transition_probabilities(known_states, unknown_states,
                             expected_frac, expected_tract_lengths):

    # doesn't depend on sequence observations but maybe it should?

    # also should we care about number of tracts rather than fraction
    # of genome? maybe theoretically, but that number is a lot more
    # suspect

    states = known_states + unknown_states

    trans = []
    for i in range(len(states)):
        state_from = states[i]
        trans.append([])
        scale_other = 1 / (1 - expected_frac[state_from])
        for j in range(len(states)):
            state_to = states[j]
            if state_from == state_to:
                trans[i].append(1 - 1./expected_length[state_from])
            else:
                trans[i].append(1./expected_length[state_from] * \
                                expected_frac[state_to] * scale_other)

    fractions = np.array([expected_frac[s] for s in states])
    lengths = 1/np.array([expected_tract_lengths[s] for s in states])

    # general case,
    # trans[i,j] = 1/ length[i] * expected[j] * 1 /(1 - fraction[i])
    transitions = np.outer(
        np.multiply(lengths, 1/(1-fractions)),
        fractions)
    # when i == j, trans[i,j] = 1 - 1/length[i]
    np.fill_diagonal(transitions, 1-lengths)

    # normalize
    return transitions / transitions.sum(axis=1)[:, None]

def initial_hmm_parameters(seq, known_states, unknown_states, \
                           expected_frac, expected_length):

    # get frequencies of individual symbols (e.g. '+') and all full
    # combinations of symbols (e.g. '+++-')
    individual_symbol_freqs, symbol_freqs, weighted_match_freqs = get_symbol_freqs(seq)

    init = initial_probabilities(known_states, unknown_states,
                                 expected_frac, weighted_match_freqs)
    emis = emission_probabilities(known_states, unknown_states, symbol_freqs.keys())

    trans = transition_probabilities(known_states, unknown_states, \
                                     expected_frac, expected_length)

    # new Hidden Markov Model
    hmm = hmm_bw.HMM()

    hmm.set_initial_p(init)
    hmm.set_emissions(emis)
    hmm.set_transitions(trans)
    return hmm


def predict_introgressed(ref_seqs, predict_seq, predict_args,
                         train=True, only_poly_sites=True,
                         return_positions=False):

    # code sequence by which reference it matches at each site
    seq_coded, positions = ungap_and_code(predict_seq, ref_seqs)
    if only_poly_sites:
        seq_coded, positions = poly_sites(seq_coded, positions)
    if return_positions:
        return positions

    # sets expected number of tracts and bases for each reference
    # based on expected length of introgressed tracts and expected
    # total fraction of genome
    set_expectations(predict_args, len(predict_seq))

    # set initial hmm parameters based on combination of (1) initial
    # expectations (length of introgressed tract and fraction of
    # genome/total number tracts and bases) and (2) number of sites at
    # which predict seq matches each reference

    hmm = initial_hmm_parameters(seq_coded,
                                 predict_args['known_states'],
                                 predict_args['unknown_states'],
                                 predict_args['expected_frac'],
                                 predict_args['expected_length'])

    # make predictions

    # set states and initial probabilties
    hmm.set_hidden_states(predict_args['states'])

    # copy before setting observations to save memory
    hmm_init = copy.deepcopy(hmm)

    # set obs
    hmm.set_observations([seq_coded])

    # optional Baum-Welch parameter estimation
    if train:
        hmm.train(predict_args['improvement_frac'])

    p = hmm.posterior_decoding()
    path, path_probs = sim_process.get_max_path(p[0], hmm.hidden_states)

    # posterior
    if type(predict_args['threshold']) is float:
        path_t = sim_process.threshold_predicted(path, path_probs,
                                                 predict_args['threshold'],
                                                 predict_args['states'][0])
        return path_t, p[0], hmm, hmm_init, positions

    else:
        hmm.set_observations([seq_coded])
        predicted = sim_predict.convert_predictions(hmm.viterbi(),
                                                    predict_args['states'])
        return predicted, p[0], hmm, hmm_init, positions


def write_positions(ps, writer, strain, chrm):
    writer.write(f'{strain}\t{chrm}\t' +
                 '\t'.join([str(x) for x in ps]) + '\n')


def convert_to_blocks(state_seq, states):
    # single individual state sequence
    blocks = {}
    for state in states:
        blocks[state] = []
    prev_species = state_seq[0]
    block_start = 0
    block_end = 0
    for i in range(len(state_seq)):
        if state_seq[i] == prev_species:
            block_end = i
        else:
            blocks[prev_species].append((block_start, block_end))
            block_start = i
            block_end = i
            prev_species = state_seq[i]
    # add last block
    if prev_species not in blocks:
        blocks[prev_species] = []
    blocks[prev_species].append((block_start, block_end))

    return blocks

def write_positions(ps, f, strain, chrm):
    sep = '\t'
    f.write(strain + sep + chrm + sep + sep.join([str(x) for x in ps]) + '\n')
    f.flush()

def read_positions(fn):
    # dictionary keyed by strain and then chromosome
    with gzip.open(fn, 'rb') as reader:
        result = defaultdict({})
        for line in reader:
            line = line.split()
            strain, chrm = line[0:2]
            ps = [int(x) for x in line[2:]]
            result[strain][chrm] = ps
    return result


def write_blocks_header(writer):
    # NOTE: num_sites_hmm represents the sites considered by the HMM,
    # so it might exclude non-polymorphic sites in addition to gaps
    writer.write('\t'.join(['strain',
                            'chromosome',
                            'predicted_species',
                            'start',
                            'end',
                            'num_sites_hmm'])
                 + '\n')


# TODO: find source of all the newlines in output!!
def write_blocks(state_seq_blocks, ps, writer, strain, chrm, species_pred):
    # file format is:
    # strain chrm predicted_species start end number_non_gap
    writer.write('\n'.join(
        ['\t'.join([strain,
                    chrm,
                    species_pred,
                    str(ps[start]),
                    str(ps[end]),
                    str(end - start + 1)])
         for start, end in state_seq_blocks]))
    if state_seq_blocks:
        writer.write('\n')


def read_blocks(fn, labeled=False):
    # return dictionary of (start, end, number_non_gap, [region_id]),
    # keyed by strain and then chromosome
    with open(fn, 'r') as reader:
        reader.readline()  # header
        result = defaultdict(lambda: defaultdict(list))
        for line in reader:
            tokens = line.split()
            if labeled:
                (region_id, strain, chrm, species,
                 start, end, number_non_gap) = tokens
                item = (region_id, int(start), int(end), int(number_non_gap))
            else:
                (strain, chrm, species,
                 start, end, number_non_gap) = tokens
                item = (int(start), int(end), int(number_non_gap))
            result[strain][chrm].append(item)
    return result


def get_emis_symbols(known_states):

    symbols = [gp.match_symbol, gp.mismatch_symbol]
    emis_symbols = [''.join(x) for x in
                    itertools.product(symbols, repeat=len(known_states))]
    emis_symbols.sort()
    return emis_symbols


def write_hmm_header(known_states, unknown_states, symbols, writer):

    writer.write('strain\tchromosome\t')

    states = known_states + unknown_states

    writer.write('\t'.join(
        [f'init_{s}' for s in states] +  # initial
        [f'emis_{s}_{symbol}'
         for s in states
         for symbol in symbols] +  # emissions
        [f'trans_{s1}_{s2}'
         for s1 in states
         for s2 in states]))  # transitions

    writer.write('\n')


def write_hmm(hmm, writer, strain, chrm, emis_symbols):
    writer.write(f'{strain}\t{chrm}\t')

    states = len(hmm.hidden_states)
    writer.write('\t'.join(
        [f'{p}' for p in hmm.initial_p] +  # initial
        [f'{hmm.emissions[i, hmm.symbol_to_ind[symbol]]}'
         if symbol in hmm.symbol_to_ind else '0.0'
         for i in range(states)
         for symbol in emis_symbols] +  # emission
        [f'{hmm.transitions[i, j]}'
         for i in range(states)
         for j in range(states)]  # transition
    ))
    writer.write('\n')


def write_state_probs(probs, writer, strain, chrm, states):
    writer.write(f'{strain}\t{chrm}\t')

    writer.write('\t'.join(
        [f'{state}:' +
         ','.join([f'{site[i]:.5f}' for site in probs])
         for i, state in enumerate(states)]))

    writer.write('\n')
