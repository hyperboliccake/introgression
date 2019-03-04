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


def process_predict_args(arg_list):

    d = {}
    i = 0

    d['tag'] = arg_list[i]
    i += 1

    d['improvement_frac'] = float(arg_list[i])
    i += 1

    d['threshold'] = 'viterbi'
    try:
        d['threshold'] = float(arg_list[i])
    except:
        pass
    i += 1

    # expected length of introgressed tracts and fraction of sequence
    # introgressed
    expected_tract_lengths = {}
    expected_frac = {}

    d['known_states'] = gp.alignment_ref_order
    for ref in gp.alignment_ref_order[1:]:
        expected_tract_lengths[ref] = float(arg_list[i])
        i += 1
        expected_frac[ref] = float(arg_list[i])
        i += 1

    d['unknown_states'] = []
    while i < len(arg_list):
        state = arg_list[i]
        d['unknown_states'].append(state)
        i += 1
        expected_tract_lengths[state] = float(arg_list[i])
        i += 1
        expected_frac[state] = float(arg_list[i])
        i += 1

    d['states'] = d['known_states'] + d['unknown_states']

    expected_frac[d['states'][0]] = 0
    expected_frac[d['states'][0]] = 1 - sum(expected_frac.values())
    d['expected_frac'] = expected_frac

    # calculate these based on remaining bases, but after we know
    # which chromosome we're looking at
    expected_tract_lengths[d['states'][0]] = 0    
    d['expected_tract_lengths'] = expected_tract_lengths
    d['expected_num_tracts'] = {}
    d['expected_bases'] = {}

    return d

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
        args['expected_bases'][species_to] /\
        args['expected_num_tracts'][species_to]


def ungap_and_code_helper(predict_seq, ref_seqs, index_ref):

    # index_ref is index of reference strain to index relative to

    ps = []  # positions that we're keeping in analysis
    seq = []  # sequence of emitted symbols, e.g. '++-++'
    ind = 0  # current position in reference strain

    # build character array
    seqs = np.array([list(predict_seq)] +
                    [list(r) for r in ref_seqs])

    # make boolean for valid characters
    isvalid = np.logical_and(seqs != gp.gap_symbol,
                             seqs != gp.unsequenced_symbol)

    # positions are where everything is valid, removing invalid reference first
    # +1 removes the predict sequence at index 0
    ps = np.where(np.all(isvalid[:, isvalid[index_ref+1, :]], axis=0))[0]

    matches = np.apply_along_axis(
        lambda s: np.where(s == seqs, gp.match_symbol, gp.mismatch_symbol),
        axis=1,
        arr=seqs
    )

    #print((matches[0, :, np.all(isvalid, axis=0)]))
    # 1: indexing removes currently examined sequence
    # ref sequences not used!? (commit then test)
    seq = [''.join(row) for row in matches[0, 1:, np.all(isvalid, axis=0)]]
    return seq, list(ps)
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
            #seq.append(symbol)
            ps.append(ind)
        if ref_seqs[index_ref][i] != gp.gap_symbol:
            ind += 1

    print('---')
    print(ps)
    print(index_ref)
    print(isvalid)
    print(ps)
    # NOTE list is for unit test comparisons


def ungap_and_code(predict_seq, ref_seqs, index_ref=0):
    seq, ps = ungap_and_code_helper(predict_seq, ref_seqs, index_ref)
    ref_seqs_coded = []
    for r in ref_seqs:
        ref_seq_coded, _ = ungap_and_code_helper(r, ref_seqs, index_ref)
        ref_seqs_coded.append(ref_seq_coded)
    return seq, ref_seqs_coded, ps


def poly_sites(seq, ref_seqs, ps):
    seq_len = len(seq[0])
    # check if seq only contains match_symbol
    retain = np.vectorize(lambda x: x.count(gp.match_symbol) != seq_len)(seq)
    indices = np.where(retain)[0]
    ps_poly = [ps[i] for i in indices]
    seq_poly = [seq[i] for i in indices]
    ref_seqs_poly = [[r[i] for i in indices] for r in ref_seqs]
    return seq_poly, ref_seqs_poly, ps_poly


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


def initial_hmm_parameters(seq, known_states, unknown_states,
                           expected_frac, expected_tract_lengths):

    # get frequencies of individual symbols (e.g. '+') and all full
    # combinations of symbols (e.g. '+++-')
    individual_symbol_freqs, symbol_freqs, weighted_match_freqs = get_symbol_freqs(seq)

    init = initial_probabilities(known_states, unknown_states,
                                 expected_frac, weighted_match_freqs)
    emis = emission_probabilities(known_states, unknown_states, symbol_freqs.keys())
    trans = transition_probabilities(known_states, unknown_states,
                                     expected_frac, expected_tract_lengths)

    # new Hidden Markov Model
    hmm = hmm_bw.HMM()

    hmm.set_initial_p(init)
    hmm.set_emissions(emis)
    hmm.set_transitions(trans)
    return hmm


def predict_introgressed(ref_seqs, predict_seq, predict_args,
                         train=True, only_poly_sites=True):

    # code sequence by which reference it matches at each site
    seq_coded, ref_seqs_coded, ps = ungap_and_code(predict_seq, ref_seqs)
    if only_poly_sites:
        seq_coded, ref_seqs_coded, ps = poly_sites(seq_coded,
                                                   ref_seqs_coded, ps)

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
                                 predict_args['expected_tract_lengths'])

    ######
    # make predictions
    ######

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
        return path_t, p[0], hmm, hmm_init, ps

    else:
        hmm.set_observations([seq_coded])
        predicted = sim_predict.convert_predictions(hmm.viterbi(),
                                                    predict_args['states'])
        return predicted, p[0], hmm, hmm_init, ps


def write_positions(ps, writer, strain, chrm):
    writer.write(f'{strain}\t{chrm}\t' + '\t'.join([str(x) for x in ps]) + '\n')


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
    writer.write('\n')  # for compatibility


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
