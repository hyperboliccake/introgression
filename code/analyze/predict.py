import copy
import gzip
import itertools
from collections import defaultdict
from hmm import hmm_bw
from sim import sim_predict
from sim import sim_process
import global_params as gp
from misc import read_fasta

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

def initial_probabilities(known_states, unknown_states, \
                          expected_frac, weighted_match_freqs):
    
    init = []
    expectation_weight = .9
    for s in range(len(known_states)):
        state = known_states[s]
        expected = expected_frac[state]
        estimated = weighted_match_freqs[s]
        init.append(expected * expectation_weight + \
                    estimated * (1 - expectation_weight))
    for s in range(len(unknown_states)):
        state = unknown_states[s]
        expected_frac = expected_frac[state]
        init.append(expected_frac)
        
    return norm_list(init)

def emission_probabilities(known_states, unknown_states, symbols):

    prob_P = .9 # matches current but not master
    prob_B = .09 # matches current and master
    prob_X = .009 # matches neither
    prob_C = .001 # matches master but not current
    
    emis = []

    # this is sort of hardcoding...fix it?
    num_per_category = 2 ** (len(known_states) - 2)
    for s in range(len(known_states)):
        state = known_states[s]
        emis.append(defaultdict(float))
        for symbol in symbols:
            match_master = symbol[0] == gp.match_symbol
            match_current = symbol[s] == gp.match_symbol
            if match_current:
                if match_master:
                    emis[s][symbol] = prob_B / num_per_category
                else:
                    emis[s][symbol] = prob_P / num_per_category
            else:
                if match_master:
                    emis[s][symbol] = prob_C / num_per_category
                else:
                    emis[s][symbol] = prob_X / num_per_category
                    
        emis[s] = norm_dict(emis[s])
    
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
                             expected_frac, expected_tract_lengths):

    # doesn't depend on sequence observations but maybe it should?

    # also should we care about number of tracts rather than fraction
    # of genome? maybe theoretically, but that number is a lot more
    # suspect

    # trans[i][j!=i] = 
    # prob(transition) * frac(j) = 
    # 1/expected tract length(i) * expected frac(j)|not i

    states = known_states + unknown_states

    trans = []
    for i in range(len(states)):
        state_from = states[i]
        trans.append([])
        scale_other = 1 / (1 - expected_frac[state_from])
        for j in range(len(states)):
            state_to = states[j]
            if state_from == state_to:
                trans[i].append(1 - 1./expected_tract_lengths[state_from])
            else:
                trans[i].append(1./expected_tract_lengths[state_from] * \
                                expected_frac[state_to] * scale_other)

        trans[i] = norm_list(trans[i])

    return trans

def initial_hmm_parameters(seq, known_states, unknown_states, \
                           expected_frac, expected_tract_lengths):

    # get frequencies of individual symbols (e.g. '+') and all full
    # combinations of symbols (e.g. '+++-')
    individual_symbol_freqs, symbol_freqs, weighted_match_freqs = get_symbol_freqs(seq)

    init = initial_probabilities(known_states, unknown_states, \
                                 expected_frac, weighted_match_freqs)
    emis = emission_probabilities(known_states, unknown_states, symbol_freqs.keys())
    trans = transition_probabilities(known_states, unknown_states, \
                                     expected_frac, expected_tract_lengths)

    return init, emis, trans

def predict_introgressed(ref_seqs, predict_seq, predict_args, \
                         train=True):

    # get rid of ++ sites?
    only_poly_sites = True

    method = 'posterior'
    if not type(predict_args['threshold']) is float:
        method = 'viterbi'

    # code sequence by which reference it matches at each site
    seq_coded, ref_seqs_coded, ps = ungap_and_code(predict_seq, ref_seqs)
    if only_poly_sites:
        seq_coded, ref_seqs_coded, ps = poly_sites(seq_coded, ref_seqs_coded, ps)

    # sets expected number of tracts and bases for each reference
    # based on expected length of introgressed tracts and expected
    # total fraction of genome
    set_expectations(predict_args, len(predict_seq))

    # set initial hmm parameters based on combination of (1) initial
    # expectations (length of introgressed tract and fraction of
    # genome/total number tracts and bases) and (2) number of sites at
    # which predict seq matches each reference
    init, emis, trans = initial_hmm_parameters(seq_coded, \
                                               predict_args['known_states'], \
                                               predict_args['unknown_states'], \
                                               predict_args['expected_frac'], \
                                               predict_args['expected_tract_lengths'])

    ######
    # make predictions
    ######

    default_state = predict_args['states'][0]

    # new Hidden Markov Model
    hmm = hmm_bw.HMM()

    # set states and initial probabilties
    hmm.set_states(predict_args['states'])
    hmm.set_initial_p(init)
    hmm.set_emissions(emis)
    hmm.set_transitions(trans)

    # set obs
    print(seq_coded)
    hmm.set_observations([seq_coded])

    hmm_init = copy.deepcopy(hmm)

    # optional Baum-Welch parameter estimation
    if train:
        hmm.go(predict_args['improvement_frac'])

    #if method == "posterior":
    predicted = {}
    all_probs = {}
    # for all obs sequences, each site is a dic with one prob for each
    # state
    p = hmm.posterior_decoding() # returns list but we're hackin this
                                 # for just one species right now
    path, path_probs = sim_process.get_max_path(p[0])

    if method == "posterior":
        path_t = sim_process.threshold_predicted(path, path_probs, \
                                                 predict_args['threshold'], \
                                                 default_state)
        return path_t, p[0], hmm, hmm_init, ps
        
    if method == 'viterbi':
        hmm.set_observations([seq_coded])
        predicted = sim_predict.convert_predictions(hmm.viterbi(), predict_args['states'])
        return predicted, p[0], hmm, hmm_init, ps

    else:
        print('invalid method')

def write_positions(ps, f, strain, chrm):
    sep = '\t'
    f.write(strain + sep + chrm + sep + sep.join([str(x) for x in ps]) + '\n')
    f.flush()

def read_positions(fn):
    # dictionary keyed by strain and then chromosome
    f = gzip.open(fn, 'rb')
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
    # NOTE: num_sites_hmm represents the sites considered by the HMM,
    # so it might exclude non-polymorphic sites in addition to gaps
    f.write('strain' + sep + 'chromosome' + sep + 'predicted_species' + sep + \
            'start' + sep + 'end' + sep + 'num_sites_hmm' + '\n')
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

def read_blocks(fn, labeled=False):
    # return dictionary of (start, end, number_non_gap, [region_id]),
    # keyed by strain and then chromosome
    f = open(fn, 'r')
    f.readline() # header
    line = f.readline()
    d = defaultdict(lambda: defaultdict(list))
    while line != '':
        if labeled: 
            region_id, strain, chrm, species, start, end, number_non_gap = \
                line.strip().split('\t')
            d[strain][chrm].append((region_id, int(start), int(end), \
                                    int(number_non_gap)))
        else:
            strain, chrm, species, start, end, number_non_gap = line.strip().split('\t')
            d[strain][chrm].append((int(start), int(end), int(number_non_gap)))
        line = f.readline()
    f.close()
    return d

def get_emis_symbols(known_states): 

    symbols = [gp.match_symbol, gp.mismatch_symbol]
    emis_symbols = [''.join(x) for x in \
                           list(itertools.product(symbols, repeat=len(known_states)))]
    emis_symbols.sort()
    return emis_symbols

def write_hmm_header(known_states, unknown_states, emis_symbols, f):

    sep = '\t'
    header_string = 'strain' + sep + 'chromosome' + sep

    states = known_states + unknown_states

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

    return emis_symbols

def write_hmm(hmm, f, strain, chrm, emis_symbols):
    
    sep = '\t'

    line_string = strain + sep + chrm + sep

    # initial
    for i in range(len(hmm.states)):
        line_string += str(hmm.initial_p[i]) + sep

    # emission
    for i in range(len(hmm.states)):
        for symbol in emis_symbols:
            line_string += str(hmm.emissions[i, hmm.symbol_to_ind[symbol]]) + sep
    # transition
    for i in range(len(hmm.states)):
        for j in range(len(hmm.states)):
            line_string += str(hmm.transitions[i][j]) + sep

    f.write(line_string[:-(len(sep))] + '\n')
    f.flush()

def write_state_probs(probs, f, strain, chrm):

    # probs is list of sites, each site dic keyed by state
    
    # file format is:
    # strain\tchrm\tcer:.1,.2,.3\tpar:.9,.8,.7
    # strain\tchrm\tcer:.4,.2,.3\tpar:.6,.8,.7

    sep  ='\t'
    f.write(strain + sep + chrm)

    for state in probs[0].keys():
        f.write(sep + state + ':')
        probs_string = ','.join(["{0:.5f}".format(site[state]) \
                                     for site in probs])
        f.write(probs_string)
    f.write('\n')
    f.flush()
