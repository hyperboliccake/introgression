import sys
import os
import copy
import itertools
sys.path.append('..')
import global_params as gp
sys.path.append('../hmm')
import hmm_bw

# add in the nonpolymorphic sites
def fill_seqs(polymorphic_seqs, polymorphic_sites, nsites, fill):
    
    seqs_filled = []
    for seq in polymorphic_seqs:
        s = ''
        poly_ind = 0
        for i in range(nsites):
            if i in polymorphic_sites:
                s += seq[poly_ind]
                poly_ind += 1
            else:
                s += fill
        seqs_filled.append(s)
    return seqs_filled

# convert sequences from bases to symbols indicating which references they match
def code_seqs(seqs, nsites, ref_seqs):

    nrefs = len(ref_seqs)
    seqs_coded = []
    for seq in seqs:
        s = []
        for i in range(nsites):
            si = ''
            for r in range(nrefs):
                if ref_seqs[r][i] == gp.unsequenced_symbol:
                    si += gp.unknown_symbol
                elif ref_seqs[r][i] == seq[i]:
                    si += gp.match_symbol
                else:
                    si += gp.mismatch_symbol
            s.append(si)
        seqs_coded.append(s)
    return seqs_coded

def get_symbol_freqs_one(seqs, args):

    if len(seqs) == 0:
        return None

    known_states = copy.deepcopy(args['states'])
    if args['unknown_species'] != None:
        known_states.remove(args['unknown_species'])

    # for a single species

    known_symbols = [gp.match_symbol, gp.mismatch_symbol]
    symbols = known_symbols + [gp.unknown_symbol]

    # for species a(0): +, -; species b(1): +, -...etc
    individual_symbol_freqs = []
    for i in range(len(known_states)):
        d_species = {}
        for symbol in known_symbols:
            num = 0
            den = 0
            for seq in seqs:
                seq_species = [x[i] for x in seq]
                num += seq_species.count(symbol)
                den += len(seq_species) - seq_species.count(gp.unknown_symbol)
            d_species[symbol] = float(num) / den
        individual_symbol_freqs.append(d_species)

    # for +++, ++-, +-+, etc
    symbol_combinations = [''.join(x) for x in \
                           list(itertools.product(symbols, repeat=len(known_states)))]
    symbol_freqs = dict(zip(symbol_combinations, [0]*len(symbol_combinations)))
    for symbol in symbol_freqs.keys():
        num = 0
        den = 0
        for seq in seqs:
            num += seq.count(symbol)
            den += len(seq)
        symbol_freqs[symbol] = float(num) / den

    # weighted matches; for first species +-- gets 1 pt, +-+ and ++-
    # each get 1/2 pt, +++ gets 1/3 pt [treat ? like -] etc
    weighted_match_freqs = {}
    total = 0
    for i in range(len(known_states)):
        current = 0
        keep_symbols = [] # for first species +**
        weights = []
        for s in symbol_combinations:
            if s[i] == gp.match_symbol:
                keep_symbols.append(s)
                weights.append(1./s.count(gp.match_symbol))
        for seq in seqs:
            for j in range(len(keep_symbols)):
                current += seq.count(keep_symbols[j]) * weights[j]
        weighted_match_freqs[known_states[i]] = current
        total += current

    # for unknown species 1 pt for --- only (or ?)
    if args['unknown_species'] != None:
        current = 0
        keep_symbols = []
        for s in symbol_combinations:
            if gp.match_symbol not in s:
                keep_symbols.append(s)
        for seq in seqs:
            for j in range(len(keep_symbols)):
                current += seq.count(keep_symbols[j])
        weighted_match_freqs[args['unknown_species']] = current
        total += current

    for state in weighted_match_freqs:
        # total can only be zero if no species has matches
        weighted_match_freqs[state] /= float(total)
    
    return individual_symbol_freqs, symbol_freqs, weighted_match_freqs

def get_symbol_freqs(seqs_coded, args):

    # for each species/state set of sequences calculate: the
    # frequencies of individual symbols for each species, the
    # frequencies of combined symbols, and weighted match frequencies

    d = {}
    for state in args['states']:
        seqs_current = []
        for i in args['species_to_indices'][state]:
            if i not in args['ref_inds']:
                seqs_current.append(seqs_coded[i])
        d[state] = get_symbol_freqs_one(seqs_current, args)

    return d

def initial_probabilities(weighted_match_freqs, args):

    # a small value to add so that no frequencies are actually 0
    epsilon = .001

    # initial probabilities are proportional to number of match
    # symbols for the state, but also factor in how often we expect to
    # be in each state
    init = []
    for state in args['states']:
        # add these two things because we want to average them instead
        # of letting one of them bring the total to zero [should we
        # really give them equal weight though?]
        init.append(weighted_match_freqs[state] + \
                        float(args['expected_tract_lengths'][state] * \
                              args['expected_num_tracts'][state]) / \
                    args['num_sites'] + epsilon)
    scale = float(sum(init))
    for i in range(len(args['states'])):
        init[i] /= scale
        assert init[i] > 0, 'my HMM methods break down with zero probabilities (just get rid of the state if you don\'t want it!)'
    return init

# unlike for other parameters, calculate probabilities from the
# appropriate species sequences instead of just the one being
# predicted
def emission_probabilities(d_freqs, own_bias, args):

    # a small value to add so that no frequencies are actually 0
    epsilon = .001

    # idea is to assign emission probabilities dependant on frequency
    # of symbols; so if we want to get prob of par emitting ++- (cer
    # par bay), then we take frequency of ++-, multiplied by own_bias
    # (because + for par); if it's +?- then we have to do some more
    # complicated guessing; also note the frequencies we're looking at
    # are just from the sequences for the current species, not the ones
    # we're trying to predict
    emis = []
    for i in range(len(args['states'])):
        emis_species = {}
        individual_symbol_freqs, symbol_freqs, weighted_match_freqs \
            = d_freqs[args['states'][i]]
        for symbol in symbol_freqs.keys():
            p = symbol_freqs[symbol] + epsilon
            if args['states'][i] == args['unknown_species']:
                # treat ? as - for now
                if gp.match_symbol in symbol:
                    p *= (1 - own_bias)
                else:
                    p *= own_bias
            elif symbol[i] == gp.match_symbol:
                p *= own_bias
            else: #symbol[i] == mismatch_symbol:
                p *= (1 - own_bias)
            # maybe this kind of guessing for '?' is worth doing? but
            # should probably also do it for the initial/transition
            # probabilities as well if doing it here
            #else:
            #    p *= ((own_bias) * individual_symbol_freqs[i][match_symbol] + \
            #              (1 - own_bias) * individual_symbol_freqs[i][mismatch_symbol])
            emis_species[symbol] = p
        norm = float(sum(emis_species.values()))
        for symbol in emis_species:
            emis_species[symbol] /= norm
            assert emis_species[symbol] > 0, 'my HMM methods break down with zero probabilities (just get rid of the state if you don\'t want it!)'
        emis.append(emis_species)
    return emis

def transition_probabilities(weighted_match_freqs, args):

    # a small value to add so that no frequencies are actually 0
    epsilon = 1./1000000

    ##states_not_predict = states
    ##states_not_predict.remove(predict_species)
    
    # if a is the species with introgression in it
    
    # transition a->b 1/length not int a * frac to b
    # transition a->c 1/length not int a * frac to c

    # frac to b - num match only b / num match only b or match only c
    # frac to c - num match only c / num match only b or match only c

    # b->a 1/length int b
    # b->c basically 0

    # frac to a - 
    # frac to c - num int c / num int c + num not int

    # num int c = bases int c / length int c

    # bases int c = [(+c-b) + (+c+b) * 1/2] * some fraction accounting for randomness

    # c->a 1/length int c
    # c->b basically 0

    # TODO should be able to calculate expected length based on expected time?

    # TODO should be able to calculate expected length based on
    # expected amount of migration?

    expected_length_not_introgressed = \
        float(args['expected_tract_lengths'][args['species_to']])

    # fraction of time we should choose given species state over
    # others based on number of sites that match it
    fracs = {}
    frac_den = 0
    for state in args['states']:
        fracs[state] = weighted_match_freqs[state]
        frac_den += fracs[state]
    # TODO: question: if A:10 and B:12 and all of those matches
    # overlap, then should B get it 12/22 of the time or all the time?
    # i.e. is it just the sites that uniquely match that are relevant?
    # i think yes, the unique ones BUT want to make sure we don't set
    # any probabilities to 0
    for state in fracs:
        fracs[state] /= float(frac_den)

    trans = []
    for state_from in args['states']:
        trans_current = {}
        total = 0

        for state_to in args['states']:
            val = 0
            # staying within same species, deal with this later by
            # figuring out what's left
            if state_to == state_from:
                pass
            # moving from non-introgressed (cer) to introgressed
            elif state_from == args['species_to']:
                if expected_length_not_introgressed > 0:
                    val = 1 / expected_length_not_introgressed * fracs[state_to]
                else:
                    # we should definitely transition if we expect
                    # entire sequence to be introgressed
                    val = 1
            # moving from introgressed to non-introgressed
            elif state_to == args['species_to']:
                if args['expected_tract_lengths'][state_from] > 0:
                    val = 1 / float(args['expected_tract_lengths'][state_from])
                else:
                    # we should definitely transition if we don't
                    # expect any introgression
                    val = 1
            else:
                # from one state to another, where neither is the one
                # we're trying to predict TODO
                val = 0
            total += val
            trans_current[state_to] = val

        trans_current[state_from] = 1 - total

        # add epsilon and normalize to avoid 0 probabilities
        total = 1 + len(args['states']) * epsilon
        for state_to in args['states']:
            trans_current[state_to] += epsilon
            trans_current[state_to] /= total

        assert trans_current[state_from] > 0, 'my HMM methods break down with zero probabilities (just get rid of the state if you don\'t want it!)'
        trans_row = []
        for state_to in args['states']:
            assert trans_current[state_to] > 0, 'my HMM methods break down with zero probabilities (just get rid of the state if you don\'t want it!)'
            trans_row.append(trans_current[state_to])
        trans.append(trans_row)

    return trans

def initial_hmm_parameters(seqs_coded, args):

    # get frequencies of all symbols (i.e. matching to each
    # reference/all combinations of references)
    d_freqs = get_symbol_freqs(seqs_coded, args)

    # using only the sequences to predict introgression in: (1) the
    # frequency of alignment columns that match/don't match each
    # reference, (2) the frequency of all symbol combinations, (3)
    # weights for each symbol combination/reference
    individual_symbol_freqs, symbol_freqs, weighted_match_freqs = \
        d_freqs[args['species_to']]

    p = {}
    p['init'] = initial_probabilities(weighted_match_freqs, args)
    p['emis'] = emission_probabilities(d_freqs, .99, args)
    p['trans'] = transition_probabilities(weighted_match_freqs, args)

    return p['init'], p['emis'], p['trans']

def convert_predictions(path, states):
    new_path = []
    for p in path:
        new_path.append(states[p])
    return new_path

def run_hmm(seqs, args, init, emis, trans, train):

    # sanity checks
    if args['unknown_species'] != None:
        assert len(seqs[0][0]) == len(args['states']) - 1
        assert set(args['index_to_species']) == set(args['states']), \
            str(set(args['index_to_species'])) + ' ' + str(args['states'])
        assert args['unknown_species'] == args['states'][-1]

    # only make predictions for sequences that are species_to (and
    # that are not the reference sequences, which we clearly can't
    # make predictions for)
    predict_inds = copy.deepcopy(args['species_to_indices'][args['species_to']])
    predict_inds.remove(args['ref_inds'][0])
    seqs_to_predict = [seqs[x] for x in predict_inds]

    # new Hidden Markov Model
    hmm = hmm_bw.HMM()

    # set obs
    hmm.set_obs(seqs_to_predict)

    # set states and initial probabilties
    hmm.set_states(args['states'])
    hmm.set_init(init)
    hmm.set_emis(emis)
    hmm.set_trans(trans)

    hmm_init = copy.deepcopy(hmm)

    # optional Baum-Welch parameter estimation
    if train:
        hmm.go()

    # make predictions!
    predicted = {}
    for i in predict_inds:
        hmm.set_obs(seqs[i])
        predicted[i] = convert_predictions(hmm.viterbi(), args['states'])

    return predicted, hmm, hmm_init

def set_up_seqs(sim, args): 

    # fill in nonpolymorphic sites
    fill_symbol = '0'
    seqs_filled = fill_seqs(sim['seqs'], sim['positions'], \
                            args['num_sites'], fill_symbol)

    # code sequences by which references they match at each position
    ref_seqs = [seqs_filled[r] for r in args['ref_inds']]
    seqs_coded = code_seqs(seqs_filled, args['num_sites'], ref_seqs)

    return seqs_coded

def predict_introgressed(sim, args, train):
    
    seqs_coded = set_up_seqs(sim, args)

    # initial values for initial, emission, and transition
    # probabilities
    init, emis, trans = initial_hmm_parameters(seqs_coded, args)

    # make predictions
    predicted, hmm, hmm_init = run_hmm(seqs_coded, args, init, emis, trans, train)

    return predicted, hmm, hmm_init

def write_hmm_headers(states, emis_symbols, f, sep):

    header_string = ''

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

def write_hmm_line(hmm, f, header = False):

    sep = '\t'
    emis_symbols = hmm.emis[0].keys()

    if header:
        write_hmm_headers(hmm.states, emis_symbols, f, sep)

    line_string = ''

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
