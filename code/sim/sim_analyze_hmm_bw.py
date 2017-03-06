import os
import sys
import copy
import itertools
from concordance_functions import *

sys.path.append('../hmm')
from hmm_bw import *

sys.path.append('..')
import global_params as gp


theory = False
theory_done = True

# to use instead of d[key] = value, so we can ensure we're not
# accidentally creating new keys
def update_value(d, key, value):
    assert key in d, key + ' not in ' + str(d.keys())
    d[key] = value
    return d
        

def mean(a):
    try: 
        return sum(a) / float(len(a))
    except:
        return -1

# calculate percentage of sites at which sequences a and b match
def match_percent(a, b, n):
    assert(len(a) == len(b))
    match = 0
    for i in xrange(len(a)):
        if a[i] == b[i]:
            match += 1
    return (float(match + n - len(a))) / n

def process_args(l):

    i = 1

    tag = sys.argv[i]
    i += 1

    topology = sys.argv[i]
    i += 1

    species = get_labels(parse_topology(topology))
    assert len(species) == 2 or len(species) == 3, species

    expected_tract_lengths = {}
    expected_num_tracts = {}

    species_to = sys.argv[i]
    i += 1
    assert species_to in species
    num_samples_species_to = int(sys.argv[i])
    i += 1
    N0_species_to = int(sys.argv[i])
    i += 1
    assert sys.argv[i] == 'to'
    i += 1

    species_from1 = sys.argv[i]
    i += 1
    assert species_from1 in species
    num_samples_species_from1 = int(sys.argv[i])
    i += 1
    N0_species_from1 = int(sys.argv[i])
    i += 1
    migration_from1 = float(sys.argv[i]) * 2 * N0_species_from1
    i += 1
    print 'migration', migration_from1
    expected_tract_lengths[species_from1] = float(sys.argv[i])
    i += 1
    expected_num_tracts[species_from1] = int(sys.argv[i])
    i += 1
    has_ref_from1 = (sys.argv[i] == 'ref')
    i += 1

    species_from2 = None
    num_samples_species_from2 = 0
    N0_species_from2 = N0_species_from1
    migration_from2 = 0
    has_ref_from2 = True
    if len(species) == 3:
        species_from2 = sys.argv[i]
        i += 1
        assert species_from2 in species
        num_samples_species_from2 = int(sys.argv[i])
        i += 1
        N0_species_from2 = int(sys.argv[i])
        i += 1
        migration_from2 = float(sys.argv[i]) * 2 * N0_species_from2
        i += 1
        expected_tract_lengths[species_from2] = float(sys.argv[i])
        i += 1
        expected_num_tracts[species_from2] = int(sys.argv[i])
        i += 1
        has_ref_from2 = (sys.argv[i] == 'ref')
        i += 1

    # only makes sense to have one unknown species at most
    assert has_ref_from1 or has_ref_from2

    # calculate these based on remaining bases
    expected_num_tracts[species_to] = sum(expected_num_tracts.values()) + 1
    expected_num_introgressed_bases = expected_tract_lengths[species_from1] * \
        expected_num_tracts[species_from1]
    if species_from2 != None:
        expected_num_introgressed_bases += expected_tract_lengths[species_from2] * \
            expected_num_tracts[species_from2]
    expected_tract_lengths[species_to] = float(expected_num_introgressed_bases) / \
        expected_num_tracts[species_to]

    assert N0_species_to == N0_species_from1 and N0_species_to == N0_species_from2

    topology = parse_topology(topology, 1/float(2 * N0_species_to))

    # 13,500 sites to get about 10% with one recombination event, .3% with
    # more than one (based on poisson(.1), 1 recombination per chromosome
    # of average length 750,000)
    num_sites = int(sys.argv[i])
    i += 1

    # parameter is recombination rate between adjacent bp per generation
    # should probably be 1/750000 + 6.1 * 10^-6 (where 750000 is average
    # chr size)
    print 'recomb rate', float(sys.argv[i])
    rho = 2 * N0_species_to * float(sys.argv[i]) * (num_sites - 1)
    i += 1

    outcross_rate = float(sys.argv[i])
    i += 1
    
    print 'outcross rate', outcross_rate
    rho *= outcross_rate
    print 'rho', rho

    theta = gp.mu * 2 * num_sites * N0_species_to
    print 'theta', theta

    num_reps = int(sys.argv[i])
    
    return tag, topology, species_to, species_from1, species_from2, \
        num_samples_species_to, num_samples_species_from1, num_samples_species_from2, \
        N0_species_to, N0_species_from1, N0_species_from2, \
        migration_from1, migration_from2, \
        expected_tract_lengths, \
        expected_num_tracts, \
        has_ref_from1, has_ref_from2, \
        rho, outcross_rate, theta, num_sites, num_reps

def process_args_old(l):

    i = 1

    tag = sys.argv[i]
    i += 1

    model = sys.argv[i]
    i += 1

    N0 = int(sys.argv[i])
    i += 1

    include_bay = (sys.argv[i] == 'include_bay')
    i += 1
    include_unk = (sys.argv[i] == 'include_unk')
    i += 1

    num_samples_par = int(sys.argv[i])
    i += 1
    num_samples_cer = int(sys.argv[i])
    i += 1
    num_samples_bay = 0
    if include_bay:
        num_samples_bay = int(sys.argv[i])
        i += 1
        
    num_samples = num_samples_par + num_samples_cer + num_samples_bay

    # migration parameter is 2 * N0 * m, where mij is fraction of i made
    # up of j each generation; need to figure out how to make migration
    # rates equivalent for different models
    par_cer_migration = 2 * N0 * float(sys.argv[i])
    i += 1

    bay_cer_migration = 2 * N0 * float(sys.argv[i])
    i += 1

    # in generations
    t_par_cer = float(sys.argv[i]) / (2 * N0)
    i += 1
    t_bay_par_cer = float(sys.argv[i]) / (2 * N0)
    i += 1

    # 13,500 sites to get about 10% with one recombination event, .3% with
    # more than one (based on poisson(.1), 1 recombination per chromosome
    # of average length 750,000)
    num_sites = int(sys.argv[i])
    i += 1

    # parameter is recombination rate between adjacent bp per generation
    # should probably be 1/750000 + 6.1 * 10^-6 (where 750000 is average
    # chr size)
    rho = 2 * N0 * float(sys.argv[i]) * (num_sites - 1)
    i += 1

    outcross_rate = float(sys.argv[i])
    i += 1
    
    rho *= outcross_rate

    # estimate from humans
    mu = 1.84 * 10 ** -10
    theta = mu * 2 * num_sites * N0

    num_reps = int(sys.argv[i])
    
    return tag, model, N0, include_bay, include_unk,\
        num_samples_cer, num_samples_par, num_samples_bay,\
        par_cer_migration, bay_cer_migration,\
        t_par_cer, t_bay_par_cer,\
        num_sites, rho, theta, outcross_rate, num_reps

def parse_ms_tree_helper(t):
    if '(' not in t:
        colon_ind = t.find(':')
        # index from 0 instead of 1
        return [int(t[:colon_ind]) - 1, float(t[colon_ind+1:])]
        #return [int(t[:colon_ind]) - 1, Decimal(t[colon_ind+1:])]
    left = ''
    right = ''
    i = -1
    if t[1] != '(':
        comma_ind = t.find(',') 
        left = t[1:comma_ind]
        right = t[comma_ind+1:t.rfind(')')]
    else:
        open_count = 1
        i = 2
        while open_count > 0:
            if t[i] == '(':
                open_count += 1
            elif t[i] == ')':
                open_count -= 1
            i += 1
        i = t[i:].find(',')+i
        left = t[1:i]
        right = t[i+1:t.rfind(')')]

    time = float(t[t.rfind(':')+1:])
    #time = Decimal(t[t.rfind(':')+1:])
    return [parse_ms_tree_helper(left), parse_ms_tree_helper(right), time]

# converts newick tree string to nested list format (assuming labels from ms output)
def parse_ms_tree(t):
    return parse_ms_tree_helper(t[:-1] + ':0')

def parse_topology_helper(t, factor=1):
    if '(' not in t:
        return t
    left = ''
    right = ''
    i = -1
    # if left subtree is just a label, easier to split into left and
    # right
    if t[1] != '(':
        comma_ind = t.find(',') 
        left = t[1:comma_ind]
        right = t[comma_ind+1:t.rfind(')')]
    # otherwise, have to figure out where left subtree ends by
    # matching parens
    else:
        open_count = 1
        i = 2
        while open_count > 0:
            if t[i] == '(':
                open_count += 1
            elif t[i] == ')':
                open_count -= 1
            i += 1
        i = t[i:].find(',')+i
        left = t[1:i]
        right = t[i+1:t.rfind(')')]

    time = float(t[t.rfind(':')+1:]) * factor
    return [parse_topology_helper(left, factor), \
                parse_topology_helper(right, factor), \
                time]

# parses a topology that also includes divergence times
# e.g. '(cer,par):375000000,bay):1125000000'
# [this is not the normal way of writing branch lengths in newick trees]
def parse_topology(t, factor=1):
    return parse_topology_helper(t, factor)

# actual introgressed with 2 species total, within a single unrecombined block
def find_introgressed_2(t, cutoff_time, species_to, index_to_species):
    # the subtrees that exist at the time the populations join
    lineages = split(t, cutoff_time)
    introgressed = copy.deepcopy(index_to_species)
    num_lineages_species_to = 0
    for l in lineages:
        not_introgressed, s = is_partial_clade(l, species_to, index_to_species)
        # if this lineage has only species_to members, then add this
        # to the number of lineages for the species
        if not_introgressed:
            num_lineages_species_to += 1
        # if there's any occurence of species_from in this clade, then
        # mark all individuals as coming from species_from (since
        # we're only allowing migration in one direction)
        else:
            for label in get_labels(l):
                introgressed[label] = s

    # only return the results for the to species
    introgressed_species_to = []
    for i in range(len(index_to_species)):
        if index_to_species[i] == species_to:
            introgressed_species_to.append(introgressed[i])
        else:
            assert introgressed[i] == index_to_species[i], \
                str(introgressed[i]) + ' '  + str(index_to_species[i])

    return introgressed_species_to, num_lineages_species_to

# actual introgressed with 3 species total, within a single unrecombined block
# THIS FUNCTION ASSUMES MIGRATION ONLY HAPPENS AFTER MOST RECENT DIVERGENCE
def find_introgressed_3(t, species_to, topology, index_to_species):
    # strategy is to look at lineages that exist at most recent join
    # time; if any are not all one species, then mark all members as
    # the species that's not the to species (there's no migration
    # between the two from species); then look at the more distant
    # join time; now we can't just check if the lineages are all one
    # species because two of them have joined; so instead...figure out
    # which species is the last to join (from the topology), then...
    # nvm, only allowing migration later on makes this less complicated
    
    join_time_species_to = topology[2]
    subtree = topology[0]
    last_to_join = topology[1]
    if type(subtree) != type([]):
        subtree = topology[1]
        last_to_join = topology[0]
    if last_to_join != species_to:
        join_time_species_to = subtree[2]

    return find_introgressed_2(t, join_time_species_to, species_to, index_to_species)

def initial_probabilities(states, weighted_match_freqs, match_symbol, \
                              expected_tract_lengths, expected_num_tracts, n):
    # initial probabilities are proportional to number of match
    # symbols for the state, but also factor in how often we expect to
    # be in each state
    init = []
    for state in states:
        # add these two things because we want to average them instead
        # of letting one of them bring the total to zero [should we
        # really give them equal weight though?]
        init.append(weighted_match_freqs[state] + \
            float(expected_tract_lengths[state] * expected_num_tracts[state]) / n)
    scale = float(sum(init))
    for i in range(len(states)):
        init[i] /= scale
    return init

# unlike for other parameters, calculate probabilities from the
# appropriate species sequences instead of just the one being
# predicted
def emission_probabilities(states, unknown_species, d_freqs, match_symbol, mismatch_symbol, own_bias):
    # idea is to assign emission probabilities dependant on frequency
    # of symbols; so if we want to get prob of par emitting ++- (cer
    # par bay), then we take frequency of ++-, multiplied by own_bias
    # (because + for par); if it's +?- then we have to do some more
    # complicated guessing; also note the frequencies we're looking at
    # are just from the sequences for the current species, not the ones
    # we're trying to predict
    emis = []
    for i in range(len(states)):
        emis_species = {}
        individual_symbol_freqs, symbol_freqs, weighted_match_freqs = d_freqs[states[i]]
        for symbol in symbol_freqs.keys():
            p = symbol_freqs[symbol]
            if states[i] == unknown_species:
                # treat ? as - for now
                if match_symbol in symbol:
                    p *= (1 - own_bias)
                else:
                    p *= own_bias
            elif symbol[i] == match_symbol:
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
        emis.append(emis_species)
    return emis

def transition_probalities(expected_length_introgressed, \
                               expected_num_introgressed_tracts, \
                               predict_species, num_bases, \
                               weighted_match_freqs):

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

    assert expected_length_introgressed >= 0, expected_length_introgressed
    assert expected_num_introgressed_tracts >= 0, expected_num_introgressed_tracts

    states = expected_length_introgressed.keys()

    # TODO should be able to calculate expected length based on expected time?
    #expected_length_introgressed = {}
    expected_num_introgressed_bases = {}
    for species in states:
        expected_num_introgressed_bases[species] = \
            expected_num_introgressed_tracts[species] * \
            expected_length_introgressed[species]
    # TODO should be able to calculate expected length based on expected amount of migration?
    #expected_num_introgressed_tracts = {}

    expected_num_not_introgressed_tracts = \
        sum(expected_num_introgressed_tracts.values()) + 1
    expected_num_not_introgressed_bases = \
        num_bases - sum(expected_num_introgressed_bases.values())
    expected_length_not_introgressed = float(expected_num_not_introgressed_bases) / \
        expected_num_not_introgressed_tracts
    
    assert expected_num_not_introgressed_tracts >= 0, \
        expected_num_not_introgressed_tracts
    assert expected_num_not_introgressed_bases >= 0, \
        expected_num_not_introgressed_bases
    assert expected_length_not_introgressed >= 0, \
        expected_length_not_introgressed

    # fraction of time we should choose given species state over
    # others based on number of sites that match it
    fracs = {}
    frac_den = 0
    for state in states:
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
    for state_from in states:
        trans_current = {}
        total = 0

        for state_to in states:
            val = 0
            if state_to == state_from:
                pass
            elif state_from == predict_species:
                val = 1 / float(expected_length_not_introgressed) * fracs[state_to]
            elif state_to == predict_species:
                val = 1 / float(expected_length_introgressed[state_from])
            else:
                # from one state to another, where neither is the one
                # we're trying to predict TODO
                val = 0
            trans_current[state_to] = val
            total += val

        trans_current[state_from] = 1 - total
        assert trans_current[state_from] >= 0
        trans_row = []
        for state_to in states:
            assert trans_current[state_to] >= 0
            trans_row.append(trans_current[state_to])
        trans.append(trans_row)

    return trans

def get_symbol_freqs(states, index_to_species, seqs, unknown_species, \
                         match_symbol, mismatch_symbol, unknown_symbol):

    # for each species/state set of sequences calculate: the
    # frequencies of individuals symbols for each species, the
    # frequencies of combined symbols, and weighted match frequencies

    d = {}
    for state in states:
        seqs_current = []
        for i in range(len(seqs)):
            if index_to_species[i] == state:
                seqs_current.append(seqs[i])
        d[state] = get_symbol_freqs_one(states, seqs_current, unknown_species, \
                                            match_symbol, mismatch_symbol, \
                                            unknown_symbol)
    return d

def get_symbol_freqs_one(states, seqs, unknown_species, \
                             match_symbol, mismatch_symbol, unknown_symbol):

    known_states = copy.deepcopy(states)
    if unknown_species != None:
        known_states.remove(unknown_species)

    # for a single species

    known_symbols = [match_symbol, mismatch_symbol]
    symbols = known_symbols + [unknown_symbol]

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
                den += len(seq_species) - seq_species.count(unknown_symbol)
            d_species[symbol] = float(num) / den
        individual_symbol_freqs.append(d_species)

    # for +++, ++-, +-+, etc
    symbol_combinations = [''.join(x) for x in list(itertools.product(symbols, repeat=len(known_states)))]
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
            if s[i] == match_symbol:
                keep_symbols.append(s)
                weights.append(1./s.count(match_symbol))
        for seq in seqs:
            for j in range(len(keep_symbols)):
                current += seq.count(keep_symbols[j]) * weights[j]
        weighted_match_freqs[known_states[i]] = current
        total += current

    # for unknown species 1 pt for --- only (or ?)
    current = 0
    keep_symbols = []
    for s in symbol_combinations:
        if match_symbol not in s:
            keep_symbols.append(s)
    for seq in seqs:
        for j in range(len(keep_symbols)):
            current += seq.count(keep_symbols[j])
    weighted_match_freqs[unknown_species] = current
    total += current
    
    for state in weighted_match_freqs:
        # total can only be zero if no species has matches
        weighted_match_freqs[state] /= float(total)

    return individual_symbol_freqs, symbol_freqs, weighted_match_freqs

def convert_predictions(path, states):
    new_path = []
    for p in path:
        new_path.append(states[p])
    return new_path

def predict_introgressed_hmm(seqs, predict_species, index_to_species, states, \
                                 unknown_species, \
                                 match_symbol, mismatch_symbol, unknown_symbol, \
                                 expected_length_introgressed, \
                                 expected_num_introgressed_tracts):

    if unknown_species != None:
        assert len(seqs[0][0]) == len(states) - 1
        assert set(index_to_species) == set(states), str(set(index_to_species)) + ' ' + str(states)
        assert unknown_species == states[-1]

    seqs_to_predict = []
    seqs_to_predict_inds = []
    predict_inds = []

    for i in range(len(seqs)):
        if index_to_species[i] == predict_species:
            predict_inds.append(i)
            seqs_to_predict.append(seqs[i])
            seqs_to_predict_inds.append(i)

    d_freqs = get_symbol_freqs(\
        states, index_to_species, seqs, unknown_species, \
            match_symbol, mismatch_symbol, unknown_symbol)

    individual_symbol_freqs, symbol_freqs, weighted_match_freqs = \
        d_freqs[predict_species]

    for s in symbol_freqs:
        print s, symbol_freqs[s]

    hmm = HMM()

    # set obs
    hmm.set_obs(seqs_to_predict)

    # set states
    hmm.set_states(states)

    # set init
    hmm.set_init(initial_probabilities(states, weighted_match_freqs, match_symbol, \
                                           expected_length_introgressed, \
                                           expected_num_introgressed_tracts, \
                                           len(seqs[0])))

    # set emis
    hmm.set_emis(emission_probabilities(states, unknown_species, d_freqs, \
                                            match_symbol, mismatch_symbol, .99))

    # set trans
    hmm.set_trans(transition_probalities(expected_length_introgressed, \
                                             expected_num_introgressed_tracts, \
                                             predict_species, len(seqs_to_predict[0]), \
                                             weighted_match_freqs))

    # Baum-Welch parameter estimation
    hmm.go()

    predicted = []
    for i in range(len(seqs)):
        if i in predict_inds:
            hmm.set_obs(seqs[i])
            predicted.append(convert_predictions(hmm.viterbi(), states))

    return predicted, hmm

def group_actual_predicted_bases(actual, predicted, states):
    # number of bases that fit into every category of actual x
    # predicted y for all x,y in states

    # returns dictionary of lists (one entry per strain)

    d = {}
    for state_actual in states:
        for state_predicted in states:
            d[(state_actual,state_predicted)] = [0] * len(actual)

    for i in range(len(actual)):
        for j in range(len(actual[i])):
            d[(actual[i][j], predicted[i][j])][i] += 1

    return d

def group_actual_predicted_blocks(blocks_actual, blocks_predicted, states, num_samples):
    # each block entry is a list with four items: the species it was
    # predicted to be from, (dictionary of) the number of bases within
    # it that are actually from each species, the total block length,
    # and the index of the strain

    # returns dictionaries of lists (one entry per strain)
    
    # first loop through actual blocks
    d_actual_predicted = {}
    d_actual_counts = {}
    for state_actual in states:
        # this is legit because species_to samples always come first 
        d_actual_counts[state_actual] = [0] * num_samples
        for state_predicted in states:
            d_actual_predicted[(state_actual,state_predicted)] = [0] * num_samples

    for b in range(len(blocks_actual)):
        state_actual = blocks_actual[b][0]
        sample_ind = blocks_actual[b][3]
        d_actual_counts[state_actual][sample_ind] += 1
        for state_predicted in states:
            add = 0
            if blocks_actual[b][1][state_predicted] > 0:
                add = 1
            d_actual_predicted[(state_actual, state_predicted)][sample_ind] += add

    # then loop through predicted blocks
    d_predicted_actual = {}
    d_predicted_counts = {}
    for state_predicted in states:
        d_predicted_counts[state_predicted] = [0] * num_samples
        for state_actual in states:
            d_predicted_actual[(state_predicted,state_actual)] = [0] * num_samples

    for b in range(len(blocks_predicted)):
        state_predicted = blocks_predicted[b][0]
        sample_ind = blocks_predicted[b][3]
        d_predicted_counts[state_predicted][sample_ind] += 1
        for state_actual in states:
            add = 0
            if blocks_predicted[b][1][state_actual] > 0:
                add = 1
            d_predicted_actual[(state_predicted, state_actual)][sample_ind] += add

    '''
    # convert dictionaries to lists
    for state in states:
        d_actual_counts[state] = [d_actual_counts[state][x] for x in sorted(d_actual_counts[state].keys())]
        d_predicted_counts[state] = [d_predicted_counts[state][x] for x in sorted(d_predicted_counts[state].keys())]
    for key in d_actual_predicted.keys():
        d_actual_predicted[key] = [d_actual_predicted[key][x] for x in sorted(d_actual_predicted[key].keys())]
        d_predicted_actual[key] = [d_predicted_actual[key][x] for x in sorted(d_predicted_actual[key].keys())]
    '''

    return d_actual_predicted, d_predicted_actual, d_actual_counts, d_predicted_counts

def evaluate_predicted_blocks(predicted, actual, species_to, all_species):

    # each block entry is a list with four items: the species it was
    # predicted to be from, (dictionary of) the number of bases within
    # it that are actually from each species, the total block length,
    # and the index of the strain
    blocks_predicted = []
    # the same as above except for true blocks (the second entry being
    # the number of bases that are predicted for that species)
    blocks_actual = []

    for i in range(len(predicted)):

        assert len(predicted[i]) == len(actual[i]), \
            str(len(predicted[i])) + ' ' + str(len(actual[i]))

        seq_predicted = predicted[i] + ['END']
        seq_actual = actual[i] + ['END']

        current_species_predicted = seq_predicted[0]
        current_species_actual = seq_actual[0]

        # in the current predicted block, the true sequence of species states
        predicted_block_actual_sequence = []
        # in the current true block, the predicted sequence of species states
        actual_block_predicted_sequence = []

        for j in range(len(seq_predicted)):
            # === first take care of predicted block ====
            # not changing to a new block
            if seq_predicted[j] == current_species_predicted:
                predicted_block_actual_sequence.append(seq_actual[j])
            # changing to a new block (or END)
            else:
                count_species = {}
                for sp in all_species:
                    count_species[sp] = predicted_block_actual_sequence.count(sp)
                current_block = \
                    [current_species_predicted,\
                         count_species,\
                         len(predicted_block_actual_sequence),\
                         i]
                blocks_predicted.append(current_block)
                current_species_predicted = seq_predicted[j]
                predicted_block_actual_sequence = [seq_actual[j]]

            # === then take care of actual block ====
            # not changing to a new block
            if seq_actual[j] == current_species_actual:
                actual_block_predicted_sequence.append(seq_predicted[j])
            # changing to a new block (or END)
            else:
                count_species = {}
                for sp in all_species:
                    count_species[sp] = actual_block_predicted_sequence.count(sp)
                current_block = \
                    [current_species_actual,\
                         count_species,\
                         len(actual_block_predicted_sequence),\
                         i]
                blocks_actual.append(current_block)
                current_species_actual = seq_actual[j]
                actual_block_predicted_sequence = [seq_predicted[j]]

    return blocks_predicted, blocks_actual

def evaluate_predicted(predicted, actual, species_to):

    # predicted is a list; for each index that wasn't the species we
    # wanted to predict introgression, the entry is None; for other
    # indices, the entry is a list of predicted species at each
    # position (actual is the same deal)
    assert len(predicted) == len(actual), str(len(predicted)) + ' ' + str(len(actual))

    # make predictions for every cer sequence
    # and also keep track of lengths of all actual and predicted introgressed tracts
    actual_lens = []
    predicted_lens = []

    # by individual sites
    num_correct = []
    num_introgressed_correct = []

    num_predicted_tracts_actual = [] # number of predicted introgressed tracts that overlap an actual one
    num_actual_tracts_predicted = [] # number of actual introgressed tracts that overlap a predicted one

    num_introgressed_tracts = []
    num_not_introgressed_tracts = []

    num_predicted_introgressed_tracts = []
    num_predicted_not_introgressed_tracts = []

    # loop through all individuals
    for i in range(len(predicted)):

        assert len(predicted[i]) == len(actual[i]), \
            str(len(predicted[i])) + ' ' + str(len(actual[i]))

        in_ai = False # in actual introgressed region
        in_pi = False # in predicted introgressed region
        c = 0 # number of sites correct
        ci = 0 # number of introgressed sites correct
        ai_count = 0 # length of current actual introgressed tract
        pi_count = 0 # length of current predicted introgressed tract
        nai = 0 # number of actual introgressed tracts
        npi = 0 # number of predicted introgressed tracts
        hit_actual = 0 # whether we've hit an actual introgressed base
                       # in the current predicted introgressed tract
                       # (0 no, 1 yes)
        hit_predicted = 0 # whether we've hit a predicted introgressed
                          # base in the current actual introgressed
                          # tract (0 no, 1 yes)
        ap = 0
        pa = 0

        # loop through all sites (including nonpolymorphic)
        for b in range(len(predicted[i])):

            # we got it right!
            if predicted[i][b] == actual[i][b]:
                # add one to correct sites count
                c += 1
                # if introgressed, add one to correct introgressed sites count
                if actual[i][b] != species_to:
                    ci += 1
                    # if actual and predicted both introgressed, then we've found
                    # each for the current tracts
                    hit_actual = 1
                    hit_predicted = 1

            # if the current position is actually introgressed
            if actual[i][b] != species_to:
                # and we were already in an introgressed region
                if in_ai:
                    ai_count += 1
                # ...or we were not
                else:
                    in_ai = True
                    ai_count = 1
            # if the current position is not introgressed, end
            # previous actual introgressed region 
            else:
                if in_ai:
                    actual_lens.append(ai_count)
                    nai += 1
                    ap += hit_predicted
                    hit_predicted = 0
                    in_ai = False
                    ai_count = 0

            # if the current position is predicted to be introgressed
            if predicted[i][b] != species_to:
                if in_pi:
                    pi_count += 1
                else:
                    in_pi = True
                    pi_count = 1
            # if the current position is not predicted introgressed, end
            # previous predicted introgressed region 
            else:
                if in_pi:
                    predicted_lens.append(pi_count)
                    npi += 1
                    pa += hit_actual
                    hit_actual = 0
                    in_pi = False
                    pi_count = 0

        # count the last tracts if we happened to end in one
        if in_ai:
            actual_lens.append(ai_count)
            nai += 1
            ap += hit_predicted
        if in_pi:
            predicted_lens.append(pi_count)
            npi += 1
            pa += hit_actual

        # count up total number of actual and predicted introgressed
        # and not introgressed tracts
        num_introgressed_tracts.append(nai)
        num_predicted_introgressed_tracts.append(npi)

        num_actual_tracts_predicted.append(ap)
        num_predicted_tracts_actual.append(pa)

        # for *not introgressed*, check ends
        nnotai = nai - 1
        if actual[i][0] == species_to:
            nnotai += 1
        if actual[i][-1] == species_to:
            nnotai += 1
        nnotpi = npi - 1
        if predicted[i][0] == species_to:
            nnotpi += 1
        if predicted[i][-1] == species_to:
            nnotpi += 1
        num_not_introgressed_tracts.append(nnotai)
        num_predicted_not_introgressed_tracts.append(nnotpi)

        # number of sites correct
        num_correct.append(c)
        num_introgressed_correct.append(ci)

    assert(sum(num_introgressed_tracts) == len(actual_lens))
    assert(sum(num_predicted_introgressed_tracts) == len(predicted_lens))
    
    return num_correct, num_introgressed_correct, actual_lens, predicted_lens, \
        num_predicted_tracts_actual, num_actual_tracts_predicted, \
        num_introgressed_tracts, num_not_introgressed_tracts, \
        num_predicted_introgressed_tracts, num_predicted_not_introgressed_tracts

def count_blocks(blocks, negative_label, kind):
    # negative label is species to predict - i.e. if it matches
    # species to predict, then not introgressed and vice versa
    x = 0

    for block in blocks:
        if ('positive' in kind and block[0] != negative_label) or \
                ('negative' in kind and block[0] == negative_label):
            if 'true' in kind:
                if block[1] > 0:
                    x += 1
            elif 'false' in kind:
                if block[1] == 0:
                    x += 1
            else:
                x += 1
    return x

def get_positives(blocks, negative):
    x = 0
    for block in blocks:
        if block[0] != negative:
            x += 1
    return x

# seqs are NOT coded, and the order in each symbol is the order to return them in
def seq_id(seqs, index_to_species, species_order):
    # [[within 1, between 1 & 2, between 1 & 3], [within 2...]...]

    d = {} 
    for i in range(len(species_order)):
        d[species_order[i]] = i

    all_nums = []
    all_dens = []
    for i in range(len(species_order)):
        all_nums.append([0 for j in range(len(species_order))])
        all_dens.append([0 for j in range(len(species_order))])

    # loop through all pairs of sequences, and count fraction of
    # non-gap sites that match
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            seq1 = seqs[i]
            seq2 = seqs[j]
            ind1 = d[index_to_species[i]]
            ind2 = d[index_to_species[j]]
            num = 0
            den = 0
            for k in range(len(seq1)):
                if seq1[k] != '-' and seq2[k]!= '-':
                    den += 1
                    if seq1[k] == seq2[k]:
                        num += 1
            all_nums[ind1][ind2] += num
            all_nums[ind2][ind1] += num
            all_dens[ind1][ind2] += den
            all_dens[ind2][ind1] += den

    # when we're counting above, either one of the species may come
    # first; here we're averaging those two possibilities
    s = []
    for i in range(len(species_order)):
        row = []
        for j in range(len(species_order)):
            if all_dens[i][j] != 0:
                row.append(all_nums[i][j]/float(all_dens[i][j]))
            else:
                row.append(-1)
        s.append(row)
    return s

# given fractional positions for snvs and length of sequence l,
# determine integer positions; if allow_multi_hit is true, easy but if
# not, shift them around to include all the snvs
def integer_positions(positions, l, allow_multi_hit = False):

    int_positions = [int(x * l) for x in positions]
    if allow_multi_hit:
        return int_positions
    
    assert(len(positions) <= l)

    # keep first position
    a = int_positions[:1]
    # go through all other positions and move them left or right if
    # they're already taken
    for p in int_positions[1:]:
        new_p = p
        n = 1
        add = True # adding or substracting (right or left)
        while new_p in a or new_p < 0 or new_p >= l:
            if add:
                new_p = p + n
                add = False
            else:
                new_p = p - n
                add = True
                n += 1
        a.append(new_p)
    return a

def read_sim(f, num_sites, num_samples):

    # read in all the trees for blocks with no recombination within them
    t_string = f.readline()
    recomb_sites = []
    trees = []
    
    while t_string[0] == '[':
        t_start = t_string.find(']') + 1
        recomb_sites.append(int(t_string[1:t_start-1]))
        t_string = t_string[t_start:-1]
        t = parse_ms_tree(t_string)
        trees.append(t)
        t_string = f.readline()

    # read next couple of lines before sequences begin
    segsites = int(t_string[len('segsites: '):-1])
    positions = [float(x) for x in f.readline()[len('positions: '):].split()]
    # convert positions to integers
    # (allow sites to be hit multiple times because that seems reasonable)
    # (zero-indexed)
    positions = integer_positions(positions, num_sites, allow_multi_hit=True)
    # read in sequences (at this point only sites that are polymorphic)
    seqs = []
    for i in range(num_samples):
        seqs.append(f.readline()[:-1])
        assert(len(seqs[-1]) > 0)    

    return trees, recomb_sites, segsites, positions, seqs

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
def code_seqs(seqs, nsites, ref_seqs, match_symbol, mismatch_symbol, unknown_symbol, unsequenced_symbol):

    nrefs = len(ref_seqs)
    seqs_coded = []
    for seq in seqs:
        s = []
        for i in range(nsites):
            si = ''
            for r in range(nrefs):
                if ref_seqs[r][i] == unsequenced_symbol:
                    si += unknown_symbol
                elif ref_seqs[r][i] == seq[i]:
                    si += match_symbol
                else:
                    si += mismatch_symbol
            s.append(si)
            #if si == '--+':
            #    print '--+', polymorphic_sites.index(i)
            #    print  
        seqs_coded.append(s)
    return seqs_coded

# TODO change stuff to account for multiple donor species better
def make_output_dic(species, species_to, unknown_state = None):

    d = {}

    for species_from in species:
        # can skip all of these because we can get them from summing
        # different things in the categories below (i.e. to get number
        # of bases actually introgressed from par to cer, take sum of
        # num_bases_actual_par_predicted_*
        """
        if species_from != species_to:
            # introgressed bases, actual and predicted
            d['num_introgressed_bases_actual_from_' + species_from + '_to_' + species_to] = None
            d['num_introgressed_bases_predicted_from_' + species_to + '_to_' + species_to] = None

            # introgressed tracts, actual and predicted
            d['num_introgressed_tracts_actual_from_' + species_to + '_to_' + species_to] = None
            d['num_introgressed_tracts_predicted_from_' + species_to + '_to_' + species_to] = None
        """
        # introgressed tract lengths, actual and predicted
        d['tract_lengths_actual_' + species_from] = None
        d['tract_lengths_predicted_' + species_from] = None

        # number of tracts, actual and predicted
        d['num_tracts_actual_' + species_from] = None
        d['num_tracts_predicted_' + species_from] = None

        # so now we need a category for each combination of
        # predicted and actual species so that we can calculate
        # accuracy and stuff (this includes both introgressed on
        # not introgressed bases/tracts)
        for species_other in species:
            d['num_bases_actual_' + species_from + '_predicted_' + species_other] = None
            # need both of these because math
            d['num_tracts_actual_' + species_from + '_predicted_' + species_other] = None
            d['num_tracts_predicted_' + species_from + '_actual_' + species_other] = None
    
    #d['num_bases_correct_' + species_to] = None
    #d['num_tracts_correct_' + species_to] = None
    d['prob_topological_concordance'] = None
    d['prob_monophyletic_concordance'] = None
    d['fraction_concordant'] = None

    for i in range(len(species)):
        for j in range(i, len(species)):
            d['avg_identity_' + species[i] + '_' + species[j]] = None


    # HMM parameters
    for i in range(len(species)):
        d['init_' + species[i]] = None
                   
    for i in range(len(species)):
        for j in range(len(species)):
            d['trans_' + species[i] + '_' + species[j]] = None

    for i in range(len(species)):
        to_states = species
        if unknown_state != None:
            to_states.append(unknown_state)
        for j in range(len(to_states)):
            d['emis_' + species[i] + '_' + species[j]] = None

    return d

def write_output_line(f, d, header_line):

    items = sorted(d.keys())
    if not header_line:
        items = [str(d[x]) for x in items]

    keys = sorted(d.keys())
    i = 0
    for item in items[:-1]:
        f.write(item + '\t')
        print keys[i]
        print item
        i += 1
    f.write(items[-1] + '\n')
