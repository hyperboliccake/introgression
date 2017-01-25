import os
import sys
import copy
import itertools

sys.path.append('../hmm')

from hmm_bw import *

theory = False
theory_done = True

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

    num_reps = int(sys.argv[11])

    return tag, model, N0, include_bay, include_unk,\
        num_samples_cer, num_samples_par, num_samples_bay,\
        par_cer_migration, bay_cer_migration,\
        t_par_cer, t_bay_par_cer,\
        num_sites, rho, outcross_rate, num_reps

def parse_tree_helper(t):
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
    return [parse_tree_helper(left), parse_tree_helper(right), time]

# converts newick tree string to nested list format
def parse_tree(t):
    return parse_tree_helper(t[:-1] + ':0')

# converts branch lengths to represent total time rather than just
# that length
def make_times_additive(t):
    if len(t) == 2:
        return t
    left = make_times_additive(t[0])
    right = make_times_additive(t[1])

    # do this slightly weird thing with the max to deal with rounding
    # (make sure times are increasing)
    prev = None
    if len(left) == 2:
        if len(right) == 2:
            prev = max(left[1], right[1])
        else:
            prev = max(left[1], right[2])
    else:
        if len(right) == 2:
            prev = max(left[2], right[1])
        else:
            prev = max(left[2], right[2])

    return [left, right, t[2] + prev]

# actual introgressed with 2 species total, within a single unrecombined block
def find_introgressed_2(t, cutoff_time, to_species, from_species, index_to_species):
    lineages = split(t, cutoff_time)
    introgressed = index_to_species
    for l in lineages:
        # if there's any occurence of the from_species in this clade,
        # then mark all individuals as coming from from_species (since
        # for now we're only allowing migration in one direction)
        if not is_partial_clade(l, to_species, index_to_species):
            for label in get_labels(l):
                introgressed[label] = from_species

    # replace the ones that were never going to be introgressed with empty list
    for i in range(len(index_to_species)):
        if index_to_species[i] == from_species:
            introgressed[i] == []

    return introgressed, len(lineages)

# actual introgressed with 3 species total
def find_introgressed_3(t, cutoff_time_1, cutoff_time_2, to_species, \
                            from_species_1, from_species_2, label_to_species):
    lineages = split(t, cutoff_time)
    non_introgressed = []
    introgressed = []
    for l in lineages:
        if is_partial_clade(l, label_to_species):
            non_introgressed += get_labels(l)
        else:
            labels = get_labels(l)
            for label in labels:
                if label_to_species[label] == from_species:
                    introgressed.append(label)
    return introgressed, len(lineages)

#TODO?
#def is_concordant

# predict which sequences are introgressed based on sequence id;
# results returned in order of indices in p1 + p2; returns predicted bases introgressed for 
def predict_introgressed(seqs_filled, window_size, window_shift, thresholds = [.7, .95, .96]):
    
    predicted = []
    
    for seq in seqs_filled:
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
                    p[i:i+window_size] = [1] * window_size
        predicted.append(p)

    return predicted

def initial_probalities(states, individual_symbol_freqs, match_symbol):
    init = []
    for state in states:
        init.append(individual_symbol_freqs[state][match_symbol])
    return init

def emission_probalities(states, symbol_freqs, match_symbol, mismatch_symbol, own_bias):

    emis = []
    for i in range(len(states)):
        emis_species = {}
        for symbol in symbols:
            p = symbol_freqs[symbol]
            if symbol[0] == match_symbol:
                p *= own_bias
            elif symbol[0] == mismatch_symbol:
                p *= (1 - own_bias)
            else:
                p *= ((own_bias) * individual_symbol_freqs[i][match_symbol] + \
                          (1 - own_bias) * individual_symbol_freqs[i][mismatch_symbol])
            emis_species[symbol] = p
        norm = float(sum(emis_species.values()))
        for symbol in emis_species:
            emis_species[symbol] /= norm
        emis.append(emis_species)
    return emis

def transition_probalities(expected_length_introgressed, expected_num_introgressed_tracts, predict_species, n, symbol_freqs, unique_match_freqs):

    states_not_predict = states
    states_not_predict.remove(predict_species)
    
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
    #expected_length_introgressed = {}
    expected_num_introgressed_bases = {}
    for species in expected_num_introgressed_tracts:
        expected_num_introgressed_bases = expected_num_introgressed_tracts[species] * \
            expected_length_introgressed[species]
    # TODO should be able to calculate expected length based on expected amount of migration?
    #expected_num_introgressed_tracts = {}

    expected_num_not_introgressed_tracts = \
        sum(expected_num_introgressed_tracts.values()) + 1
    expected_num_not_introgressed_bases = \
        n - sum(expected_num_introgressed_bases.values())
    expected_length_not_introgressed = float(expected_num_not_introgressed_bases) / \
        expected_num_not_introgressed_tracts

    # fraction of time we should choose given species state over
    # others based on number of sites that match it
    fracs = {}
    frac_den = 0
    for state in expected_length_introgressed:
        fracs[state] = unique_match_freqs[state]
        frac_den += fracs[state]
    # TODO: question: if A:10 and B:12 and all of those matches
    # overlap, then should B get it 12/22 of the time or all the time?
    # i.e. is it just the sites that uniquely match that are relevant?
    # i think yes, the unique ones
    for state in fracs:
        fracs[state] /= float(frac_den)

    trans = []
    for state_from in states:
        trans_current = {}
        total = 0

        for state_to in states:
            val = -1
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

        trans_current[trans_from] = 1 - total
        trans_row = []
        for state_to in states:
            trans_row.append(trans_current[state_to])
        trans.append(trans_row)

    return trans

def symbol_freqs(states, seqs, match_symbol, mismatch_symbol, unknown_symbol):

    known_symbols = [match_symbol, mismatch_symbol]
    symbols = known_symbols + [unknown_symbol]

    # for species a: +, -, ?; species b: +, -, ?...etc
    individual_symbol_freqs = {}
    for i in range(len(states)):
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
    symbol_combinations = [''.join(x) for x in list(itertools.combinations_with_replacement(symbols, len(states)))]
    symbol_freqs = dict(zip(symbol_combinations, [0]*len(symbol_combinations)))
    for symbol in symbols:
        num = 0
        den = 0
        for seq in seqs:
            num += seq.count(symbol)
            den += len(seq)
        symbol_freqs[symbol] = float(num) / den

    # for each species, fraction of unique matches to that species,
    # i.e. for the first species count all +--, +-?, +?-, +??
    unique_match_freqs = {}
    total = 0
    for i in range(len(states)):
        current = 0
        symbols_to_count = []
        for s in symbol_combinations:
            keep = True
            for j in range(len(states)):
                if i == j:
                    if s[j] != match_symbol:
                        keep = False
                        break
                else:
                    if s[j] == match_symbol:
                        keep = False
                        break
            if keep:
                symbols_to_count.append(s)
        for seq in seqs:
            for symbol in symbols_to_count:
                current += seq.count(symbol)
        unique_match_freqs[states[i]] = current
        total += current
    for state in unique_match_freqs:
        unique_match_freqs /= float(total)

    return individual_symbol_freqs, symbol_freqs, unique_match_freqs

def convert_predictions(path, states):
    new_path = []
    for p in path:
        new_path.append(states[p])
    return new_path

def predict_introgressed_hmm(seqs, predict_species, index_to_species, states, match_symbol, mismatch_symbol, unknown_symbol):

    seqs_to_predict = []
    seqs_to_predict_inds = []
    predict_inds = []
    for i in range(len(seqs)):
        if index_to_species[i] == predict_species:
            predict_inds.append(i)
            seqs_to_predict.append(seqs[i])
            seqs_to_predict_inds.append(i)

    individual_symbol_freqs, symbol_freqs, unique_match_freqs = symbol_freqs(states, seqs_to_predict, \
                                                                                 match_symbol, mismatch_symbol, unknown_symbol)

    hmm = HMM()

    # set obs
    hmm.set_obs(seqs_to_predict)

    # set states
    hmm.set_states(states)

    # set init
    hmm.set_init(initial_probabilities(states, individual_symbol_freqs, match_symbol))

    # set emis
    hmm.set_emis(emission_probabilites(states, symbol_freqs, match_symbol, mismatch_symbol, .99))

    # set trans
    hmm.set_trans(transition_probalities(expected_length_introgressed, expected_num_introgressed_tracts, \
                                             predict_species, n, unique_match_freqs))

    # Baum-Welch parameter estimation
    hmm.go()

    predicted = []
    for i in range(len(seqs)):
        if i in predict_inds:
            hmm.set_obs(seqs[i])
            predicted.append(convert_predictions(hmm.viterbi()))
        else:
            predicted.append([])

    return predicted, hmm

def evaluate_predicted_blocks(predicted, actual, index_to_species):

    # each block entry is a list with four items: the species it was
    # predicted to be from, the number of bases within it that are
    # actually from that species, the total block length, and the index of
    # the strain
    blocks_predicted = []
    # the same as above except for true blocks (the second entry being
    # the number of bases that are predicted for that species)
    blocks_actual = []

    # include only the individuals we made predictions for
    predicted_indices = []
    for i in range(len(predicted)):
        if predicted[i] != []:
            predicted_indices.append(i)

    for i in predicted_indices:

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
                current_block = \
                    [current_species_predicted,\
                         predicted_block_actual_sequence.count(current_species_predicted),\
                         len(predicted_block_actual_sequence),\
                         i]
                blocks_predicted.append(current_block)
                current_species_predicted = seq_predicted[j]
                predicted_block_actual_sequence = [actual_seq[j]]

            # === then take care of actual block ====
            # not changing to a new block
            if seq_actual[j] == current_species_actual:
                actual_block_predicted_sequence.append(seq_predicted[j])
            # changing to a new block (or END)
            else:
                current_block = \
                    [current_species_actual,\
                         actual_block_predicted_sequence.count(current_species_actual),\
                         len(actual_block_predicted_sequence),\
                         i]
                blocks_actual.append(current_block)
                current_species_actual = seq_actual[j]
                actual_block_predicted_sequence = [predicted_seq[j]]

    return blocks_predicted, blocks_actual

def evaluate_predicted(predicted, actual):

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

        assert len(predicted[i]) == len(actual[i]), str(len(predicted[i])) + ' ' + str(len(actual[i]))

        in_ai = False # in actual introgressed region
        in_pi = False # in predicted introgressed region
        c = 0 # number of sites correct
        ci = 0 # number of introgressed sites correct
        ai_count = 0 # length of current actual introgressed tract
        pi_count = 0 # length of current predicted introgressed tract
        nai = 0 # number of actual introgressed tracts
        npi = 0 # number of predicted introgressed tracts
        hit_actual = 0 # whether we've hit an actual introgressed base in the current predicted introgressed tract (0 no, 1 yes)
        hit_predicted = 0 # whether we've hit a predicted introgressed base in the current actual introgressed tract (0 no, 1 yes)
        ap = 0
        pa = 0

        # loop through all sites (including nonpolymorphic)
        for b in range(len(predicted[i])):

            # we got it right!
            if predicted[i][b] == actual[i][b]:
                # add one to correct sites count
                c += 1
                # if introgressed, add one to correct introgressed sites count
                if actual[i][b] == 1:
                    ci += 1
                    # if actual and predicted both introgressed, then we've found
                    # each for the current tracts
                    hit_actual = 1
                    hit_predicted = 1

            # if the current position is actually introgressed
            if actual[i][b] == 1:
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
            if predicted[i][b] == 1:
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

        # count up total number of actual and predicted introgressed and not introgressed tracts
        num_introgressed_tracts.append(nai)
        num_predicted_introgressed_tracts.append(npi)

        num_actual_tracts_predicted.append(ap)
        num_predicted_tracts_actual.append(pa)

        # for *not introgressed*, check ends
        nnotai = nai - 1
        if actual[i][0] == 0:
            nnotai += 1
        if actual[i][-1] == 0:
            nnotai += 1
        nnotpi = npi - 1
        if predicted[i][0] == 0:
            nnotpi += 1
        if predicted[i][-1] == 0:
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

    s = []
    for i in range(len(species_order)):
        s.append([[] for j in range(len(species_order))])

    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            seq1 = seqs[i]
            seq2 = seqs[j]
            ind1 = d[index_to_species[i]]
            ind2 = d[index_to_species[i]]
            num = 0
            den = 0
            for k in range(len(seq1)):
                if seq1[k] != '-' and seq2[k]!= '-':
                    den += 1
                    if seq1[k] == seq2[k]:
                        num += 1
            t = float(num) / den
            s[ind1][ind2].append(t)
            
    for i in range(len(species_order)):
        row = []
        for j in range(i, len(species_order)):
            num = sum(s[i][j]) + sum(s[j][i])
            den = len(s[i][j]) + len(s[j][i])
            t = float(num) / den
            s[i][j] = t
            s[j][i] = t

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

def read_sim(f):

    # read in all the trees for blocks with no recombination within them
    t_string = f.readline()
    while t_string[0] == '[':
        t_start = t_string.find(']') + 1
        recomb_sites.append(int(t_string[1:t_start-1]))
        t_string = t_string[t_start:-1]
        t = parse_tree(t_string)
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
    polymorphic_sites = set(polymorphic_sites)
    for seq in polymorphic_seqs:
        s = ''
        for i in range(nsites):
            poly_ind = 0
            if i in polymorphic_sites:
                s += seq[poly_ind]
                poly_ind += 1
            else:
                s += fill
        seqs_filled.append(s)
    return seqs_filled

# convert sequences from bases to symbols indicating which references they match
def code_seqs(seqs, polymorphic_sites, nsites, ref_seqs, match_symbol, mismatch_symbol, unknown_symbol, unsequenced_symbol):

    nrefs = len(ref_seqs)
    seqs_coded = []
    for seq in seqs:

        s = []
        for i in range(nsites):
            si = ''
            poly_ind = 0
            if i in polymorphic_sites:
                for r in range(nrefs):
                    if ref_seqs[r][i] == unsequenced_symbol:
                        si += unknown_symbol
                    elif ref_seqs[r][i] == seq[poly_ind]:
                        si += match_symbol
                    else:
                        si += mismatch_symbol
                poly_ind += 1
            else:
                for r in range(nrefs):
                    if ref_seqs[r][i] == unsequenced_symbol:
                        si += unknown_symbol
                    else:
                        si += match_symbol
            s.append(si)
        seqs_coded.append(s)
    return seqs_coded

# TODO change stuff to account for multiple donor species better
def make_output_dic(species, species_to_predict, unknown_state = None):

    d = {'num_introgressed_bases_actual_' + species_to_predict:None,\
             'num_introgressed_tracts_actual_' + species_to_predict:None,\
             'num_not_introgressed_tracts_actual_' + species_to_predict:None,\
             'num_introgressed_bases_predicted_' + species_to_predict:None,\
             'num_introgressed_tracts_predicted_' + species_to_predict:None,\
             'num_not_introgressed_tracts_predicted_' + species_to_predict:None,\
             'num_bases_correct_' + species_to_predict:None,\
             'num_actual_introgressed_bases_predicted_' + species_to_predict:None,\
             'num_actual_introgressed_tracts_predicted_' + species_to_predict:None,\
             'introgressed_tract_lengths_actual_' + species_to_predict:None,\
             'introgressed_tract_lengths_predicted_' + species_to_predict:None,\
             'prob_topological_concordance':None,\
             'prob_monophyletic_concordance':None,\
             'fraction_concordant':None)

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
        for j in range(to_states):
            d['emis_' + species[i] + '_' + species[j]] = None

def write_output_line(f, d, header_line):

    items = sorted(d.keys())
    if not header_line:
        items = [str(d[x]) for x in items]

    for item in items[:-1]:
        f.write(item + '\t')
    f.write(items[-1] + '\n')
