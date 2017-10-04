def find_in_blocks(i, blocks):
    # return the state that has i in one of its blocks
    for state in blocks:
        for block in blocks[state]:
            if i >= block[0] and i <= block[1]:
                return state

def count_bases_one(d1, d2, args):
    # for one individual
    d = {}
    for state1 in args['states']:
        for state2 in args['states']:
            d[(state1, state2)] = 0

    for i in range(args['num_sites']):
        state1 = find_in_blocks(i, d1)
        state2 = find_in_blocks(i, d2)
        d[(state1, state2)] += 1
    return d

def count_bases(d1, d2, args):
    # separate counts for all individuals
    d = {}
    # average counts across all individuals
    d_avg = {}
    for state1 in args['states']:
        for state2 in args['states']:
            d_avg[(state1, state2)] = 0
    # loop through all individuals
    num_inds = 0
    for ind in d1.keys():
        if d2.has_key(ind):
            d[ind] = count_bases_one(d1[ind], d2[ind], args)
            for pair in d[ind]:
                d_avg[pair] += d[ind][pair]
            num_inds += 1
    for pair in d_avg:
        d_avg[pair] = d_avg[pair]/float(num_inds)
    return d, d_avg

def write_compare_header(f, states, suffix1, suffix2, sep='\t'):

    header_string = ''
    for state1 in states:
        for state2 in states:
            header_string += 'bases_' + suffix1 + '_' + state1 + \
                             '_' + suffix2 + '_' + state2 + sep
    f.write(header_string[:-len(sep)] + '\n')

def write_compare_line(avg_base_counts, f, states, suffix1, suffix2, sep='\t'):

    line_string = ''
    for state1 in states:
        for state2 in states:
            line_string += str(avg_base_counts[(state1, state2)]) + sep

    f.write(line_string[:-len(sep)] + '\n')


###############


def group_actual_predicted_blocks(blocks_actual, blocks_predicted, states, inds_to_predict):
    # each block entry is a list with four items: the species it was
    # predicted to be from, (dictionary of) the number of bases within
    # it that are actually from each species, the total block length,
    # and the index of the strain

    # returns dictionaries of lists (one entry per strain)
    
    num_samples = len(inds_to_predict)

    # first loop through actual blocks
    d_actual_predicted = {}
    d_actual_counts = {}
    for state_actual in states:
        # this is legit because species_to samples always come first 
        d_actual_counts[state_actual] = dict(zip(inds_to_predict, [0] * num_samples))
        for state_predicted in states:
            d_actual_predicted[(state_actual,state_predicted)] = dict(zip(inds_to_predict, [0] * num_samples))

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
        d_predicted_counts[state_predicted] = dict(zip(inds_to_predict, [0] * num_samples))
        for state_actual in states:
            d_predicted_actual[(state_predicted,state_actual)] = dict(zip(inds_to_predict, [0] * num_samples))

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

    # convert to dictionaries containing lists for each individual

    for state_actual in states:
        d_actual_counts[state_actual] = [d_actual_counts[state_actual][i] for i in inds_to_predict]
        # looks dumb but correct
        d_predicted_counts[state_actual] = [d_predicted_counts[state_actual][i] for i in inds_to_predict]
        for state_predicted in states:
            d_actual_predicted[(state_actual, state_predicted)] = \
                [d_actual_predicted[(state_actual, state_predicted)][i] for i in inds_to_predict]
            d_predicted_actual[(state_predicted, state_actual)] = \
                [d_predicted_actual[(state_predicted, state_actual)][i] for i in inds_to_predict]

    return d_actual_predicted, d_predicted_actual, d_actual_counts, d_predicted_counts

def evaluate_predicted_blocks(predicted, actual, species_to, all_species, inds_to_predict):

    # each block entry is a list with four items: the species it was
    # predicted to be from, (dictionary of) the number of bases within
    # it that are actually from each species, the total block length,
    # and the index of the strain
    blocks_predicted = []
    # the same as above except for true blocks (the second entry being
    # the number of bases that are predicted for that species)
    blocks_actual = []

    for i in range(len(inds_to_predict)):
        strain_id = inds_to_predict[i]
        
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
                         strain_id]
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
                         strain_id]
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
