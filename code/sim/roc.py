import compare_introgressed

def reformat_probs(probs):
    # convert probabilities from {1:{cer:.9,.9,..., par:.1,.1,...}}
    # format to {1:[{cer:.9, par:.1},{cer:.9, par:.1},...]} format
    p = {}
    for ind in probs.keys():
        p_ind = []
        for i in range(len(probs[ind][probs[ind].keys()[0]])):
            d = {}
            for state in probs[ind].keys():
                d[state] = probs[ind][state][i]
            p_ind.append(d)
        p[ind] = p_ind
    return p

def get_stats(actual, predicted, sim_args):

    d, d_avg = compare_introgressed.count_bases(actual, predicted, \
                                                sim_args, 'actual', 'predicted', \
                                                sim_args['species_from1'])
    return d_avg

def threshold_probs(probs, threshold, default_state):

    # when there's just one non-default state, choose that one if
    # it's over the threshold (though not necessarily over .5);
    # when there's more than one, choose the highest-probability
    # non-default state (?) -> this is weird when it's more
    # complicated than binary

    predicted_thresholded = []

    positive_states = probs[0].keys()
    positive_states.remove(default_state)

    for i in range(len(probs)):
        # only call a negative if total of all positive states
        # don't exceed threshold
        if (1 - probs[i][default_state]) < threshold:
            predicted_thresholded.append(default_state)
        else:
            max_positive_state = None
            max_positive_prob = -1
            for state in positive_states:
                if probs[i][state] > max_positive_prob:
                    max_positive_state = state
                    max_positive_prob = probs[i][state]
            predicted_thresholded.append(max_positive_state)
    return predicted_thresholded
    
def write_roc_header(f, stats, sep):

    f.write('threshold')
    for key in sorted(stats.keys()):
        if type(key) == type(()):
            key = '_'.join(key)
        f.write(sep + str(key))
        f.flush()
    f.write('\n')

def write_roc_line(f, threshold, stats, header):

    sep = '\t'
    if header:
        write_roc_header(f, stats, sep)
    f.write(str(threshold))
    for key in sorted(stats.keys()):
        f.write(sep + str(stats[key]))
    f.write('\n')
    f.flush()
