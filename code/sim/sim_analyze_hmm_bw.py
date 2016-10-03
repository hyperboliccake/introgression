import os
import sys
import copy

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

# finds the lineages that exist at a given time
def split(t, cutoff_time):
    t = make_times_additive(t)
    current_lineages = [t]
    final_lineages = []
    while len(current_lineages) > 0:
        l = current_lineages.pop()
        if len(l) == 2 or l[2] < cutoff_time:
            final_lineages.append(l)
        else:
            left = l[0]
            right = l[1]
            prev_time = -1
            if len(left) == 2:
                prev_time = left[1]
            else:
                prev_time = left[2]
            if prev_time < cutoff_time:
                final_lineages.append(l)
            else:
                current_lineages.append(left)
                current_lineages.append(right)
    return final_lineages

# return list of lineages that exist at given time
def get_labels(t):
    if len(t) == 2:
        return [t[0]]
    return get_labels(t[0]) + get_labels(t[1])

# checks whether t consists _only_ of labels in A (but does not
# necessarily include all of them)
def is_partial_clade(t, A):
    labels = get_labels(t)
    for l in labels:
        if l not in A:
            return False
    return True

# checks whether t consists _only_ of labels in A or _only_ of
# labels not in A
def is_one_species(t, A):
    labels = get_labels(t)
    is_A = False
    if labels[0] in A:
        is_A = True
    for l in labels:
        if (l in A) != is_A:
            return False
    return True

# checks whether t contains all labels in A _and_ only labels in A
def is_whole_species(t, A):
    labels = get_labels(t)
    if len(labels) != len(A):
        return False
    for l in labels:
        if l not in A:
            return False
    return True

# actual introgressed from set A according to coalescent tree
def find_introgressed(t, cutoff_time, A, B = []):
    lineages = split(t, cutoff_time)
    non_introgressed = []
    introgressed = []
    for l in lineages:
        if is_partial_clade(l, A + B):
            non_introgressed += get_labels(l)
        else:
            labels = get_labels(l)
            for label in labels:
                if label in A:
                    introgressed.append(label)
    return introgressed, len(lineages)

def get_internal_nodes(t):
    assert(type(t) == type([]))
    assert(len(t) == 3)

    nodes = [t]
    if len(t[0]) == 3:
        nodes += get_internal_nodes(t[0])
    if len(t[1]) == 3:
        nodes += get_internal_nodes(t[1])
    return nodes

def get_species(t, label_to_species):
    labels = get_labels(t)
    species = []
    for l in labels:
        s = label_to_species[l]
        if s not in species:
            species.append(s)
    return species

def collapse_tree(t, label_to_species):

    # collapse gene tree by proceeding backwards in time until a
    # coalescence occurs between two different species. group the two
    # species into a clade. continue backwards until a coalescence has
    # occurred between two clades. if both clades involve have already
    # experienced inter-clade coalescences, ignore the
    # event. otherwise, group the two clades into a larger
    # clade. proceed backwards until all species have experience inter
    # clade coalescences.

    t = make_times_additive(t)
    # get all coalescent events
    internal_nodes = get_internal_nodes(t)
    assert(len(internal_nodes) == len(label_to_species) - 1)
    # sort by time that they occurred (most recent at beginning of list)
    # note that we have to go down one level to get the time of the node
    def key_function(x):
        left = None
        if len(x[0]) == 2:
            left = x[0][1]
        else:
            left = x[0][2]
        right = None
        if len(x[1]) == 2:
            right = x[1][1]
        else:
            right = x[1][2]
        return max(left, right)
    internal_nodes.sort(key = key_function)

    for it in internal_nodes:
        s = get_species(it, label_to_species)
        if len(s) == 2:
            for l in label_to_species:
                if label_to_species[l] in s:
                    label_to_species[l] = s
        else:
            assert(len(s) == 1)
    return label_to_species.values()[0]

def sort_recursively(a):
    if type(a) != type([]):
        return a
    left = sort_recursively(a[0])
    right = sort_recursively(a[1])
    if left < right:
        return [left, right]
    else:
        return [right, left]

def equivalent_topologies(a, b):
    return sort_recursively(a) == sort_recursively(b)

def is_concordant(t, A, B):

    return (is_whole_species(t[0], A) and is_whole_species(t[1], B)) or \
        (is_whole_species(t[1], A) and is_whole_species(t[0], B))

def is_monophyletically_concordant(t, A, B, C, split, label_to_species = None, species_topology = None):
    if split:
        if not is_split_concordant(t, B, C):
            return False
    else: 
        if not is_topologically_concordant(t, label_to_species, species_topology):
            return False
    nodes = get_internal_nodes(t)
    clade_A = False
    clade_B = False
    clade_C = False
    for node in nodes:
        if is_whole_species(node, A):
            clade_A = True
        elif is_whole_species(node, B):
            clade_B = True
        elif is_whole_species(node, C):
            clade_C = True
    return clade_A and clade_B and clade_C

def is_topologically_concordant(t, label_to_species, species_topology):

    gene_topology = collapse_tree(t, label_to_species)
    if equivalent_topologies(species_topology, gene_topology):
        return True
    return False

# can do this simpler thing for the special case of 3 species
def is_split_concordant(t, A, B = []):

    return (is_one_species(t[0], A + B) and is_one_species(t[1], A + B))

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

def predict_introgressed_hmm(seqs_filled):

    # create a hidden markov model and determine which reference genome we
    # are most likely to be in at each variant site
    hmm = HMM()
        
    hmm.set_obs(seqs_filled)
    #hmm.set_obs([seqs_filled[0], seqs_filled[1]])
    #hmm.set_obs(['0000001000000000010000000','000011111111111111122222'])
    #hmm.set_obs(['012', '000'])
    #hmm.set_obs(['00', '01'])
    hmm.set_init([.5,.5])
    """
    fo = open('obs.txt', 'w')
    for os in seqs_filled:
        for o in os:
            fo.write(o)
        fo.write('\n')
    fo.close()
    
    fo = open('int.txt', 'w')
    for os in seqs_filled:
        for o in os:
            fo.write(o)
        fo.write('\n')
    fo.close()
    """
    # 0,1 would also work here
    hmm.set_states(['cer', 'par'])
    # a little weird because we know recombination rate we used
    hmm_trans = .0005
    
    # hmm.set_trans({'cer':{'cer':1-hmm_trans, 'par':hmm_trans},
    # 'par':{'cer':hmm_trans, 'par':1-hmm_trans}})
    hmm.set_trans([[1-hmm_trans,hmm_trans],[hmm_trans,1-hmm_trans]])

    # combined error in parent sequences and observed sequence -
    # in this case 'error' comes from amout of difference
    # reasonable to see between reference and other strains of
    # same species; remember we've coded 0 for cer (non
    # introgressed) and 1 for par (introgressed)
    # hmm.set_emis({'cer':{0:.95, 1:.05},'par':{0:.05, 1:.95}}) 
    hmm.set_emis([{'0':.5, '1':.0001, '2':.4998, '3':.0001},{'0':.0001, '1':.5, '2':.4998, '3':.0001}])

    # Baum-Welch parameter estimation
    hmm.go()

    predicted = []
    for i in range(len(seqs_filled)):
        obs_seq = seqs_filled[i]
        hmm.set_obs(obs_seq)
        predicted.append(hmm.viterbi())

    return predicted, hmm

def evaluate_predicted_blocks(predicted, actual):
    # in the sequences, 1 indicates introgressed (actual or predicted)
    # and 0 indicates not

    # each block entry is a list with four items: whether it is
    # predicted to be introgressed (n for no, i for yes), the number
    # of bases within it that are actually introgressed, the total
    # block length, and the index of the strain
    blocks_predicted = []
    # the same as above except for true blocks
    blocks_actual = []

    for i in range(len(predicted)):

        assert len(predicted[i]) == len(actual[i]), \
            str(len(predicted[i])) + ' ' + str(len(actual[i]))

        in_predicted_introgressed_block = (predicted[i][0] == 1)
        in_actual_introgressed_block = (actual[i][0] == 1)

        predicted_block_actual_sequence = []
        actual_block_predicted_sequence = []

        predicted_seq = predicted[i] + ['END']
        actual_seq = actual[i] + ['END']

        for j in range(len(predicted_seq)):
            # === first take care of predicted block ====
            # base not predicted to be introgressed
            if predicted_seq[j] == 0 or \
                    (predicted_seq[j] == 'END' and in_predicted_introgressed_block):
                # end predicted introgressed block
                if in_predicted_introgressed_block:
                    current_predicted_block = \
                        ['i', predicted_block_actual_sequence.count(1), \
                             len(predicted_block_actual_sequence), i]
                    blocks_predicted.append(current_predicted_block)
                    in_predicted_introgressed_block = False
                    predicted_block_actual_sequence = [actual_seq[j]]
                # add on to current nonintrogressed block
                else:
                    predicted_block_actual_sequence.append(actual_seq[j])
            # base predicted to be introgressed
            else:
                # add on to current introgressed block
                if in_predicted_introgressed_block:
                    predicted_block_actual_sequence.append(actual_seq[j])
                # end predicted nonintrogressed block
                else:
                    current_predicted_block = \
                        ['n', predicted_block_actual_sequence.count(1), \
                             len(predicted_block_actual_sequence), i]
                    blocks_predicted.append(current_predicted_block)
                    in_predicted_introgressed_block = True
                    predicted_block_actual_sequence = [actual_seq[j]]
                
            # === then take care of actual block ====
            # base not actually introgressed
            if actual_seq[j] == 0 or \
                    (actual_seq[j] == 'END' and in_actual_introgressed_block):
                # end actual introgressed block
                if in_actual_introgressed_block:
                    current_actual_block = \
                        ['i', actual_block_predicted_sequence.count(1), \
                             len(actual_block_predicted_sequence), i]
                    blocks_actual.append(current_actual_block)
                    in_actual_introgressed_block = False
                    actual_block_predicted_sequence = [predicted_seq[j]]
                # add on to current nonintrogressed block
                else:
                    actual_block_predicted_sequence.append(predicted_seq[j])
            # base actually introgressed
            else:
                # add on to introgressed block
                if in_actual_introgressed_block:
                    actual_block_predicted_sequence.append(predicted_seq[j])
                # add on to introgressed block
                else:
                    current_actual_block = \
                        ['n', actual_block_predicted_sequence.count(1), \
                             len(actual_block_predicted_sequence), i]
                    blocks_actual.append(current_actual_block)
                    in_actual_introgressed_block = True
                    actual_block_predicted_sequence = [predicted_seq[j]]

    return blocks_predicted, blocks_actual

def evaluate_predicted(predicted, actual):


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


tag = sys.argv[1]
outdir = '../../results/sim/'
outfilename = 'sim_out_' + tag + '.txt'
results_filename = 'sim_out_' + tag + '_summary.txt'

model = sys.argv[2]

N0 = int(sys.argv[3])

num_samples_par = int(sys.argv[4])
num_samples_cer = int(sys.argv[5])
num_samples = num_samples_par + num_samples_cer

# migration parameter is 2 * N0 * m, where mij is fraction of i made
# up of j each generation; need to figure out how to make migration
# rates equivalent for different models
par_cer_migration = 2 * N0 * float(sys.argv[6])

# in generations
t_par_cer = float(sys.argv[7]) / (2 * N0)

# 13,500 sites to get about 10% with one recombination event, .3% with
# more than one (based on poisson(.1), 1 recombination per chromosome
# of average length 750,000)
num_sites = int(sys.argv[8])

# parameter is recombination rate between adjacent bp per generation
# should probably be 1/750000 + 6.1 * 10^-6 (where 750000 is average
# chr size)
rho = 2 * N0 * float(sys.argv[9]) * (num_sites - 1)

outcross_rate = float(sys.argv[10])

rho *= outcross_rate

# estimate from humans
mu = 1.84 * 10 ** -10
theta = mu * 2 * num_sites * N0

num_reps = int(sys.argv[11])

# take first index from each population to be reference sequence
# [s288c is part of 41-strain group, which we're saying has migration happen]
ref_ind_par = 0
ref_ind_cer = num_samples_par 

#####
# stats to keep track of
#####

# for each sim, number of introgressed bases in each strain
num_introgressed_cer = [[0 for i in range(num_samples_cer)] for n in range(num_reps)]

# success of HMM predictions
fraction_correct_all = [[] for i in range(num_reps)]

fraction_introgressed_correct_all = [[] for i in range(num_reps)]

predicted_lens_all = [[] for i in range(num_reps)]

# success of windowed id predictions
fraction_correct_all_window = [[] for i in range(num_reps)]

fraction_introgressed_correct_all_window = [[] for i in range(num_reps)]

predicted_lens_all_window = [[] for i in range(num_reps)]

# actual lengths don't depend on prediction strategy
actual_lens_all = [[] for i in range(num_reps)]

fpr_cer = [0] * num_reps
fpr_cer_only_par = [0] * num_reps

# id between all cer and cer ref
avg_id_cer = [0] * num_reps

# id between all cer and par ref 
avg_id_cer_par = [0] * num_reps

num_lineages_at_join = [[] for n in range(num_reps)]

concordant = [[] for n in range(num_reps)]

prob_topological_concordance = 0

# windows for looking for introgressed sequence
window_size = 1000
window_shift = 500

#####
# loop through all reps
#####

# write results
fout = open(outdir + results_filename, 'w')
print 'writing to', outdir + results_filename
fout.write('num_introgressed_cer\t' + \
               'num_introgressed_tracts_cer\t' + \
               'num_not_introgressed_tracts_cer\t' + \
               'num_predicted_introgressed_cer\t' + \
               'num_predicted_introgressed_tracts_cer\t' + \
               'num_predicted_not_introgressed_tracts_cer\t' + \
               'num_predicted_introgressed_cer_window\t' + \
               'num_predicted_introgressed_tracts_cer_window\t' + \
               'num_predicted_not_introgressed_tracts_cer_window\t' + \
               'num_correct\t' + \
               'num_introgressed_correct\t' + \
               'num_predicted_tracts_actual\t' + \
               'num_actual_tracts_predicted\t' + \
               'num_correct_window\t' + \
               'num_introgressed_correct_window\t' + \
               'num_predicted_tracts_actual_window\t' + \
               'num_actual_tracts_predicted_window\t' + \
               'actual_lens\t' + \
               'predicted_lens\t' + \
               'predicted_lens_window\t' + \
               'avg_id_cer\t' + \
               'avg_id_cer_par\t' + \
               'num_lineages\t' + \
               'fraction_concordant\t' + \
               'init_cer\t' + \
               'init_par\t' + \
               'trans_cer_cer\t' + \
               'trans_cer_par\t' + \
               'trans_par_cer\t' + \
               'trans_par_par\t' + \
               'emis_cer_cer\t' + \
               'emis_cer_par\t' + \
               'emis_cer_np\t' + \
               'emis_par_cer\t' + \
               'emis_par_par\t' + \
               'emis_par_np\n')

prob_topological_concordance = -1
prob_monophyletic_concordance = -1

outdir = '../../results/sim/'
outfilename = 'sim_out_' + tag + '.txt'
results_filename = 'sim_out_' + tag + '_summary.txt'

f_blocks_predicted = open(outdir + 'sim_out_' + tag + '_blocks_predicted.txt', 'w')
f_blocks_actual = open(outdir + 'sim_out_' + tag + '_blocks_actual.txt', 'w')
f_blocks_predicted_window = \
    open(outdir + 'sim_out_' + tag + '_blocks_predicted_window.txt', 'w')
f_blocks_actual_window = \
    open(outdir + 'sim_out_' + tag + '_blocks_actual_window.txt', 'w')

assert not theory, 'theory calculations not yet implemented for single cer population'
if theory:
    if theory_done:
        print 'reading theory results from file'
        f = open('out/' + results_filename, 'r')
        line = f.readline().split()
        assert(line[12] == 'prob_topological_concordance')
        assert(line[13] == 'prob_monophyletic_concordance')
        line = f.readline().split()
        prob_topological_concordance = float(line[12])
        prob_monophyletic_concordance = float(line[13])
        f.close()


f = open(outdir + outfilename, 'r')
line = f.readline()
n = 0
while line != '' and n < num_reps:
    if line == '//\n':
        print n

        # the main difference in adding in recombination is that we
        # now have to consider multiple trees, one for each block that
        # has not experienced recombination within it
        trees = []
        recomb_sites = []
        t_string = f.readline()
        # read in all the trees for blocks with no recombination within them
        while t_string[0] == '[':
            t_start = t_string.find(']') + 1
            recomb_sites.append(int(t_string[1:t_start-1]))
            t_string = t_string[t_start:-1]
            t = parse_tree(t_string)
            trees.append(t)
            t_string = f.readline()
        # read next couple of lines before sequences begin
        segsites = int(t_string[len('segsites: '):-1])
        #print 'segsites:', segsites
        positions = [float(x) for x in f.readline()[len('positions: '):].split()]
        #print 'positions:', positions
        #print recomb_sites
        # convert positions to integers
        positions = integer_positions(positions, num_sites)
        #print positions
        # read in sequences (at this point only sites that are polymorphic)
        seqs = []
        for i in range(num_samples):
            seqs.append(f.readline()[:-1])
            #print(seqs[-1])
            assert(len(seqs[-1]) > 0)

        # coding:
        # 0 -> matches cer ref but not par ref
        # 1 -> matches par ref but not cer ref
        # 2 -> matches cer ref and par ref
        # 3 -> doesn't match cer ref or par ref
        seqs_filled = []
        for seq in seqs:
            # nonpolymorphic sites will match both references
            s = ['2'] * num_sites
            # polymorphic sites
            for i in range(segsites):
                if seq[i] == seqs[ref_ind_cer][i]:
                    if seq[i] == seqs[ref_ind_par][i]:
                        s[positions[i]] = '2'
                    else:
                        s[positions[i]] = '0'
                elif seq[i] == seqs[ref_ind_par][i]:
                        s[positions[i]] = '1'
                else:
                        s[positions[i]] = '3'
                
            seqs_filled.append(s)

        ########
        # figure out which sites are actually introgressed by looking at trees for
        # all regions without recombination
        ########

        # store the sites that are introgressed for each cer strain
        introgressed_actual = [[] for i in range(num_samples_cer)]

        # sequence of states (cer or par), one for each site and strain (mostly like above)
        actual_state_seq = [[] for i in range(num_samples)]

        site_ind = 0
        # loop through the trees for all blocks with no recombination within them
        for ti in range(len(trees)):
            t = trees[ti]
            #print 'tree', ti, 'out of', len(trees)
            #print t
            
            # should this be 0 or 1 indexed??
            if is_concordant(t, range(1, num_samples_par), range(num_samples_par + 1, num_samples + 1)):
                concordant[n].append(1)
            else:
                concordant[n].append(0)
            
            # identify cerevisiae sequences that are actually introgressed from
            # paradoxus (based on coalescent tree); indexed from 0
            # only return for cer strains, not par
            introgressed = find_introgressed(t, t_par_cer, range(num_samples_par, num_samples))

            # number of lineages that were present when cerevisiae and paradoxus populations joined
            num_lineages_at_join[n].append(introgressed[1])

            # number of sites in the current block of sequence
            num_sites_t = recomb_sites[ti]

            # for all strains that have this block introgressed, at
            # the length of the block to the total number of
            # introgressed sites across all strains
            for i in introgressed[0]:
                num_introgressed_cer[n][i - num_samples_par] += num_sites_t
                introgressed_actual[i - num_samples_par] += range(site_ind, site_ind + num_sites_t)
            site_ind += num_sites_t

            for s in range(num_samples_par, num_samples):
                if s in introgressed[0]:
                    # 1 for par (introgressed)
                    actual_state_seq[s] += [1] * num_sites_t
                else:
                    # 0 for cer (not introgressed)
                    actual_state_seq[s] += [0] * num_sites_t

        ########
        # predict whether each site in each cer strain is
        # introgressed with a hidden markov model
        ########

        print 'HMM'

        predicted, hmm = predict_introgressed_hmm(seqs_filled[num_samples_par:])
        assert set(predicted[0]) == set([0]), predicted[0]
        if 1 in predicted[ref_ind_cer]:
            print 'PROBLEM\n' + 'HMM\n ' + str(predicted[ref_ind_cer - num_samples_par]) + '\n' + str(seqs_filled[ref_ind_cer])
        num_predicted_introgressed = [sum(x) for x in predicted]
        num_correct, num_introgressed_correct, actual_lens, predicted_lens, \
            num_predicted_tracts_actual, num_actual_tracts_predicted, \
            num_introgressed_tracts, num_not_introgressed_tracts, \
            num_predicted_introgressed_tracts, num_predicted_not_introgressed_tracts = \
            evaluate_predicted(predicted, actual_state_seq[num_samples_par:])
        blocks_predicted, blocks_actual = \
            evaluate_predicted_blocks(predicted, actual_state_seq[num_samples_par:])

        ########
        # predict whether each site in each cer strain is
        # introgressed by windowing (each window is in each strain
        # is either completely introgressed or not introgressed)
        ########

        print 'windowing'

        predicted = predict_introgressed(seqs_filled[num_samples_par:], window_size, window_shift)
        if 1 in predicted[ref_ind_cer]:
            print 'PROBLEM\n' + 'window\n ' + str(predicted[ref_ind_cer - num_samples_par]) + '\n' + str(seqs_filled[ref_ind_cer])
        num_predicted_introgressed_window = [sum(x) for x in predicted]
        num_correct_window, num_introgressed_correct_window, actual_lens_window, predicted_lens_window, \
            num_predicted_tracts_actual_window, num_actual_tracts_predicted_window, \
            num_introgressed_tracts_window, num_not_introgressed_tracts_window, \
            num_predicted_introgressed_tracts_window, num_predicted_not_introgressed_tracts_window = \
            evaluate_predicted(predicted, actual_state_seq[num_samples_par:])
        blocks_predicted_window, blocks_actual_window = \
            evaluate_predicted_blocks(predicted, actual_state_seq[num_samples_par:])

        ########
        # sequence identities
        ########

        total_len = 0

        for seq in seqs_filled[num_samples_par:]:

            count0 = seq.count('0')
            count1 = seq.count('1')
            count2 = seq.count('2')
            count3 = len(seq) - count0 - count1 - count2

            avg_id_cer[n] += count0 + count2
            avg_id_cer_par[n] += count1 + count2
            total_len += len(seq)
 
        avg_id_cer[n] /= float(total_len)
        avg_id_cer_par[n] /= float(total_len)

        #####
        # write results to file
        #####

        # LIST, one for each individual
        fout.write(str(num_introgressed_cer[n][:num_samples_cer]) + '\t')

        # LIST, one for each individual
        # introgressed tracts (positives) - just take an arbitrary one
        # for each reference because they're the same
        fout.write(str(num_introgressed_tracts) + '\t')

        # LIST, one for each individual
        # not introgressed tracts (negatives)
        fout.write(str(num_not_introgressed_tracts) + '\t')

        # LIST, one for each individual
        fout.write(str(num_predicted_introgressed) + '\t')

        # LIST, one for each individual
        fout.write(str(num_predicted_introgressed_tracts) + '\t')

        # LIST, one for each individual
        fout.write(str(num_predicted_not_introgressed_tracts) + '\t')

        # LIST, one for each individual
        fout.write(str(num_predicted_introgressed_window) + '\t')

        # LIST, one for each individual
        fout.write(str(num_predicted_introgressed_tracts_window) + '\t')

        # LIST, one for each individual
        fout.write(str(num_predicted_not_introgressed_tracts_window) + '\t')

        # LIST, one for each individual
        fout.write(str(num_correct) + '\t')
        fout.write(str(num_introgressed_correct) + '\t')
        fout.write(str(num_predicted_tracts_actual) + '\t')
        fout.write(str(num_actual_tracts_predicted) + '\t')
        fout.write(str(num_correct_window) + '\t')
        fout.write(str(num_introgressed_correct_window) + '\t')
        fout.write(str(num_predicted_tracts_actual_window) + '\t')
        fout.write(str(num_actual_tracts_predicted_window) + '\t')

        # LIST, one for each tract across all individuals
        fout.write(str(actual_lens) + '\t')
        fout.write(str(predicted_lens) + '\t')
        fout.write(str(predicted_lens_window) + '\t')

        fout.write(str(avg_id_cer[n]) + '\t')
        fout.write(str(avg_id_cer_par[n]) + '\t')

        fout.write(str(mean(num_lineages_at_join[n])) + '\t')

        fout.write(str(mean(concordant[n])) + '\t')

        fout.write(str(hmm.init[0]) + '\t')
        fout.write(str(hmm.init[1]) + '\t')

        fout.write(str(hmm.trans[0][0]) + '\t')
        fout.write(str(hmm.trans[0][1]) + '\t')
        fout.write(str(hmm.trans[1][0]) + '\t')
        fout.write(str(hmm.trans[1][1]) + '\t')

        fout.write(str(hmm.emis[0]['0']) + '\t')
        fout.write(str(hmm.emis[0]['1']) + '\t')
        fout.write(str(hmm.emis[0]['2']) + '\t')
        fout.write(str(hmm.emis[1]['0']) + '\t')
        fout.write(str(hmm.emis[1]['1']) + '\t')
        fout.write(str(hmm.emis[1]['2']) + '\n')

        fout.flush()

        sys.stdout.flush()

        for b in blocks_predicted:
            f_blocks_predicted.write(str(n) + ' ' +  ' '.join([str(x) for x in b]) + '\n')
        for b in blocks_actual:
            f_blocks_actual.write(str(n) + ' ' + ' '.join([str(x) for x in b]) + '\n')
        for b in blocks_predicted_window:
            f_blocks_predicted_window.write(str(n) + ' ' + ' '.join([str(x) for x in b]) + '\n')
        for b in blocks_actual_window:
            f_blocks_actual_window.write(str(n) + ' ' + ' '.join([str(x) for x in b]) + '\n')

        f_blocks_predicted.flush()
        f_blocks_actual.flush()
        f_blocks_predicted_window.flush()
        f_blocks_actual_window.flush()

        n += 1

    line = f.readline()


f.close()
fout.close()
f_blocks_predicted.close()
f_blocks_actual.close()
f_blocks_predicted_window.close()
f_blocks_actual_window.close()


