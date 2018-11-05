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

def read_one_sim(f, num_sites, num_samples):

    sim = {}

    line = f.readline()
    while line != '//\n':
        if line == '':
            return None
        line = f.readline()

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
    sim['recomb_sites'] = recomb_sites
    sim['trees'] = trees

    # read next couple of lines before sequences begin
    sim['segsites'] = int(t_string[len('segsites: '):-1])
    positions = [float(x) for x in f.readline()[len('positions: '):].split()]
    # convert positions to integers
    # (allow sites to be hit multiple times because that seems reasonable)
    # (zero-indexed)
    sim['positions'] = integer_positions(positions, num_sites, allow_multi_hit=True)
    # read in sequences (at this point only sites that are polymorphic)
    seqs = []
    for i in range(num_samples):
        seqs.append(f.readline()[:-1])
        assert(len(seqs[-1]) > 0)
    sim['seqs'] = seqs

    return sim

def convert_to_blocks_one(state_seq, states):
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
    # make first block extend to first position, and last to last
    # (for the case that we're only considering a subset of sites)
    #first_block = blocks[state_seq[0]][0]
    #blocks[state_seq[0]][0] = (seq_start, first_block[1])
    return blocks


def convert_to_blocks(state_seqs, states):
    blocks = {}
    for ind in state_seqs.keys():
        blocks[ind] = convert_to_blocks_one(state_seqs[ind], states)
    return blocks

def write_introgression(state_seq, f, rep, states):
    # file format is:
    # rep 0
    # 2\tpar:3,4,5,9,10,12\tbay:15,16
    # 3\tpar:3,4,5,9,10,12\tbay:15,16
    # rep 1 ...

    d = {}
    # loop through all individuals (in species_to)
    for ind in state_seq.keys():
        # loop through all positions in the sequence
        d[ind] = {}
        for p in range(len(state_seq[ind])):
            current_species = state_seq[ind][p]
            if current_species != species:
                if not d[ind].has_key(current_species):
                    d[ind][current_species] = []
                d[ind][current_species].append(p)

    f.write('rep ' + str(rep) + '\n')

    for ind in d.keys():
        f.write(str(ind))
        for s in d[ind]:
            f.write('\t' + s + ':')
            f.write(','.join([str(x) for x in d[ind][s]]))
        f.write('\n')
    f.flush()

def read_introgression(f, line, states):
    # inverse of write_introgression function above
    
    d = {} # keyed by individual, then strain, list of sites
    assert line.startswith('rep'), line
    rep = int(line[len('rep '):-1])
    # process each individual
    line = f.readline()
    while line != '' and not line.startswith('rep'):
        x = line[:-1].split('\t')
        ind = line[0]
        d_ind = {}
        for species in states:
            d_ind[species] = []
        for s in x[1:]:
            species, sites = s.split(':')
            d_ind[species] = [int(i) for i in sites.split(',')]
        d[ind] = d_ind
        line = f.readline()
                
    return d, rep, line

def write_introgression_blocks(state_seq_blocks, f, rep, states):
    # file format is:
    # rep 0
    # 2\tcer:0-2,6-8\tpar:3-5,9-12\tbay:15-16
    # 3\tcer:8-9,16-16\tpar:3-7,10-12\tbay:15-15
    # rep 1 ...

    f.write('rep ' + str(rep) + '\n')

    for ind in state_seq_blocks:
        f.write(str(ind))
        for state in states:
            f.write('\t' + state + ':')
            blocks_string = ''
            for block in state_seq_blocks[ind][state]:
                blocks_string += str(block[0]) + '-' + str(block[1]) + ','
            f.write(blocks_string[:-1])
        f.write('\n')
    f.flush()

def read_introgression_blocks(f, line, states):
    # inverse of write_introgression_blocks function above
    
    d = {} # keyed by individual, then species, list of (start, end)
    assert line.startswith('rep'), line
    rep = int(line[len('rep '):-1])
    # process each individual
    line = f.readline()
    while line != '' and not line.startswith('rep'):
        x = line[:-1].split('\t')
        ind = int(line[0])
        d_ind = {}
        for state in states:
            d_ind[state] = []
        for s in x[1:]:
            species, blocks = s.split(':')
            blocks = blocks.split(',')
            blocks = filter(lambda x: x != '', blocks) # a lil janky
            for block in blocks:
                block_start, block_end = block.split('-')
                d_ind[species].append((int(block_start), int(block_end)))
        d[ind] = d_ind
        line = f.readline()
                
    return d, rep, line

def unblock(blocks, num_sites):
    # blocks is keyed by individual, then species, list of (start, end)
    d = {} # keyed by individual, list of species, one for each site
    for ind in blocks.keys():
        d[ind] = ['None' for x in range(num_sites)]
        for state in blocks[ind].keys():
            for block in blocks[ind][state]:
                for i in range(block[0], block[1]+1):
                    d[ind][i] = state
    return d


def write_state_probs(probs, f, rep):

    # probs is keyed by individual, list of sites, each site dic keyed
    # by state
    
    # file format is:
    # rep 0
    # 2\tcer:.1,.2,.3\tpar:.9,.8,.7
    # 3\tcer:.4,.2,.3\tpar:.6,.8,.7
    # rep 1 ...

    f.write('rep ' + str(rep) + '\n')

    for ind in probs.keys():
        f.write(str(ind))
        for state in probs[ind][0].keys():
            f.write('\t' + state + ':')
            probs_string = ','.join(["{0:.5f}".format(site[state]) \
                                     for site in probs[ind]])
            f.write(probs_string)
        f.write('\n')
    f.flush()

def read_state_probs(f, line):

    d = {} # keyed by individual, then species, list of probs
    assert line.startswith('rep'), line
    rep = int(line[len('rep '):-1])
    # process each individual
    line = f.readline()
    while line != '' and not line.startswith('rep'):
        x = line[:-1].split('\t')
        ind = int(line[0])
        d_ind = {}
        for s in x[1:]:
            species, probs = s.split(':')
            probs = [float(x) for x in probs.split(',')]
            d_ind[species] = probs
        d[ind] = d_ind
        line = f.readline()
                
    return d, rep, line

def threshold_predicted(predicted, probs, threshold, default_state):

    predicted_thresholded = []
    for i in range(len(predicted)):
        if probs[i] >= threshold:
            predicted_thresholded.append(predicted[i])
        else:
            predicted_thresholded.append(default_state)
    return predicted_thresholded

def fill_seq(seq, polymorphic_sites, nsites, fill):
    s = [fill for x in range(nsites)]
    for i in range(len(polymorphic_sites)):
        s[polymorphic_sites[i]] = seq[i]
    return s

# add in the nonpolymorphic sites
def fill_seqs(polymorphic_seqs, polymorphic_sites, nsites, fill):
    
    # note that polymorphic sites can have duplicates
    seqs_filled = []
    for seq in polymorphic_seqs:
        s = fill_seq(seq, polymorphic_sites, nsites, fill)
        seqs_filled.append(''.join(s)) # TODO should return this as list instead
    return seqs_filled

def get_max_path(p):
    # p is a list of dictionaries, one per site; each dict has keys
    # for each state, with associated probability
    max_path = []
    max_probs = []
    for site_probs in p:
        max_state = None
        max_prob = -1
        for state in site_probs:
            if site_probs[state] > max_prob:
                max_prob = site_probs[state]
                max_state = state
        max_path.append(max_state)
        max_probs.append(max_prob)
    return max_path, max_probs

