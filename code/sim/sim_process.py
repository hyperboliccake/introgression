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

def write_introgression(state_seq, f, rep, species):

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

                

