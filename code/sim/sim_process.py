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

def read_one_sim(f, num_sites, num_samples):

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

    return [trees, recomb_sites, segsites, positions, seqs]
