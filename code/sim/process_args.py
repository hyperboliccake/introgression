import sys

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

def process_args(print_args=True):

    # store all arguments in dictionary
    d = {}

    # index of current argument being processed
    i = 1

    ########
    # arguments!
    ########

    # simulation tag
    d['tag'] = sys.argv[i]
    i += 1

    # topology of the species
    d['topology'] = sys.argv[i]
    i += 1

    # species names
    d['species'] = get_labels(parse_topology(d['topology']))
    assert len(species) == 2 or len(species) == 3, species

    # expected length and number of tracts...
    expected_tract_lengths = {}
    expected_num_tracts = {}

    # ...for the species with introgression
    d['species_to'] = sys.argv[i]
    i += 1
    assert d['species_to'] in d['species']
    d['num_samples_species_to'] = int(sys.argv[i])
    i += 1
    d['N0_species_to'] = int(sys.argv[i])
    i += 1
    assert sys.argv[i] == 'to'
    i += 1

    # ...for one species introgression is coming from
    d['species_from1'] = sys.argv[i]
    i += 1
    assert d['species_from1'] in d['species']
    d['num_samples_species_from1'] = int(sys.argv[i])
    i += 1
    d['N0_species_from1'] = int(sys.argv[i])
    i += 1
    d['migration_from1'] = float(sys.argv[i]) * 2 * d['N0_species_from1']
    i += 1
    expected_tract_lengths[d['species_from1']] = float(sys.argv[i])
    i += 1
    expected_num_tracts[d['species_from1']] = int(sys.argv[i])
    i += 1
    d['has_ref_from1'] = (sys.argv[i] == 'ref')
    i += 1

    # ...for second species introgression is coming from (optional)
    d['species_from2'] = None
    d['num_samples_species_from2'] = 0
    d'N0_species_from2'] = d['N0_species_from1']
    d['migration_from2'] = 0
    d['has_ref_from2'] = False
    if len(d['species']) == 3:
        d['species_from2'] = sys.argv[i]
        i += 1
        assert d['species_from2'] in d['species']
        d['num_samples_species_from2'] = int(sys.argv[i])
        i += 1
        d'[N0_species_from2'] = int(sys.argv[i])
        i += 1
        d['migration_from2'] = float(sys.argv[i]) * 2 * d['N0_species_from2']
        i += 1
        expected_tract_lengths[d['species_from2']] = float(sys.argv[i])
        i += 1
        expected_num_tracts[d['species_from2']] = int(sys.argv[i])
        i += 1
        d['has_ref_from2'] = (sys.argv[i] == 'ref')
        i += 1

    # only makes sense to have one unknown species at most
    assert d['has_ref_from1'] or d['has_ref_from2']

    # calculate these based on remaining bases
    expected_num_tracts[d['species_to']] = sum(expected_num_tracts.values()) + 1
    expected_num_introgressed_bases = expected_tract_lengths[d['species_from1']] * \
        expected_num_tracts[d['species_from1']]
    if d['species_from2'] != None:
        expected_num_introgressed_bases += expected_tract_lengths[d['species_from2']] * \
            expected_num_tracts[d['species_from2']]
    expected_tract_lengths[d['species_to']] = float(expected_num_introgressed_bases) / \
        expected_num_tracts[d['species_to']]

    assert d'N0_species_to'] == d['N0_species_from1'] and \
        d['N0_species_to'] == d['N0_species_from2']

    d['topology'] = parse_topology(d['topology'], 1/float(2 * d['N0_species_to']))

    # 13,500 sites to get about 10% with one recombination event, .3% with
    # more than one (based on poisson(.1), 1 recombination per chromosome
    # of average length 750,000)
    d['num_sites'] = int(sys.argv[i])
    i += 1

    # parameter is recombination rate between adjacent bp per generation
    # should probably be 1/750000 + 6.1 * 10^-6 (where 750000 is average
    # chr size)
    # recombination rate
    d['rho'] = 2 * f['N0_species_to'] * float(sys.argv[i]) * (d['num_sites'] - 1)
    i += 1

    d['outcross_rate'] = float(sys.argv[i])
    i += 1
    
    d['rho'] *= d['outcross_rate']

    d['theta'] = gp.mu * 2 * d['num_sites'] * d['N0_species_to']

    d['num_reps'] = int(sys.argv[i])

    #####
    # reference stuff
    #####

    d['num_samples'] = d['num_samples_species_to'] + \
        d['num_samples_species_from1'] + d['num_samples_species_from2']

    # species_to always comes first
    index_to_species = [d['species_to']] * d['num_samples_species_to'] + \
        [d['species_from1']] * d['num_samples_species_from1'] + \
        [d['species_from2']] * d['num_samples_species_from2']

    species_to_indices = {}
    for i in index_to_species:
        species = index_to_species[i]
        if species not in species_to_indices:
            species_to_indices[species] = []
        species_to_indices[species].append(i)
    d['species_to_indices'] = species_to_indices

    # take first index from each population to be reference sequence
    ref_ind_species_to = 0
    ref_ind_species_from1 = num_samples_species_to
    ref_ind_species_from2 = num_samples_species_to + num_samples_species_from1
    ref_inds = [ref_ind_species_to]
    states = [species_to, species_from1]
    unknown_species = None
    if has_ref_from1:
        ref_inds.append(ref_ind_species_from1)
    else:
        unknown_species = species_from1
    if species_from2 != None:
        states.append(species_from2)
        if has_ref_from2:
            ref_inds.append(ref_ind_species_from2)
        else:
            unknown_species = species_from2

    d['unknown_species'] = unknown_species
    d['states'] = states
    d['ref_inds'] = ref_inds
    # if there are three species and the second is the species that has no
    # reference, flip the order of the states so that the species with no
    # reference always comes last; this will ensure that the indices of
    # the species in states correspond to the indices of the references
    # (and the sequence codings later); ACTUALLY just force the unknown
    # species to come last
    if species_from2 != None:
        assert has_ref_from1

    if print_args:
        for key in d.keys():
            print key, d[key]

    return d
    
    #return tag, topology, species_to, species_from1, species_from2, \
    #    num_samples_species_to, num_samples_species_from1, num_samples_species_from2, \
    #    N0_species_to, N0_species_from1, N0_species_from2, \
    #    migration_from1, migration_from2, \
    #    expected_tract_lengths, \
    #    expected_num_tracts, \
    #    has_ref_from1, has_ref_from2, \
    #    rho, outcross_rate, theta, num_sites, num_reps



