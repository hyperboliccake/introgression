import sys
import concordance_functions
sys.path.append('..')
import global_params as gp

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

def process_args(arg_list, i=1, print_args=True):

    print arg_list

    # store all arguments in dictionary
    d = {}

    ########
    # arguments!
    ########

    # simulation tag
    d['tag'] = arg_list[i]
    i += 1

    # topology of the species
    d['topology'] = arg_list[i]
    i += 1

    # species names
    d['species'] = concordance_functions.get_labels(parse_topology(d['topology']))
    print d['topology']
    print d['species']
    assert len(d['species']) == 2 or len(d['species']) == 3, d['species']

    # ...for the species with introgression
    d['species_to'] = arg_list[i]
    i += 1
    assert d['species_to'] in d['species']
    d['num_samples_species_to'] = int(arg_list[i])
    i += 1
    d['N0_species_to'] = int(arg_list[i])
    i += 1

    # ...for one species introgression is coming from
    d['species_from1'] = arg_list[i]
    i += 1
    assert d['species_from1'] in d['species']
    d['num_samples_species_from1'] = int(arg_list[i])
    i += 1
    d['N0_species_from1'] = int(arg_list[i])
    i += 1
    d['migration_from1'] = float(arg_list[i]) * 4 * d['N0_species_from1']
    i += 1

    # ...for second species introgression is coming from (optional)
    d['species_from2'] = None
    d['num_samples_species_from2'] = 0
    d['N0_species_from2'] = d['N0_species_from1']
    d['migration_from2'] = 0

    if len(d['species']) == 3:
        d['species_from2'] = arg_list[i]
        i += 1
        assert d['species_from2'] in d['species']
        d['num_samples_species_from2'] = int(arg_list[i])
        i += 1
        d['N0_species_from2'] = int(arg_list[i])
        i += 1
        d['migration_from2'] = float(arg_list[i]) * 4 * d['N0_species_from2']
        i += 1

    # only supporting same population size for all species right now,
    # but could change this later
    assert d['N0_species_to'] == d['N0_species_from1'] and \
        d['N0_species_to'] == d['N0_species_from2']

    d['topology'] = parse_topology(d['topology'], 1/float(4 * d['N0_species_to']))

    # 13,500 sites to get about 10% with one recombination event, .3% with
    # more than one (based on poisson(.1), 1 recombination per chromosome
    # of average length 750,000)
    d['num_sites'] = int(arg_list[i])
    i += 1

    # parameter is recombination rate between adjacent bp per
    # generation should probably be 1/750000 + 6.1 * 10^-6 = 7.425 *
    # 10^-6 (where 750000 is average chr size) recombination rate
    d['rho'] = 4 * d['N0_species_to'] * float(arg_list[i]) * (d['num_sites'] - 1)
    i += 1

    d['outcross_rate'] = float(arg_list[i])
    i += 1
    
    d['rho'] *= d['outcross_rate']

    d['theta'] = gp.mu * 4 * d['num_sites'] * d['N0_species_to']

    d['num_reps'] = int(arg_list[i])

    #####
    # reference stuff
    #####

    d['num_samples'] = d['num_samples_species_to'] + \
        d['num_samples_species_from1'] + d['num_samples_species_from2']

    # species_to always comes first
    d['index_to_species'] = [d['species_to']] * d['num_samples_species_to'] + \
        [d['species_from1']] * d['num_samples_species_from1'] + \
        [d['species_from2']] * d['num_samples_species_from2']

    species_to_indices = {}
    for ind in range(len(d['index_to_species'])):
        species = d['index_to_species'][ind]
        if species not in species_to_indices:
            species_to_indices[species] = []
        species_to_indices[species].append(ind)
    d['species_to_indices'] = species_to_indices

    if print_args:
        for key in d.keys():
            print key, d[key]

    return d, i
    
def process_args_by_tag(fn, tag):
    
    f = open(fn, 'r')
    line = f.readline()
    args = None
    while line != '':
        args, last_read = process_args(line.strip().split(' '), i=0)
        if args['tag'] == tag:
            break
        line = f.readline()
    f.close()
    return args
