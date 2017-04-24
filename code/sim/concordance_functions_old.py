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
def is_partial_clade(t, species, index_to_species):
    labels = get_labels(t)
    for l in labels:
        if index_to_species[l] != species:
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

# actual introgressed with 2 species total
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

