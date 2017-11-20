import sys
import os
import copy
import concordance_functions

def seq_id(a, b, l = -1, use_gaps = False):
    assert len(a) == len(b)
    ndiff = 0
    ntotal = 0
    for i in xrange(len(a)):
        # use sequence a as denominator
        if a[i] != '-':
            if b[i] != '-' or use_gaps:
                ntotal += 1
                if a[i] != b[i]:
                    ndiff += 1
    if l == -1:
        if ntotal == 0:
            return -1, 0
        return 1 - float(ndiff) / ntotal, ntotal
    return 1 - float(ndiff) / l

def sim_stats(sim, args):

    # give some summary stats about simulated sequences, such as:
    # - average sequence identity within and between species

    num_species = len(args['species'])
    ids = {}
    for i in range(num_species):
        for j in range(i, num_species):
            ids_ij = []
            species1 = args['species'][i]
            species2 = species1
            inds1 = args['species_to_indices'][species1]
            # within species
            if i == j:
                for k in range(len(inds1)):
                    for l in range(k+1, len(inds1)):
                        s1 = inds1[k]
                        s2 = inds1[l]
                        ids_ij.append(\
                            seq_id(sim['seqs'][s1], sim['seqs'][s2], args['num_sites']))
                        
            # between species
            else:
                species2 = args['species'][j]
                inds2 = args['species_to_indices'][species2] 
                for s1 in inds1:
                    for s2 in inds2:
                        ids_ij.append(\
                            seq_id(sim['seqs'][s1], sim['seqs'][s2], args['num_sites']))

            ids[(species1, species2)] = ids_ij

    stats = {'seq_ids':ids}

    return stats

def weighted_average(values, weights):
    assert len(values) == len(weights)

    total = 0
    for i in range(len(values)):
        total += values[i] * weights[i]
    total = float(total) / sum(weights)

    return total

def calculate_ils(sim, args):

    num_trees = len(sim['trees'])
    concordant = [0] * num_trees
    for ti in range(num_trees):
        t = sim['trees'][ti]

        if concordance_functions.is_concordant(\
            t, args['index_to_species'], args['species_to']):
            concordant[ti] = 1

    stats = {'concordant_site_freq':weighted_average(concordant, sim['recomb_sites']),\
                 'concordant_tree_frac':float(sum(concordant)) / num_trees}
            
    return stats


def find_introgressed_2(t, cutoff_time, species_to, index_to_species):
    # find species that are introgressed within a single unrecombined
    # block, given there is one species migrating into one other
    # species

    # strategy: find all lineages that exist at the time the two
    # populations join; then for each of those lineages, if any are
    # mixed-species, then there's introgression

    # the subtrees that exist at the time the populations join
    lineages = concordance_functions.split(t, cutoff_time)
    # make this a dictionary instead of just a list of the
    # introgressed individuals because in other cases we'll want to
    # keep track of which species introgression came from
    introgressed = {}
    num_lineages_species_to = 0
    for l in lineages:
        not_introgressed, s = \
            concordance_functions.is_partial_clade(l, species_to, index_to_species)
        # if this lineage has only species_to members, then add this
        # to the number of lineages for the species
        if not_introgressed:
            num_lineages_species_to += 1
        # if there's any occurence of species_from in this clade, then
        # mark all individuals as coming from species_from (since
        # we're only allowing migration in one direction)
        else:
            for label in concordance_functions.get_labels(l):
                # note that species labels/indices start at 0 (see
                # sim_process.parse_ms_tree_helper())
                if index_to_species[label] == species_to:
                    introgressed[label] = s

    return introgressed, num_lineages_species_to

def find_introgressed_3(t, topology, species_to, index_to_species):
    # find species that are introgressed within a single unrecombined
    # block, given there are two species migrating into one other
    # species

    # IMPORTANT: this function assumes migration only happens after
    # most recent divergence

    # strategy: look at lineages that exist at most recent join time;
    # if any are not all one species, then mark all members as the
    # species that's not the to species (there's no migration between
    # the two from species); then look at the more distant join time;
    # now we can't just check if the lineages are all one species
    # because two of them have joined; so instead...figure out which
    # species is the last to join (from the topology), then...  nvm,
    # only allowing migration later on makes this less complicated
    
    join_time_species_to = topology[2]
    subtree = topology[0]
    last_to_join = topology[1]
    if type(subtree) != type([]):
        subtree = topology[1]
        last_to_join = topology[0]
    if last_to_join != species_to:
        join_time_species_to = subtree[2]

    return find_introgressed_2(t, join_time_species_to, species_to, index_to_species)

def find_introgressed(sim, args):

    ##======
    # figure out which sites are actually introgressed by separately
    # looking at the tree for each stretch without recombination
    ##======

    # indices of all species_to individuals (the ones for which we
    # care about introgression)
    inds = args['species_to_indices'][args['species_to']]
    # sequence of states, one for each site and strain
    actual_state_seq = \
        dict(zip(inds, [[] for i in range(args['num_samples_species_to'])]))
    # how many bases are introgressed in total in each strain
    num_introgressed = \
        dict(zip(inds, [0] * args['num_samples_species_to']))
    # how many bases are introgressed in total in each strain but NOT
    # in the reference (i.e. the bases we'd have a chance of detecting)
    num_introgressed_non_ref = \
        dict(zip(inds, [0] * args['num_samples_species_to']))
    # loop through the trees for all blocks with no recombination
    # within them
    num_trees = len(sim['trees'])
    for ti in range(num_trees):

        # note that species indices/labels are shifted to start at
        # 0 instead of 1
        t = sim['trees'][ti]
        # identify sequences that are introgressed from the one or
        # two other species, based on coalescent tree; could clean
        # this up a little
        introgressed = None
        # two species total
        if args['species_from2'] == None:
            introgressed, num_lineages = \
                find_introgressed_2(t, args['topology'][2], args['species_to'], \
                                        args['index_to_species'])
        # three species total
        else:
            introgressed, num_lineages = \
                find_introgressed_3(t, args['topology'], args['species_to'], \
                                        args['index_to_species'])

        # number of sites in the current block of sequence
        num_sites_t = sim['recomb_sites'][ti]
        # for all strains that have this block introgressed, add
        # the length of the block to the total number of
        # introgressed sites across all strains; also update the
        # state sequence
        
        # just call ref ind the first index of the species
        ref_ind = min(args['species_to_indices'][args['species_to']]) 
        for i in inds:
            if introgressed.has_key(i):
                num_introgressed[i] += num_sites_t
                actual_state_seq[i] += [introgressed[i]] * num_sites_t
                if not introgressed.has_key(ref_ind):
                    num_introgressed_non_ref[i] += num_sites_t
            else:
                actual_state_seq[i] += [args['species_to']] * num_sites_t

    stats = {'num_introgressed':num_introgressed, \
             'num_introgressed_non_ref':num_introgressed_non_ref}

    return stats, actual_state_seq

def one_output_header_chunk(d, sep):

    s = ''
    for key in sorted(d.keys()):
        value = d[key]
        if type(value) == type({}):
            for subkey in sorted(value.keys()):
                s += str(key) + '_' + str(subkey) + sep
        elif type(value) == type([]):
            for i in range(len(value)):
                s += str(key) + '_' + str(i) + sep
        else:
            s += str(key) + sep
    return s

def write_output_headers(summary_info, concordance_info, introgression_info, f, sep):

    line_string = ''

    line_string += one_output_header_chunk(summary_info, sep)
    line_string += one_output_header_chunk(concordance_info, sep)
    line_string += one_output_header_chunk(introgression_info, sep)
    
    f.write(line_string[:-1] + '\n')
    

def one_output_chunk(d, sep):

    s = ''
    for key in sorted(d.keys()):
        value = d[key]
        if type(value) == type({}):
            for subkey in sorted(value.keys()):
                s += str(value[subkey]) + sep
        elif type(value) == type([]):
            for i in range(len(value)):
                s += str(value[i]) + sep
        else:
            s += str(value) + sep
    return s

def write_output_line(summary_info, concordance_info, introgression_info, f, \
                          header = False):

    sep = '\t'

    if header:
        write_output_headers(summary_info, concordance_info, introgression_info, f, sep)

    line_string = ''

    line_string += one_output_chunk(summary_info, sep)
    line_string += one_output_chunk(concordance_info, sep)
    line_string += one_output_chunk(introgression_info, sep)
    
    f.write(line_string[:-1] + '\n')

    
