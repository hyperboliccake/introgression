import sys
import os
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

    num_species = len(args['states'])
    ids = {}
    for i in range(num_species):
        for j in range(i, num_species):
            ids_ij = []
            species1 = args['states'][i]
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
                species2 = args['states'][j]
                inds2 = args['species_to_indices'][species2] 
                for s1 in inds1:
                    for s2 in inds2:
                        ids_ij.append(\
                            seq_id(sim['seqs'][s1], sim['seqs'][s2], args['num_sites']))

            ids[(species1, species2)] = ids_ij

    return ids

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

    return weighted_average(concordant, sim['recomb_sites']), \
        float(sum(concordant)) / num_trees

def find_introgressed(sim, args):

    ##======
    # figure out which sites are actually introgressed by separately
    # looking at the tree for each stretch without recombination
    ##======

    # sequence of states, one for each site and strain
    actual_state_seq = [[] for i in range(args['num_samples_species_to'])]
    # how many bases are introgressed in total in each strain
    num_introgressed = [0] * args['num_samples_species_to']
    # loop through the trees for all blocks with no recombination
    # within them
    num_trees = len(sim['trees'])
    for ti in range(num_trees):

        # note that species indices/labels are shifted to start at
        # 0 instead of 1
        t = trees[ti]

        # identify sequences that are introgressed from the one or
        # two other species, based on coalescent tree; could clean
        # this up a little
        introgressed = None

        # two species
        if sim['species_from2'] == None:
            find_introgressed_2(sim, args)
        # three species
        else:
            find_introgressed_3(sim_args)


        print introgressed[0]

        # number of sites in the current block of sequence
        num_sites_t = recomb_sites[ti]

        # for all strains that have this block introgressed, add
        # the length of the block to the total number of
        # introgressed sites across all strains; also update the
        # state sequence
        for i in range(num_samples_species_to):
            if introgressed[i] != species_to:
                num_introgressed[i] += num_sites_t
            actual_state_seq[i] += [introgressed[i]] * num_sites_t


