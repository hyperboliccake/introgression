import sys
import os

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
    for i in range(num_species):
        for j in range(i, num_species):
            species1 = args['states'][i]
            inds1 = args['species_to_indices'][species1]
            if i == j:
                for k in range(len(inds1)):
                    for l in range(k+1, len(inds1)):
                        seq_id(sim['seqs'][k], sim['seqs'][l], args['num_sites'])
            else:
                species2 = args['states'][j]
                inds2 = args['species_to_indices'][species2]
            


def find_introgressed(sim, ):

    ##======
    # figure out which sites are actually introgressed by
    # separately looking at the tree for each stretch without
    # recombination
    ##======

    # sequence of states, one for each site and strain
    actual_state_seq = [[] for i in range(num_samples_species_to)]
    # keep track of whether each tree is concordant with the species tree
    concordant = []
    # and how many lineages of the to species are left when it
    # first joins another species
    num_lineages_at_join = []
    # and how many bases are introgressed in total in each strain
    num_introgressed = [0] * num_samples_species_to
    # loop through the trees for all blocks with no recombination
    # within them
    for ti in range(len(trees)):

        # note that species indices/labels are shifted to start at
        # 0 instead of 1
        t = trees[ti]
            
        # is this tree concordant with the species tree? (only
        # checks whether the to species is monophyletic, which
        # indicates that ILS not possible)
        if is_concordant(t, index_to_species, species_to):
            concordant.append(True)
        else:
            concordant.append(False)

        # identify sequences that are introgressed from the one or
        # two other species, based on coalescent tree; could clean
        # this up a little
        introgressed = None
        num_lineages_at_join_current = None
        # three species
        if num_from_species == 2:
            introgressed, num_lineages_at_join_current = \
                find_introgressed_3(t, species_to, topology, index_to_species)
        # two species
        else:
            assert num_from_species == 1
            # introgressed is a list of species (one entry for
            # each individual in to species)
            introgressed, num_lineages_at_join_current = \
                find_introgressed_2(t, topology[2], species_to, index_to_species)
        print introgressed[0]
        # number of lineages that were present when all
        # populations joined
        num_lineages_at_join.append(num_lineages_at_join_current)

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

