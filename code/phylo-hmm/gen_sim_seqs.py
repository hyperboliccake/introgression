# THIS IS BS
##############
# based on ms output, make file with one tree per site
# generate sequence by jukes cantor model with seq-gen


import os
import sys
sys.path.insert(0, '../sim')
import sim_analyze_hmm_bw
sys.path.insert(0, '..')
import global_params as gp


#####
# read parameters
#####

tag, topology, species_to, species_from1, species_from2, \
    num_samples_species_to, num_samples_species_from1, num_samples_species_from2, \
    N0_species_to, N0_species_from1, N0_species_from2, \
    migration_from1, migration_from2, \
    expected_length_introgressed, \
    expected_num_introgressed_tracts, \
    has_ref_from1, has_ref_from2, \
    rho, outcross_rate, theta, num_sites, num_reps = \
    sim_analyze_hmm_bw.process_args(sys.argv)

num_samples = num_samples_species_to + num_samples_species_from1 + num_samples_species_from2

# TODO don't hardcode this
labels = ['1', '2', '3', '4']
names = ['C1', 'C2', 'P','OUTGROUP']
label_to_name = dict(zip(labels, names))

#####
# read appropriate ms output file
#####

gp_dir = '../'
ms_out_fn = gp_dir + gp.sim_out_dir + '/ms/' + gp.sim_out_prefix + tag + gp.sim_out_suffix
f = open(ms_out_fn, 'r')

# format is //\n[num sites]tree\n[num sites]tree\n etc
line = f.readline()
trees = []
while line != '':
    # start next rep block
    if line == '//\n':
        trees_rep = []
        line = f.readline()
        while line[0] == '[':
            close_bracket = line.find(']')
            n = int(line[1:close_bracket])
            tree = line[close_bracket+1:-1]
            # TODO make this less dumb
            for label in labels:
                tree = tree.replace('(' + label + ':', '(' + label_to_name[label] + ':')
                tree = tree.replace(',' + label + ':', ',' + label_to_name[label] + ':')
            trees_rep.append((n, tree))
            line = f.readline()
        trees.append(trees_rep)
    line = f.readline()
f.close()
assert len(trees) == num_reps, str(len(trees)) + ' ' + str(num_reps)

#####
# make file with one tree per site
#####

tree_fns = []
for rep in range(num_reps):
    individual_trees_fn = gp_dir + gp.sim_out_dir + '/seq-gen/temp/' + \
        gp.sim_out_prefix + \
        'individual_trees_' +  tag + '_rep' + str(rep) + gp.sim_out_suffix
    tree_fns.append(individual_trees_fn)
    f = open(individual_trees_fn, 'w')
    for n, tree in trees[rep]:
        for i in range(n):
            f.write(tree + '\n')
    f.close()

#####
# run seq-gen
#####

seq_fns = []
for rep in range(len(tree_fns)):
    fn = tree_fns[rep]
    seq_gen_out_fn = gp_dir + gp.sim_out_dir + '/seq-gen/temp/' + gp.sim_out_prefix + \
        'seq_gen_' +  tag + '_rep' + str(rep) + gp.sim_out_suffix
    seq_fns.append(seq_gen_out_fn)
    # JC69 model equivalent to HKY model but with base frequencies set
    # equal and transition/transversion rates set equal (which are the
    # default options for -c and -t); length is same as given to ms;
    # the -l1 option is so that we get one site per tree (since we've
    # made a file with one tree per site)
    seq_gen_command = \
        '~/software/Seq-Gen.v1.3.3/source/seq-gen -mHKY -l1' + \
        ' < ' + fn + ' > ' + seq_gen_out_fn
    print seq_gen_command
    os.system(seq_gen_command)

#####
# convert seq-gen output to fasta format (one file for each rep)
#####


# one file for each rep

sample_labels = [str(x) for x in range(1, num_samples+1)]
for rep in range(len(seq_fns)):
    fn = seq_fns[rep]
    f = open(fn, 'r')
    fasta_fn = seq_gen_out_fn = gp_dir + gp.sim_out_dir + '/seq-gen/' + \
        gp.sim_out_prefix + \
        'seq_gen_' +  tag + '_rep' + str(rep) + '.fasta'
    fasta_f = open(fasta_fn, 'w')
    seqs = dict(zip(sample_labels, [''] * num_samples))
    for i in range(num_sites):
        a, sites_per_tree = f.readline().strip().split(' ')
        assert int(a) == num_samples, a + ' ' + str(num_samples)
        for s in sample_labels:
            line = f.readline().strip()
            label = line[:line.find(' ')]
            seq = line[-int(sites_per_tree):]
            seqs[s] += seq
    for s in sample_labels:
        fasta_f.write('> ' + label_to_name[s] + '\n')
        fasta_f.write(seqs[s] + '\n')
fasta_f.close()
