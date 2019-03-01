import sys
import os
import global_params as gp

def flatten(l):
    return [item for sublist in l for item in sublist]

def get_strains(dirs):
    # get all non-reference strains of cerevisiae and paradoxus; could
    # generalize this someday...

    # s is list of tuples: (strain_name, directory it came from)
    s = []
    for d in dirs:
        fns = os.listdir(d)
        # only look at fasta files in the directory
        fns = filter(lambda x: x.endswith(gp.fasta_suffix), fns)
        # only look at files containing '_chr' which should be chromosome
        # sequence files
        fns = list(filter(lambda x: '_chr' in x, fns))
        num_files = len(fns)
        if num_files == 0:
            print('found no chromosome sequence files in', d, '(perhaps you should check the _chr naming convention?)')
        fns = list(set([x[:x.find('_chr')] for x in fns]))
        num_strains = len(fns)
        # might be greater because of masked sequence files
        assert num_files >= num_strains * len(gp.chrms), \
            'some strains in ' + d + ' are missing chromosome sequence files'
        entries = [(x, d) for x in fns]
        s += entries
    return sorted(s)

def concatenate_fasta(fns, names, f):
    
    f = open(f, 'w')
    for i in range(len(fns)):
        fn = fns[i]
        f_current = open(fn, 'r')
        header = f_current.readline()[:-1]
        f.write('> ' + names[i] + ' ' + fn + '\n')
        f.write(f_current.read())
        f_current.close()
        f.write('\n')
    f.close()
