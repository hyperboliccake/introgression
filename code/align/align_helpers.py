import sys
import os
sys.path.insert(0, '..')
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
        fns = filter(lambda x: '_chr' in x, fns)    
        num_files = len(fns)
        if num_files == 0:
            print 'found no chromosome sequence files in', d, '(perhaps you should check the _chr naming convention?)'
        fns = list(set([x[:x.find('_chr')] for x in fns]))
        num_strains = len(fns)
        assert num_files == num_strains * len(gp.chrms), \
            'some strains in ' + d + ' have the wrong number of chromosome sequence files'
        entries = [(x, d) for x in fns]
        s += entries
    return s

def concatenate_fasta(fns, f):
    
    f = open(f, 'w')
    for fn in fns:
        f_current = open(fn, 'r')
        header = f_current.readline()[:-1]
        f.write(header + ' ' + fn + '\n')
        f.write(f_current.read())
        f_current.close()
        f.write('\n')
    f.close()
