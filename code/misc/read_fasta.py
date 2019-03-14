import gzip
import numpy as np


def read_fasta(fn, gz=False):

    headers = []
    seqs = []

    f = None
    if gz:
        f = gzip.open(fn, 'rb')
    else:
        f = open(fn, 'r')
    line = f.readline()

    while line[0] != '>':
        line = f.readline()

    while True:
        h = line[:-1]
        s = []
        line = f.readline()
        while line != '' and line[0] != '>':
            s += list(line[:-1])
            line = f.readline()
        headers.append(h)
        seqs.append(s)
        if line == '':
            break
    f.close()

    return headers, np.asarray(seqs)
