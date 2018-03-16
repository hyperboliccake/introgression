import sys
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc')
import read_fasta
import write_fasta

def read_intervals(fn):

    f = open(fn, 'r')
    f.readline() # header
    line = f.readline()
    intervals = []
    while line != '':
        start, end = [int(x) for x in line[:-1].split(' - ')]
        intervals.append((start, end))
        line = f.readline()
    f.close()
    return intervals

def mask(fn, masked_fn, intervals_fn):

    headers, seqs = read_fasta.read_fasta(fn)
    seq = list(seqs[0])
    intervals = read_intervals(intervals_fn)
    for start, end in intervals:
        for i in range(start, end + 1):
            seq[i] = gp.unsequenced_symbol
    seq = ''.join(seq)
    write_fasta.write_fasta(headers, [seq], masked_fn)

