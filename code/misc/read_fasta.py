
def read_fasta(fn):

    headers = []
    seqs = []

    f = open(fn, 'r')
    line = f.readline()
    while line[0] != '>':
        line = f.readline()
    while True:
        h = line[:-1]
        s = ''
        line = f.readline()
        while line != '' and line[0] != '>':
            s += line[:-1]
            line = f.readline()
        headers.append(h)
        seqs.append(s)
        if line == '':
            break

    return headers, seqs
