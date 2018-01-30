import gzip

def write_fasta(headers, seqs, fn, gz=False):
    
    f = None
    if gz:
        f = gzip.open(fn + '.gz', 'wb')
    else:
        f = open(fn, 'w')
    for i in range(len(headers)):
        if headers[i][0] != '>':
            headers[i] = '>' + headers[i]
        if headers[i][-1] != '\n':
            headers[i] += '\n'
        f.write(headers[i])
        f.write(seqs[i] + '\n')
    f.close()
