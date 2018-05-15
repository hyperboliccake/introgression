import gzip

def write_fasta(headers, seqs, fn, gz=False):
    
    f = None
    if gz:
        f = gzip.open(fn + '.gz', 'wb')
    else:
        f = open(fn, 'w')
    for i in range(len(headers)):
        header = headers[i]
        if header[0] != '>':
            header = '>' + header
        if header[-1] != '\n':
            header += '\n'
        f.write(header)
        f.write(seqs[i] + '\n')
    f.close()
