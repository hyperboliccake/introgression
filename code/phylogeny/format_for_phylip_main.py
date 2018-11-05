import sys
import os
sys.path.insert(0, '../misc/')
import read_fasta

def pad(s, n):
    s = s.strip()
    return s[:n] + (n - len(s)) * ' '

headers, seqs = read_fasta.read_fasta(sys.argv[1])

fp = open(sys.argv[2], 'w')
fp.write(str(len(headers)) + ' ' + str(len(seqs[0])) + '\n')
for i in range(len(headers)):
    h = pad(headers[i][1:], 10)
    fp.write(h + seqs[i] + '\n')

fp.close()
