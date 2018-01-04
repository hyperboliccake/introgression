import os
import sys

fns = filter(lambda x: 'mafft' in x, os.listdir('../../alignments/genbank/'))
for fn in fns:
    print fn
    f = open('../../alignments/genbank/' + fn, 'r')
    g = open('../../alignments/genbank/' + fn + '.fixed', 'w')
    line = f.readline()
    while line != '':
        if line[0] == '>':
            if 'CBS432' in line:
                line_new = '>CBS432 ' + line.strip().split(' ')[-1] + '\n'
                line = line_new
        g.write(line)
        line = f.readline()
    g.close()
    f.close()
    os.system('mv ' + '../../alignments/genbank/' + fn + '.fixed ' + '../../alignments/genbank/' + fn) 
