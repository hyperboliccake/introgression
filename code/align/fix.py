import os
import sys

fns = filter(lambda x: 'mafft' in x and '1385' not in x, os.listdir('../../alignments/genbank/'))
for fn in fns:
    print fn
    f = open('../../alignments/genbank/' + fn, 'r')
    g = open('../../alignments/genbank/' + fn + '.fixed', 'w')
    line = f.readline()
    while line != '':
        if line[0] == '>':
            if len(line.split(' ')) == 2:
                s = line.split(' ')[0][1:]
                c = line.split(' ')[1][:-1]
                x = ''
                if s == 'S288c':
                    x = '/net/akey/vol2/aclark4/nobackup/100_genomes/genomes/S288c_SGD-R64/S288c_SGD-R64_' + c + '.fa'
                elif s == 'CBS432':
                    x = '/net/gs/vol1/home/aclark4/projects/introgression/data/CBS432/CBS432_' + c + '.fa'
                else:
                    assert s[:3] == 'yjm', s
                    x = '/net/akey/vol2/aclark4/nobackup/100_genomes/genomes_gb/' + s + '_' + c + '.fa'
                g.write(line[:-1] + ' ' + x + '\n')
            else:
                assert len(line.split(' ')) == 3, line.split(' ')
                g.write(line)                
        else:
            g.write(line)
        line = f.readline()
    g.close()
    f.close()
    os.system('mv ' + '../../alignments/genbank/' + fn + '.fixed ' + '../../alignments/genbank/' + fn) 
