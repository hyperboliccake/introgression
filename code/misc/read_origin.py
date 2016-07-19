import sys
import re

f = open(sys.argv[1], 'r')
strain = sys.argv[2]
chrm = sys.argv[3]
line = f.readline()
while line != '':

    m = re.search('Saccharomyces cerevisiae (?P<strain>[a-zA-Z0-9]+) chromosome (?P<chrm>[IVXM]+)', line)
    if m != None:
        current_strain = m.group('strain')
        current_chrm = m.group('chrm')
        if current_strain.lower() == strain.lower() and current_chrm == chrm:
            line = f.readline()
            while line[:6] != 'ORIGIN':
                line = f.readline()
            line = f.readline()
            seq = ''
            while line != '//\n':
                line = line.split()
                for x in line[1:]:
                    seq += x
                line = f.readline()
            fo = open(sys.argv[4], 'w')
            fo.write(seq.upper() + '\n')
            fo.close()
            break
    line = f.readline()

