import re

chrms_roman = ['0', 'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XIV']


f = open('../../results/introgressed_id.txt', 'r')
strain = 'yjm456'
line = f.readline()
region_count = 1
f = open('../../results/introgressed_id.txt', 'r')
f_out = open('../../results/strains/' + strain + '_introgressed.bed', 'w')
while line != '':
    if strain in line:
        line = line.strip().split(',')
        
        m = re.search('chr(?P<chrm>[A-Z]+)\.', line[0])
        chrm = m.group('chrm')
        chrm = 'chr' + str(chrms_roman.index(chrm))
        
        d = line[2].find('-')
        start = line[2][:d]
        end = line[2][d+1:]
        name = 'region' + str(region_count)
        region_count += 1
        score = '0'
        strand = line[1][line[1].find(' strand') - 1]
        f_out.write(chrm + '\t' + start + '\t' + end + '\t' + name + '\t' + score + '\t' + strand + '\n')
    line = f.readline()
f.close()
f_out.close()
        
