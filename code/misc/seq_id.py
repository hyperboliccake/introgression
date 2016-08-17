# calculate sequence identity between all pairs of cerevisiae in 100
# genomes set, and also with paradoxus reference

# do this from the alignments? just the three-way portions?

import os
import re

chrms = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XIV']

alignment_dir = '../../alignments/genbank/'
prefix = 'S288c_CBS432_'
fns = os.listdir(alignment_dir)
strains = set()
for fn in fns:
    m = re.match(prefix + '(?P<strain>)_chr')
    strains.add(m.group('strain'))

for strain in strains:
    for chrm in chrms:
        f = open(alignment_dir + prefix + strain + '_chr' + chrm + '.maf')
    
