# columns:
# position in strain
# position in cer ref
# gene
# in ORF?

import re
import sys
import os
import copy
import gzip
from annotate_positions import *
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../align/')
import align_helpers

##======
# get strains
##======

i = int(sys.argv[1])
s = align_helpers.get_strains(align_helpers.flatten(gp.non_ref_dirs.values()))
strain, d = s[i]

##======
# get genes on each chromosome
##======

genes_by_chrm = {}
for chrm in gp.chrms:
    fn = gp.analysis_out_dir_absolute + gp.master_ref + '_chr' + chrm + \
         '_genes.txt'
    genes_by_chrm[chrm] = get_genes(fn)

##======
# loop through all strains and chromosomes, generating annotated
# position file for each
##======

coord_dir = gp.analysis_out_dir_absolute + 'coordinates/'
if not os.path.exists(coord_dir + 'annotated'):
    os.makedirs(coord_dir + 'annotated')

for chrm in gp.chrms:

    print strain, chrm

    fn = strain + '_to_' + gp.master_ref + '_chr' + chrm + '.txt.gz'

    fn_orfs = d + 'orfs/' + strain + '_chr' + chrm + \
              '_orfs' + gp.fasta_suffix
    orfs = get_orfs(fn_orfs)

    fn_out = coord_dir + 'annotated/' + fn
    coords = [float(line) for line in gzip.open(coord_dir + fn, 'rb').readlines()]
    write_annotated_file(coords, genes_by_chrm[chrm], orfs, fn_out)
    




#for strain, d in s:

    #m = re.search('(?P<strain1>[a-zA-Z0-9]+)_to_(?P<strain2>[a-zA-Z0-9]+)_chr(?P<chrm>[IVXM]+)', fn)
    #if m == None:
    #    continue
    #strain1 = m.group('strain1')
    #strain2 = m.group('strain2')
    #chrm = m.group('chrm')
    
    #if strain1 == gp.master_ref:
     #   continue

    # don't deal with paradoxus just for now
    #if strain1 in gp.alignment_ref_order or strain2 != gp.master_ref:
    #    continue

    #print fn

