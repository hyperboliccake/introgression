import sys
import os
import gzip
import math
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../align/')
import align_helpers
sys.path.insert(0, '../misc/')
import read_fasta
import read_table

shared_regions, l = \
    read_table.read_table_rows('shared_introgression_nonsingleton_polymorphism.txt', \
                               '\t')
f = open('shared_introgression_nonsingleton_polymorphism_3strains.txt', 'w')
f.write('region_number\tchromosome\tstart\tend\tin_3strains\tin_only_3strains\tpi\tfrac_poly\tnum_poly\tnum_total\n')
s3 = set(['yjm1078', 'yjm1242', 'yjm248'])
for region_number in shared_regions.keys():
    strains = set(shared_regions[region_number]['strain_list'].split(','))
    f.write(region_number + '\t')
    f.write(shared_regions[region_number]['chromosome'] + '\t')
    f.write(shared_regions[region_number]['start'] + '\t')
    f.write(shared_regions[region_number]['end'] + '\t')
   
    if strains.intersection(s3) == set(): # none of three strains present
        f.write('FALSE\tFALSE\t')
    elif strains - s3 == set(): # no other strain present
        f.write('TRUE\tTRUE\t')
    else:
        f.write('TRUE\tFALSE\t')

    f.write(shared_regions[region_number]['pi'] + '\t')
    f.write(shared_regions[region_number]['frac_poly'] + '\t')
    f.write(shared_regions[region_number]['num_poly'] + '\t')
    f.write(shared_regions[region_number]['num_total'] + '\n')
    f.flush()
