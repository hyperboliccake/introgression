# optional, probably unnecessary step---linked loci can be okay if
# they're snps from across the genome

import sys
import os
import gzip
import predict
from collections import defaultdict
import gene_predictions
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc/')
import read_fasta
import read_table
import seq_functions

args = predict.process_predict_args(sys.argv[2:])

chrm = gp.chrms[int(sys.argv[1])]

# maybe getting strains should be simpler
strains = [line.split('\t')[0] for line in \
           open(gp.analysis_out_dir_absolute + args['tag'] + \
                '/state_counts_by_strain.txt', 'r').readlines()[1:]]

nucs = set(['a', 't', 'g', 'c'])

out_dir = gp.analysis_out_dir_absolute + args['tag'] + '/structure/'
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

gp_dir = '../'


##======
# use program ldselect to find set of tag snps all in low LD for
# specified chromosome
##======

# input file for ldselect is formatted so that each row is a snp and
# each column is the genotype for a strain, e.g.

fn = out_dir + 'ldselect_input_chr' + chrm + '.tsv'
f = open(fn, 'w')
snps = defaultdict(list)
# loop through all the strains
for strain in strains:
    print '-', strain
    # read multiple alignment file for this strain with the master
    # reference (and other references which we don't care about
    # here)
    headers, seqs = read_fasta.read_fasta(gp_dir + gp.alignments_dir + \
                                          '_'.join(gp.alignment_ref_order) + \
                                          '_' + strain + '_chr' + chrm + \
                                          '_mafft.maf')
    # look at all alignment columns, keeping track of the index in
    # the master reference
    i = 0
    for c in range(len(seqs[0])):
        # if the master reference doesn't have a gap in this
        # column, then store the allele that the current strain
        # has at this site
        if seqs[0][c] != gp.gap_symbol and seqs[0][c] != gp.unsequenced_symbol:
            snps[i].append(seqs[-1][c])
            i += 1

# get reference sequence (unaligned, without gaps)
# TODO correct alignment file location
ref_seq = read_fasta.read_fasta(gp_dir + gp.alignments_dir + \
                                '_'.join(gp.alignment_ref_order) + \
                                '_' + strains[0] + '_chr' + chrm + \
                                '_mafft.maf')[1][0].replace(gp.gap_symbol, '')
open(out_dir + 'chromosome_lengths.txt', 'a').write(chrm + '\t' + \
                                                    str(len(ref_seq)) + '\n')

# loop through all the sites we collected above
for snp in snps.keys():
    # and only keep the ones that are polymorphic but have no gaps
    # or unsequenced symbols (including in master reference)
    snp_set = set(snps[snp])
    if len(snp_set - nucs) == 0 and ref_seq[snp] in nucs and len(snp_set) > 1:
        # naming convention of snps: chromosome_position
        # TODO do names have to be integers and/or equal in length?
        snp_id = str(snp)
        # write row for master reference
        f.write(snp_id + '\t' + \
                gp.alignment_ref_order[0] + '\t' + \
                ref_seq[snp] + '\n')
        # and one row for each of the other strains
        for si in range(len(strains)):
            f.write(snp_id + '\t' + \
                    strains[si] + '\t' + \
                    snps[snp][si] + '\n')

f.close()

"""
# run ldselect on this input file
fn_out = fn.replace('input', 'output')
os.system('perl ' + gp.ldselect_install_path + 'ldSelect.pl -pb ' + fn + ' > ' + fn_out)

# extract one tag snp from each set of equivalent tag snps from
# ldselect output file
# TODO correct output file name/location
f = open(fn_out, 'r')
line = f.readline()
all_tag_snps = set([])
while line != '':
    if "TagSnps:" in line:
        tag_snps = line[line.find(':')+1:].split()
        all_tag_snps.add(int(tag_snps[0]))
    line = f.readline()
f.close()

f = open('structure/tag_snps_chr' + chrm + '.txt', 'w')
f.write('\n'.join([str(x) for x in sorted(list(all_tag_snps))]))
f.close()
"""
