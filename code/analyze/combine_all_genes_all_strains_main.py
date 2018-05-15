import re
import sys
import os
import math
import Bio.SeqIO
import copy
from combine_all_strains import *
import gene_predictions
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../align/')
import align_helpers
sys.path.insert(0, '../misc/')
import read_table
import read_fasta
import write_fasta
import mystats
import subprocess

tag = sys.argv[1]
chrm = sys.argv[2]

strains_for_each_gene_fn = gp.analysis_out_dir_absolute + tag + '/' + \
                           'strains_for_each_gene_chr' + chrm + '_' + tag + '.txt'
genes = [line.split()[0] for line in open(strains_for_each_gene_fn, 'r').readlines()]
i = 1
for gene in genes:
    gene = genes[i-1]
    print gene + ' (' + str(i) + '/' + str(len(genes)) + ')'
    sys.stdout.flush()
    i += 1
    if os.path.isfile(gp.analysis_out_dir_absolute + tag + '/genes/' + \
                      gene + '/' + gene + '_introgressed' + gp.alignment_suffix):
        print 'already done'
        sys.stdout.flush()
    else:
        cmd = 'python combine_gene_all_strains_main.py ' + tag + ' ' + \
              gene + ' ' + chrm
        r = subprocess.call(cmd, shell=True)
        sys.stdout.flush()
        print r

