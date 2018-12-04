import re
import sys
import os
import math
import Bio.SeqIO
import copy
import gzip
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

tag = 'u3_i.001_tv_l1000_f.01'
suffix = '_filtered'
region_id = 'r2793' #'r2259' #'r2402'
make_dbs = False

gp_dir = '../'
s = align_helpers.get_strains(gp.non_ref_dirs[gp.master_ref])
s = dict(s)

# read in filtered regions
fn_regions = gp.analysis_out_dir_absolute + tag + '/' + \
             'introgressed_blocks' + suffix + '_par_' + tag + '_summary_plus.txt'
regions, l = read_table.read_table_rows(fn_regions, '\t')

region_info = regions[region_id]

coord_fn = gp.analysis_out_dir_absolute + 'coordinates/' + \
           gp.master_ref + '_to_' + region_info['strain'] + \
           '_chr' + region_info['chromosome'] + '.txt.gz'
f_coord = gzip.open(coord_fn, 'rb')
ref_ind_to_strain_ind = [float(line[:-1]) for line in f_coord.readlines()]
strain_start = int(math.ceil(ref_ind_to_strain_ind[int(region_info['start'])]))
strain_end = int(math.floor(ref_ind_to_strain_ind[int(region_info['end'])]))

strain_fn = s[region_info['strain']] + region_info['strain'] + '_chr' + \
            region_info['chromosome'] + gp.fasta_suffix
strain_seq = read_fasta.read_fasta(strain_fn)[1][0][strain_start:strain_end+1]
query_fn = 'query.txt'
f = open(query_fn, 'w')
f.write(strain_seq)
f.close()

par_strains = os.listdir(gp.non_ref_dirs['CBS432'][0])

if make_dbs:
    for par_strain in par_strains:
        db_fn = gp.non_ref_dirs['CBS432'][0] + par_strain + '/assembly/genome.fa'
        cmd_string = gp.blast_install_path + 'makeblastdb' + \
                     ' -in ' + db_fn + \
                     ' -dbtype nucl'
        print cmd_string
        os.system(cmd_string)
            

outfmt = '"6 sseqid slen evalue bitscore sstart send length pident"'
for par_strain in par_strains:
    db_fn = gp.non_ref_dirs['CBS432'][0] + par_strain + '/assembly/genome.fa'
    out_fn = 'blast_temp.txt'

    cmd_string = gp.blast_install_path + 'blastn' + \
                 ' -db ' + db_fn + \
                 ' -query ' + query_fn + \
                 ' -out ' + out_fn + \
                 ' -outfmt ' + outfmt
    #print cmd_string
    os.system(cmd_string)

    lines = open(out_fn, 'r').readlines()
    print par_strain, lines[0]

"""
headers, strain_chrm_seqs = read_fasta.read_fasta(strain_fn)
strain_seq = ''
for i in range(len(headers)):

    current_chrm = int(headers[i][headers[i].find('chr')+3:])
    if current_chrm == chrms_ara[region_info['chromosome']]:
        strain_seq = strain_chrm_seqs[strain_start:strain_end+1]
        break

s[region_info['strain']] + region_info['strain'] + '/assembly/' + \
            'genome.fa'

"""
