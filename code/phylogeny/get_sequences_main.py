# creates fasta file of genomes of all strains (chromosomes appended),
# with position for every site in cerevisiae reference
# creates same file but with only introgressed sites
# creates same file but with only nonintrogressed sites
# creates same above three files but with alignment columns with any gaps removed

import sys
import os
from get_gene_set import *
import math
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../align/')
import align_helpers
sys.path.insert(0, '../analyze/')
import gene_predictions
sys.path.insert(0, '../misc/')
import read_table
import read_fasta
import write_fasta
import mystats

def get_inds_from_alignment(fn, rind, sind):
    headers, seqs = read_fasta.read_fasta(fn)
    n = len(seqs[0])
    ri = -1
    si = -1
    ps = []
    for i in range(n):
        s_gap = True
        if seqs[sind][i] != gp.gap_symbol:
            si += 1
            s_gap = False
        if seqs[rind][i] != gp.gap_symbol:
            ri += 1
            if s_gap:
                ps.append(None)
            else:
                ps.append(str(si))
    return ps

tag = 'u3_i.001_tv_l1000_f.01'
suffix = '_filtered'
gp_dir = '../'

gap_sites = set([])
int_sites = set([])
chrm_offsets = {}

fn_all = gp.analysis_out_dir_absolute + tag + \
         '/seqs_for_phylogeny_all' + gp.fasta_suffix
fn_all_ungapped = gp.analysis_out_dir_absolute + tag + \
                  '/seqs_for_phylogeny_all_ungapped' + gp.fasta_suffix
fn_int = gp.analysis_out_dir_absolute + tag + \
         '/seqs_for_phylogeny_int' + gp.fasta_suffix
fn_int_ungapped = gp.analysis_out_dir_absolute + tag + \
                  '/seqs_for_phylogeny_int_ungapped' + gp.fasta_suffix
fn_nonint = gp.analysis_out_dir_absolute + tag + \
            '/seqs_for_phylogeny_nonint' + gp.fasta_suffix
fn_nonint_ungapped = gp.analysis_out_dir_absolute + tag + \
                     '/seqs_for_phylogeny_nonint_ungapped' + gp.fasta_suffix

#======
# write all sites, including gaps
#======

print 'writing file with all sites'

f_all = open(fn_all, 'w')

# master reference (cerevisiae)
print '*', gp.master_ref
f_all.write('>' + gp.master_ref + '\n')
chrm_offset = 0
for chrm in gp.chrms:
    seq = read_fasta.read_fasta(gp.ref_dir[gp.master_ref] + \
                                gp.ref_fn_prefix[gp.master_ref] + \
                                '_chr' + chrm + gp.fasta_suffix)[1][0]
    f_all.write(seq)
    chrm_offsets[chrm] = chrm_offset
    chrm_offset += len(seq)
f_all.write('\n')

# other reference (paradoxus)
other_ref_strain = gp.ref_fn_prefix[gp.alignment_ref_order[1]]
print '*', other_ref_strain
f_all.write('>' + other_ref_strain + '\n')
for chrm in gp.chrms:
    align_fn = gp_dir + gp.alignments_dir + \
               '_'.join(gp.alignment_ref_order) + '_chr' + chrm + \
               '_mafft' + gp.alignment_suffix
    ps = get_inds_from_alignment(align_fn, 0, 1)
    chrom_seq = read_fasta.read_fasta(gp.ref_dir[other_ref_strain] + \
                                      gp.ref_fn_prefix[other_ref_strain] + \
                                      '_chr' + chrm + gp.fasta_suffix)[1][0]
    for i in range(len(ps)):
        try:
            x = int(ps[i])
            f_all.write(chrom_seq[x])
        except:
            f_all.write(gp.gap_symbol)
            gap_sites.add(chrm_offsets[chrm] + i)

f_all.write('\n')

# all other strains
s = align_helpers.get_strains(align_helpers.flatten(gp.non_ref_dirs.values()))

for strain, d in s:
    print '*', strain
    f_all.write('>' + strain + '\n')
    for chrm in gp.chrms:
        align_fn = gp_dir + gp.alignments_dir + \
                   '_'.join(gp.alignment_ref_order) + '_' + \
                   strain + '_chr' + chrm + \
                   '_mafft' + gp.alignment_suffix
        ps = get_inds_from_alignment(align_fn, 0, 2)
        chrom_seq = read_fasta.read_fasta(d + '/' + strain + \
                                          '_chr' + chrm + gp.fasta_suffix)[1][0]
        for i in range(len(ps)):
            try:
                x = int(ps[i])
                f_all.write(chrom_seq[x])
            except:
                f_all.write(gp.gap_symbol)
                gap_sites.add(chrm_offsets[chrm] + i)
                
    f_all.write('\n')

f_all.close()

#======
# write all the other files with/without gaps and introgressed sites
#======

print 'writing all other files'

# get introgressed sites
fn_regions = gp.analysis_out_dir_absolute + tag + '/' + \
             'introgressed_blocks' + suffix + '_par_' + tag + '_summary_plus.txt'
regions, l = read_table.read_table_rows(fn_regions, '\t')
for rid in regions:
    region = regions[rid]
    for i in range(int(region['start']), int(region['end'])+1):
        int_sites.add(chrm_offsets[region['chromosome']] + i)

# write files
f = open(fn_all, 'r')
f_all_ungapped = open(fn_all_ungapped, 'w')
f_int = open(fn_int, 'w')
f_int_ungapped = open(fn_int_ungapped, 'w')
f_nonint = open(fn_nonint, 'w')
f_nonint_ungapped = open(fn_nonint_ungapped, 'w')
for i in range(len(s) + 2):
    header = f.readline()[:-1]
    print '*', header[1:]
    seq = f.readline()[:-1]

    f_all_ungapped.write(header + '\n')
    f_int.write(header + '\n')
    f_int_ungapped.write(header + '\n')
    f_nonint.write(header + '\n')
    f_nonint_ungapped.write(header + '\n')

    for i in range(len(seq)):
        gap = i in gap_sites
        intd = i in int_sites
        a = seq[i]
        if not gap:
            f_all_ungapped.write(a)
        if intd:
            f_int.write(a)
            if not gap:
                f_int_ungapped.write(a)
        else:
            f_nonint.write(a)
            if not gap:
                f_nonint_ungapped.write(a)
        
    f_all_ungapped.write('\n')
    f_int.write('\n')
    f_int_ungapped.write('\n')
    f_nonint.write('\n')
    f_nonint_ungapped.write('\n')

f.close()
f_all_ungapped.close()
f_int.close()
f_int_ungapped.close()
f_nonint.close()
f_nonint_ungapped.close()




