# Loop through all introgressed genes (might be just a small part)
# that have paralogs
# Extract introgressed portion of gene 
# Blast that portion against:
# - Cerevisiae gene
# - Paradoxus gene (region aligned to cerevisiae gene)
# - Cerevisiae paralog
# - Paradoxus paralog (region aligned to cervisiae paralog)
# If best match is
# - Cerevisiae gene -> whatever, probably just not a great call
# - Paradoxus gene -> looks introgressed, as expected
# - Cerevisiae paralog -> paralogous gene conversion
# - Paradoxus paralog -> interesting...


import re
import sys
import os
import math
import Bio.SeqIO
import copy
import gzip
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

postprocess = False

tag = 'u3_i.001_tv_l1000_f.01'
query_fn = 'blastquerytemp.txt'
db_fn = 'blastdbtemp' + gp.fasta_suffix
out_fn = 'blastouttemp.txt'

outfmt = '"6 sseqid slen evalue bitscore sstart send pident"'

gp_dir = '../'

s = align_helpers.get_strains(align_helpers.flatten(gp.non_ref_dirs.values()))
strain_dirs = dict(s)

# dict of dicts keyed by region id and column names; includes filtered
# and unfiltered regions
region_to_genes = {}
f = open(gp.analysis_out_dir_absolute + tag + \
         '/genes_for_each_region_' + tag + '.txt', 'r')
line = f.readline()
while line != '':
    line = line.split('\t')
    region_id = line[0]
    genes = line[2::2]
    region_to_genes[region_id] = genes
    line = f.readline()
f.close()

# dict of lists keyed by region id
t_regions_filtered, l = \
    read_table.read_table_rows(gp.analysis_out_dir_absolute + tag + \
                               '/introgressed_blocks_filtered_par_' + tag + \
                               '_summary_plus.txt', \
                               '\t', header=True)


# get all genes introgressed in at least one strain
gene_to_regions = {}
gene_to_regions_filtered = {}
for region_id in region_to_genes:
    genes = region_to_genes[region_id]
    for gene in genes:
        if not gene_to_regions.has_key(gene):
            gene_to_regions[gene] = []
        gene_to_regions[gene].append(region_id)
    if region_id in t_regions_filtered:
        for gene in genes:
            if not gene_to_regions_filtered.has_key(gene):
                gene_to_regions_filtered[gene] = []
            gene_to_regions_filtered[gene].append(region_id)

# read in paralogs
fn_paralogs = '../../data/S288c_paralogs.tsv'
f_paralogs = open(fn_paralogs, 'r')
lines = [line.strip().split('\t') for line in f_paralogs.readlines()]
f_paralogs.close()
paralogs = {}
for line in lines:
    if line[0] != '""' and line[3] != '""':
        paralogs[line[0]] = line[3]

# read in all gene coordinates
gene_coords = {}
for chrm in gp.chrms:
    f = open(gp.analysis_out_dir_absolute + \
             'S288c_chr' + chrm + '_genes.txt', 'r')
    lines = [line.strip().split('\t') for line in f.readlines()]
    f.close()
    for line in lines:
        gene_coords[line[0]] = (chrm, int(line[1]), int(line[2]))

# loop through all introgressed genes
keys = ['cer_gene', 'par_gene', 'cer_paralog', 'par_paralog', 'none']
all_rankings = dict(zip(keys, [[] for key in keys]))

genes_to_analyze = gene_to_regions_filtered.keys()
if postprocess:
    genes_to_analyze = [line.split('\t')[0] for line in \
                        open('check_paralogs_out_cer_paralog.tsv', 'r').readlines()]
    genes_to_analyze = list(set(genes_to_analyze))

ip = 0
for gene in genes_to_analyze:
    if gene not in paralogs:
        continue

    print ip
    ip += 1

    chrm, ref_gene_start, ref_gene_end = gene_coords[gene]

    gene_headers, gene_seqs = \
        read_fasta.read_fasta(gp.analysis_out_dir_absolute + tag + '/genes/' + \
                              gene + '/' + gene + '_from_alignment.fa')
    gene_headers = [x[1:].strip() for x in gene_headers]
    strain_seqs = dict(zip(gene_headers, gene_seqs))

    cer_seq = strain_seqs['S288c']
    par_seq = strain_seqs['CBS432']

    paralog = paralogs[gene]
    gene_headers, gene_seqs = \
        read_fasta.read_fasta(gp.analysis_out_dir_absolute + tag + '/genes/' + \
                              paralog + '/' + paralog + '_from_alignment.fa')
    gene_headers = [x[1:].strip() for x in gene_headers]
    strain_paralog_seqs = dict(zip(gene_headers, gene_seqs))

    cer_paralog_seq = strain_paralog_seqs['S288c']
    par_paralog_seq = strain_paralog_seqs['CBS432']

    f = open(db_fn, 'w')
    f.write('> cer_gene\n')
    f.write(cer_seq + '\n')
    f.write('> par_gene\n')
    f.write(par_seq + '\n')
    f.write('> cer_paralog\n')
    f.write(cer_paralog_seq + '\n')
    f.write('> par_paralog\n')
    f.write(par_paralog_seq + '\n')
    f.close()

    cmd_string = gp.blast_install_path + 'makeblastdb' + \
                 ' -in ' + db_fn + \
                 ' -dbtype nucl'
    os.system(cmd_string)

    strain_intd_seqs = {}
    for region in gene_to_regions_filtered[gene]:
        strain = t_regions_filtered[region]['strain']

        ref_region_start = int(t_regions_filtered[region]['start'])
        ref_region_end = int(t_regions_filtered[region]['end'])

        ref_to_strain_coords = [float(x[:-1]) for x in \
                                gzip.open(gp.analysis_out_dir_absolute + \
                                          'coordinates/S288c_to_' + strain + \
                                          '_chr' + chrm + '.txt.gz').readlines()]

        gene_start = int(max(0, math.ceil(ref_to_strain_coords[ref_gene_start])))
        gene_end = int(math.floor(ref_to_strain_coords[ref_gene_end]))
        
        region_start = int(max(0, math.ceil(ref_to_strain_coords[ref_region_start])))
        region_end = int(math.floor(ref_to_strain_coords[ref_region_end]))

        start = max(gene_start, region_start)
        end = min(gene_end, region_end)

        chrom_seq = read_fasta.read_fasta(strain_dirs[strain] + strain + '_chr' + \
                                          chrm + gp.fasta_suffix)[1][0]
        seq = chrom_seq[start:end+1]

        if not strain_intd_seqs.has_key(strain):
            strain_intd_seqs[strain] = chrom_seq[gene_start:gene_end+1].lower()
        relative_start = start - gene_start
        relative_end = end - gene_start
        strain_intd_seqs[strain] = \
            strain_intd_seqs[strain][:relative_start] + \
            strain_intd_seqs[strain][relative_start:relative_end+1].upper() + \
            strain_intd_seqs[strain][relative_end+1:]

        f = open(query_fn, 'w')
        f.write(seq + '\n')
        f.close()

        cmd_string = gp.blast_install_path + 'blastn' + \
                     ' -db ' + db_fn + \
                     ' -query ' + query_fn + \
                     ' -out ' + out_fn + \
                     ' -outfmt ' + outfmt
        print cmd_string
        os.system(cmd_string)

        if os.stat(out_fn).st_size == 0:
            cmd_string = gp.blast_install_path + 'blastn' + \
                         ' -db ' + db_fn + \
                         ' -query ' + query_fn + \
                         ' -out ' + out_fn + \
                         ' -task "blastn-short"' + \
                         ' -outfmt ' + outfmt
            print cmd_string
            os.system(cmd_string)
            
        lines = open(out_fn, 'r').readlines()
        best_key = 'none'
        if len(lines) != 0:
            best_key = lines[0].split('\t')[0]
        pidents = dict(zip(keys[:-1], [[] for k in keys[:-1]]))
        for line in lines:
            line = line.strip().split('\t')
            key = line[0]
            pident = line[-1]
            pidents[key].append(pident)

        r = [gene, region, strain]
        for k in keys[:-1]:
            if pidents[k] == []:
                r.append('NA')
            else:
                r.append(','.join(pidents[k]))
        r = tuple(r)
        all_rankings[best_key].append(r)

    # write reference genes and paralogs and all introgressed
    # genes to file and then align
    fn = gp.analysis_out_dir_absolute + tag + '/paralogs/' + \
         gene + gp.fasta_suffix
    headers = ['S288c ' + gene, 'CBS432 ' + gene, \
               'S288c ' + paralog, 'CBS432 ' + paralog]
    seqs = [cer_seq.lower(), par_seq.lower(), \
            cer_paralog_seq.lower(), par_paralog_seq.lower()]
    for strain in strain_intd_seqs:
        headers.append(strain + ' ' + gene)
        seqs.append(strain_intd_seqs[strain])
    write_fasta.write_fasta(headers, seqs, fn)

    aligned_fn = fn.replace(gp.fasta_suffix, gp.alignment_suffix)
    cmd_string = gp.mafft_install_path + '/mafft ' + \
                 ' --quiet --reorder --preservecase ' + \
                 fn + ' > ' + aligned_fn
    os.system(cmd_string)
        
f = open('check_paralogs_out.tsv', 'w')
f.write('category\tnum_total_genes\tnum_unique_genes\n')
for key in keys:
    f.write(key + '\t')
    f.write(str(len(all_rankings[key])) + '\t')
    num_unique_genes = len(set([x[0] for x in all_rankings[key]]))
    f.write(str(num_unique_genes) + '\n')

    fk = open('check_paralogs_out_' + key + '.tsv', 'w')
    fk.write('gene\tregion\tstrain\t' + '\t'.join(keys[:-1]) + '\n')
    for item in all_rankings[key]:
        fk.write('\t'.join(item) + '\n')
    fk.close()
    
f.close()


