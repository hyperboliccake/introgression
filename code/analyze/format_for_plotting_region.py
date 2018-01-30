import gene_predictions
import sys
import os
import gzip
sys.path.insert(0, '../misc/')
import read_fasta
import global_params as gp
sys.path.insert(0, '../sim/')

def read_annotated_alignment(fn, nstrains):
    f = gzip.open(fn, 'rb')
    lines = f.readlines()
    f.close()
    strains = [l[:-1] for l in lines[:nstrains]]
    genes = lines[nstrains + 2][len('genes:'):-1].split()
    
    x = 11
    match_cer = ''
    match_par = ''
    gene = ''
    gene_ind = -1
    intd = ''

    while x < len(lines):
        cline = lines[x][:-1]
        match_cer += cline

        pline = lines[x+1][:-1]
        match_par += pline

        gline = lines[x+2][:-1]
        gene += gline

        iline = lines[x+3][:-1]
        intd += iline

        x += 8

    return match_cer, match_par, gene, genes, intd

def write_ps_annotated(match_cer, match_par, gene, glist, intd, region, fn):

    f = open(fn, 'w')
    f.write('ps\tmatch\tintd\tgene\n')

    block_start = int(region['start']) - intd.index('I') 
    block_end = len(intd) - intd.rindex('I') + int(region['end'])

    out_of_gene = True
    gene_ind = -1
    for i in range(len(match_cer)):
        f.write(str(block_start + i) + '\t')
        f.write(match_cer[i] + match_par[i] + '\t')
        f.write(intd[i] + '\t')
        g = ''
        if gene[i] == '=':
            if out_of_gene:
                gene_ind += 1
                out_of_gene = False
            g = glist[gene_ind]
        else:
            out_of_gene = True
        f.write(g)
        f.write('\n')
    f.close()

tag = sys.argv[1]
region = sys.argv[2]

blocks_fn = gp.analysis_out_dir_absolute + tag + '/' + \
            'introgressed_blocks_' + 'par' + '_' + tag + '_summary.txt'
r = gene_predictions.read_region_summary(blocks_fn)
strain = r[region]['strain']
chrm = r[region]['chromosome']


align_fn = gp.analysis_out_dir_absolute + tag + '/' + \
           'regions/' + region + '_annotated.txt.gz'
match_cer, match_par, gene, glist, intd = read_annotated_alignment(align_fn, 3)


plot_dir = gp.analysis_out_dir_absolute + '/' + tag + '/plots/' + region + '/'
if not os.path.isdir(plot_dir):
    os.makedirs(plot_dir)
fn_out = plot_dir + 'ps_annotations.txt'

write_ps_annotated(match_cer, match_par, gene, glist, intd, r[region], fn_out)

#probs_f = gzip.open(gp.analysis_out_dir_absolute + tag + '/' + \
#                    'probs_' + tag + '.txt.gz', 'rb')

