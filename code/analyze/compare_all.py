# generate a file in ../../results/gene_alignments/ for each
# introgressed gene, which contains one threeway alignment for each
# strain in which the gene was called introgressed; the format of
# this file is:
# one line for each of the three strains in the alignment. each of
# these lines contains the strain name, chromosome # (in roman numerals),
# start position, end position, gene name
# next there is a blank line, followed by the alignment.
# the alignment is broken across lines for readability, and there is a
# fourth row indicating whether each position matches the cerevisaie ('c') or
# paradoxus ('p') reference or both (' ') or neither ('-').
# alignments for other strains, if any, follow in the same format after
# a blank line

import sys

lines = [x.split(',') for x in open('../../data/Table_S5_introgressed_genes.csv', 'r').readlines()]
genes = []
genes_verified = []
for i in range(2, len(lines)):
    genes.append(lines[i][2])
    if lines[i][4] == 'Verified':
        genes_verified.append(lines[i][2])

lines = [x.strip().split(' ') for x in open('../../results/introgressed_id_genes.txt', 'r').readlines()]
my_genes = [x[0] for x in lines]
my_genes_strains = {}
for line in lines:
    my_genes_strains[line[0]] = line[1:]

pm = []
pnm = []
npm = []
for g in genes:
    if g in my_genes:
        pm.append(g)
    else:
        pnm.append(g)
for g in my_genes:
    if g not in genes:
        npm.append(g)

outf = open('../../results/genes_1.txt', 'w')
#outf.write('genes I identify that Strope et al. didn\'t:\n')
for x in sorted(npm):
    outf.write(x)
    for s in my_genes_strains[x]:
        outf.write(' ' + s)
    outf.write('\n')
outf.close()

outf = open('../../results/genes_2.txt', 'w')
#outf.write('\ngenes I identify that Strope et al. also did:\n')
for x in sorted(pm):
    outf.write(x)
    for s in my_genes_strains[x]:
        outf.write(' ' + s)
    outf.write('\n')
outf.close()

outf = open('../../results/genes_3.txt', 'w')
#outf.write('\ngenes Strope et al. identified that I don\'t:\n')
for x in sorted(pnm):
    if x != '':
        outf.write(x + '\n')
outf.close()


lines = [x.strip().split(' ') for x in open('../../results/introgressed_id_genes_fns.txt', 'r').readlines()]
gene_to_fns = {}
for line in lines:
    gene_to_fns[line[0]] = line[1:]

for gene in gene_to_fns:
    outf = open('../../results/gene_alignments/' + gene + '_alignments.txt', 'w')    
    for fn in gene_to_fns[gene]:
        f = open(fn, 'r')

        print fn[:-(len('_annotated.maf'))] + '_alignments.txt'
        hc = f.readline()
        seqc = f.readline().strip()
        hp = f.readline()
        seqp = f.readline().strip()
        hx = f.readline()
        seqx = f.readline().strip()

        seq = ''
        for i in range(len(seqx)):
            if seqx[i] == '|':
                seq += '|'
            elif seqc[i] == seqx[i]:
                if seqp[i] == seqx[i]:
                    seq += ' '
                else:
                    seq += 'c'
            elif seqp[i] == seqx[i]:
                seq += 'p'
            else:
                seq += '-'
        line_length = 80
        m = fn[fn.rfind('/'):].split('_')
        outf.write(hc[1:])
        outf.write(hp[1:])
        outf.write(hx[1:] + '\n')
        for i in range(0, len(seq), line_length):
            outf.write(seqc[i:i+line_length] + '\n')
            outf.write(seqp[i:i+line_length] + '\n')
            outf.write(seqx[i:i+line_length] + '\n')
            outf.write(seq[i:i+line_length] + '\n\n')

    outf.close()
# TODO: get alignments for genes found in paper but not by me; print positions in each genome before alignments
