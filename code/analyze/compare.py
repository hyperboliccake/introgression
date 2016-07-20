lines = [x.split(',') for x in open('../Table_S5_introgressed_genes.csv', 'r').readlines()]
genes = []
genes_verified = []
for i in range(2, len(lines)):
    genes.append(lines[i][2])
    if lines[i][4] == 'Verified':
        genes_verified.append(lines[i][2])

lines = [x.split(' ') for x in open('../../results/introgressed_id_genes.txt', 'r').readlines()]
my_genes = [x[0] for x in lines]

print len(genes), 'genes from paper'
print len(genes_verified), 'verified genes from paper'
print len(my_genes), '(verified) genes I identify'

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

print 'genes found in paper that I found (', len(pm), '):'
for x in pm:
    print x
print 'genes found in paper that I didn\'t find (', len(pnm), '):'
for x in pnm:
    print x
print 'genes that I found not in paper(', len(npm), '):'
for x in npm:
    print x

lines = [x.strip().split(' ') for x in open('../../results/introgressed_id_genes_fns.txt', 'r').readlines()]
gene_to_fns = {}
for line in lines:
    gene_to_fns[line[0]] = line[1:]

while True:
    gene = raw_input('=========================================\nwhich gene? ')
    try:
        gene_to_fns[gene]
    except:
        print 'that gene wasn\'t one i found'
        continue
    for fn in gene_to_fns[gene]:
        f = open(fn)

        f.readline()
        seqc = f.readline().strip()
        f.readline()
        seqp = f.readline().strip()
        f.readline()
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
        print '==========', fn
        line_length = 10000
        for i in range(0, len(seq), line_length):
            print seqc[i:i+line_length]
            print seqp[i:i+line_length]
            print seqx[i:i+line_length]
            print seq[i:i+line_length]
            print
    print
    raw_input('')


# TODO: get alignments for genes found in paper but not by me; print positions in each genome before alignments
