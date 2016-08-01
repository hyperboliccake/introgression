# make site frequency spectrum
# histogram of sequences found in 1, 2, ..., n strains
# - genes
# - individual sites (indexed by cer reference)

# before running this, need to have run process.py to generate files
# containing alignments for all introgressed regions, as well as file
# of which genes are introgressed in which strains

# then, given these, easy to make histogram of introgressed genes

# for individual sites, run through all introgressed region files and
# keep track of all cerevisiae indices introgressed in one or more
# strains

#####
# sfs of genes
#####

f_genes = open('../../results/introgressed_id_genes.txt', 'r')
lines = [x.split() for x in f_genes.readline().strip()]
f_genes.close()
gene_hist = {}
for l in lines:
    x = len(l[1:])
    if x not in gene_hist:
        gene_hist[x] = 0
    gene_hist[x] += 1
f = open('../../results/sfs_genes.txt', 'w')
k = gene_hist.keys()
for i in range(1, max(k)):
    f.write(str(i) + ' ')
    if i in k:
        f.write(str(gene_hist[i]))
    else:
        f.write('0')
    f.write('\n')
f.close()

#####
# sfs of individual sites
#####

prefix = '../../results/regions/'
fns = [prefix + os.listdir(prefix)]
chrms_roman = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XIV']
site_counts = dict(zip(chrms_roman, [{} for x in chrms_roman]))
for fn in fns:
    f = open(fn, 'r')
    line = f.readline()[1:].split()
    assert line[0] == 'S288c'
    f.close()
    chrm, start, end = line[1:]
    assert chrm in chrms_roman
    line = f.readline()
    # inclusive lower bound
    i = 0
    while line[i] != '|':
        if line[i] != '-':
            start += 1
        i += 1
    line = line[i+1:line.rfind('|')]
    # exclusive upper bound
    i = 0
    end = start
    while line[i] != '|':
        if line[i] != '-':
            end += 1
        i += 1
    for i in range(start, end):
        if i not in site_counts[chrm]:
            site_counts[chrm][i] = 0
        site_counts[chrm][i] += 1
    
site_hist = {}
for chrm in site_counts:
    for site in site_counts[chrm]:
        count = site_counts[chrm][site]
        if count not in site_hist:
            site_hist[count] = 0
        site_hist[count] += 1
f = open('../../results/sfs_genes.txt', 'w')
k = site_hist.keys()
for i in range(1, max(k)):
    f.write(str(i) + ' ')
    if i in k:
        f.write(str(site_hist[i]))
    else:
        f.write('0')
    f.write('\n')
f.close()
    


