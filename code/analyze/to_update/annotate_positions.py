import sys
import re
import gzip
sys.path.insert(0, '../misc/')
import overlap
import read_fasta

def get_genes(fn):

    f = open(fn, 'r')
    genes = {}
    for line in f.readlines():
        line = line[:-1].split('\t')
        genes[(int(line[1]), int(line[2]))] = line[0]
    f.close()
    return genes

def get_orfs(fn):
    headers, seqs = read_fasta.read_fasta(fn)
    orfs = {}
    for h in headers:
        m = re.search(' (?P<name>[a-zA-Z0-9]+)_(?P<strain>[a-zA-Z0-9\.]+):(?P<start>[0-9]+):(?P<end>[0-9]+)', h)
        orfs[(int(m.group('start')), int(m.group('end')))] = m.group('name')
    return orfs

def write_annotated_file(coords, genes, orfs, fn):
    # could definitely do this all way more efficiently
    sep = '\t'
    f = gzip.open(fn, 'wb')
    f.write('position' + sep + 'ref_position' + sep + 'gene' + sep + 'ORF\n')
    for i in range(len(coords)):
        f.write(str(i) + sep)
        if int(coords[i]) == coords[i]:
            f.write(str(int(coords[i])) + sep)
            gene = overlap.contained_any_named(coords[i], genes)
            if gene != None:
                f.write(gene)
            f.write(sep)
        else:
            f.write(str(coords[i]) + sep + sep)
        orf = overlap.contained_any_named(i, orfs)
        if orf != None:
            f.write(orf)
        f.write('\n')
    f.close()

