import sys
import os
import re
import process_helpers
import summary_stats_helpers
sys.path.insert(0, '../')
import global_params as gp
sys.path.insert(0, '../misc/')
import read_maf

def pad(s, n, c = ' '):
    end = min(len(s), n)
    return s[:end] + (n - end) * c

def get_gene_coords(gene, fn, strain, chrm, get_seq=False):

    f = open(fn, 'r')

    line = f.readline()
    while True:
        if line.startswith('DEFINITION'):
            m = re.search(' (?P<strain>[a-zA-Z0-9]+) chromosome (?P<chrm>[IVXM]+)', line)
            current_strain = m.group('strain')
            current_chrm = m.group('chrm')
            if (current_strain == strain or current_strain.lower() == strain) \
                    and current_chrm == chrm:
                break

        line = f.readline()

    eof = False
    start = -1
    end = -1
    while line != '':

        while not line.startswith('     gene'):
            if line == '':
                eof = True
                break
            line = f.readline()
        if eof:
            break

        # starting with new gene
        #assert line.strip().startswith('gene'), line
        skip_this_gene = False

        # regex for finding coordinates
        m = re.search(r'[><]?(?P<start>[0-9]+)[.><,0-9]*\.\.[><]?(?P<end>[0-9]+)', line)
        # subtract one to index from zero TODO is this correct? end is
        # inclusive

        start = int(m.group('start')) - 1
        end = int(m.group('end')) - 1

        # look for the name of the gene in the lines following the
        # start of the entry
        line = f.readline()
        while not line.strip().startswith('/gene="'):
            # sometimes we never run into a gene name for whatever
            # reason, and in that case, we'll just skip over this
            # entry that we found coordinates for
            if line == '':
                eof = True
                break
            if line.startswith('     gene'):
                skip_this_gene = True
                break
            line = f.readline()

        if not skip_this_gene:
            gene_name = line[line.find('/gene="')+7:-2]
            if gene_name == gene:
                break

    if not get_seq:
        f.close()    
        return start, end

    assert not eof

    while not line.startswith('ORIGIN'):
        line = f.readline()

    seq = ''
    line = f.readline()
    while not line.startswith('//'):
        seq += ''.join(line.strip().split()[1:])
        line = f.readline()

    return start, end, seq[start:end+1]

gp_dir = '../'
tag = sys.argv[1]
fn = gp.analysis_out_dir_absolute + 'introgressed_hmm_' + tag + '_genes.txt'
regions = process_helpers.read_regions_with_genes(fn)
regions = summary_stats_helpers.remove_duplicates(regions)

line_length = 80

goi = ['LCB3']

for strain in regions:
    for chrm in regions[strain]:
        # a gene might show up multiple times on this chromosome but
        # we only need to deal with it once, but keep track of all
        # introgressed portions of it
        genes = {}
        gene_block_labels = {}
        label = None
        for entry in regions[strain][chrm]:
            for gene in entry['genes']:
                if gene not in genes:
                    genes[gene] = []
                genes[gene].append(entry)
                gene_block_labels[gene] = entry['block_label']

        for gene in genes:
            if gene not in goi:
                continue
            
            print gene, strain, chrm

            gene_seqs = {}
            for ref in gp.alignment_ref_order[:-1]:
                fn = gp.gb_master_dir + ref + '/' + ref + '_chr' + chrm + '.gb'
                start, end, seq = get_gene_coords(gene, fn, ref, chrm, get_seq=True)
                gene_seqs[ref] = seq

            fn = gp.gb_all
            start, end, seq = get_gene_coords(gene, fn, strain, chrm, get_seq=True)
            gene_seqs[strain] = seq

            fn_gene = gp.genes_out_dir_absolute + gene + '/' + \
                tag + '/' + gene + '_' + \
                '_'.join(gp.alignment_ref_order[:-1]) + '_' + strain + '_unaligned.txt'
            if not os.path.exists(os.path.dirname(fn_gene)):
                try:
                    os.makedirs(os.path.dirname(fn_gene))
                except:
                    pass
            f_gene = open(fn_gene, 'w')

            for s in gp.alignment_ref_order[:-1] + [strain]:
                f_gene.write(gene_seqs[s] + '\n')

            f_gene.close()
