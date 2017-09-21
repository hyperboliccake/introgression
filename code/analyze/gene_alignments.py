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

def get_gene_coords(gene, fn):

    f = open(fn, 'r')

    line = f.readline()
    eof = False
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
                f.close()
                return start, end

    f.close()    

gp_dir = '../'
tag = sys.argv[1]
fn = gp.analysis_out_dir_absolute + 'introgressed_hmm_' + tag + '_genes.txt'
regions = process_helpers.read_regions_with_genes(fn)
regions = summary_stats_helpers.remove_duplicates(regions)

line_length = 80

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
            print gene
            fn = gp.gb_master_dir + gp.master_ref + '/' + gp.master_ref + \
                '_chr' + chrm + '.gb'
            gene_start, gene_end = get_gene_coords(gene, fn)
            fn = gp_dir + gp.alignments_dir + '_'.join(gp.alignment_ref_order) + \
                '_' + strain + '_chr' + chrm + '.maf'
            block = read_maf.read_mugsy_block(gene_block_labels[gene], fn)
            seq_master_ref = block['strains'][gp.master_ref]['sequence']
            alignment_block_start = block['strains'][gp.master_ref]['start']
            i_start = process_helpers.index_ignoring_gaps(\
                seq_master_ref, gene_start, alignment_block_start)
            i_end = process_helpers.index_ignoring_gaps(\
                seq_master_ref, gene_end, alignment_block_start)
            if i_start == -1:
                i_start = 0
            if i_end == len(seq_master_ref):
                i_end = len(seq_master_ref) - 1

            gene_seqs = {}
            for s in gp.alignment_ref_order + [strain]:
                gene_seqs[s] = block['strains'][s]['sequence'][i_start:i_end+1]
                assert len(gene_seqs[s]) == i_end - i_start + 1, \
                    str(len(gene_seqs[s])) + ' ' + \
                    str(i_start) + ' ' + str(i_end) + ' ' + \
                    str(block['strains'][gp.master_ref]['start']) + ' ' + \
                    str(block['strains'][gp.master_ref]['length'])

            matches = []
            for i in range(i_end - i_start + 1):
                # TODO is there a way to generalize this?
                if gene_seqs[strain][i] == gp.gap_symbol or \
                        gene_seqs[gp.alignment_ref_order[0]][i] == gp.gap_symbol or \
                        gene_seqs[gp.alignment_ref_order[1]][i] == gp.gap_symbol:
                    matches.append('-')
                elif gene_seqs[strain][i] == gene_seqs[gp.alignment_ref_order[0]][i]:
                    if gene_seqs[strain][i] == gene_seqs[gp.alignment_ref_order[1]][i]:
                        matches.append(' ')
                    else:
                        matches.append('C')
                else:
                    if gene_seqs[strain][i] == gene_seqs[gp.alignment_ref_order[1]][i]:
                        matches.append('P')
                    else:
                        matches.append('N')

            # this length includes gaps
            annotations = [' '] * (i_end - i_start + 1)
            for entry in genes[gene]:
                region_start = process_helpers.index_ignoring_gaps(\
                    gene_seqs[gp.master_ref], entry['region_start'], gene_start)
                region_end = process_helpers.index_ignoring_gaps(\
                    gene_seqs[gp.master_ref], entry['region_end'], gene_start)

                if region_start == -1:
                    region_start = 0
                if region_end == len(gene_seqs[gp.master_ref]):
                    region_end = len(gene_seqs[gp.master_ref]) - 1

                l = region_end - region_start + 1
                annotations[region_start:region_end+1] = ['i'] * l

            fn_gene = gp.genes_out_dir_absolute + gene + '/' + \
                tag + '/' + gene + '_' + \
                '_'.join(gp.alignment_ref_order) + '_' + strain + '.txt'
            if not os.path.exists(os.path.dirname(fn_gene)):
                try:
                    os.makedirs(os.path.dirname(fn_gene))
                except:
                    pass
            f_gene = open(fn_gene, 'w')

            l = 0
            while l < len(matches):
                for s in gp.alignment_ref_order + [strain]:
                    f_gene.write(gene_seqs[s][l:l+line_length] + '\n')
                f_gene.write(''.join(matches[l:l+line_length]) + '\n')            
                f_gene.write(''.join(annotations[l:l+line_length]) + '\n')
                f_gene.write('\n')
                l += line_length

            f_gene.close()
