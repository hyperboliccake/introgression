# Given a list of introgressed regions in this format:
# Sigma1278b.chrX, + strand, regionStart-regionEnd, blockStart, blockEnd
#
# generate a set of annotations in the ../../results/analyze/ folder:
#
# for each region called introgressed, generate a file in
# results/analyze/regions/ that is the alignment for the two references and
# the strain with the introgressed region,
# S288c_CBS432_strain_chr_start-end.fasta. In this file, the aligned
# bases within coding sequence are upper case. In addition, there is a
# corresponding file S288c_CBS432_strain_chr_start-end.genes.txt
# listing the genes that overlap with this region, and the indices of
# the bases they overlap, in this format: 
# gene_name, 0-149, 25236-25385 
# gene_name, 200-600, ....
# 
# also generate a file in results/gene_alignments/ for each
# introgressed gene, which contains one threeway alignment for each
# strain in which the gene was called introgressed...followed by all
# the genes that weren't
#
# todo in future:
# for each gene that is called introgressed in at least one strain,
# create folder gene/ containing the alignment
# of cerevisiae and paradoxus references to all the introgressed
# versions (gene_introgressed.fasta), and also to all of the versions
# (gene_all.fasta).


import re
import sys
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc/')
import read_maf

def reverse_complement(s):
    t = ''
    for i in s:
        if i == 'A':
            t += 'T'
        elif i == 'T':
            t += 'A'
        elif i == 'G':
            t += 'C'
        elif i == 'C':
            t += 'G'
        else:
            t += i
    return t[::-1]

translate = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
       "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
       "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
       "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
       "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
       "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
       "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
       "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
       "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
       "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
       "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

def starts_with_any(s, l):
    for item in l:
        if s.startswith(item):
            return True
    return False

def read_regions(fn):
    regions = {} 
    f = open(fn, 'r')
    line = f.readline()
    i = 0
    while line != '':
        line = [x.strip() for x in line.split(',')]
        strain = line[0][:line[0].find('_')].lower()
        chrm = re.search(r'chr(?P<chrm>[IVXM]+)', line[0]).group('chrm')

        # note that end is INCLUSIVE
        entry = {'strand':line[1][0], \
                     'start':int(line[2]), \
                     'end':int(line[3]), \
                     'block_start':int(line[4]), \
                     'block_length':int(line[5]), \
                     'id':'r' + str(i)}
        if strain not in regions:
            regions[strain] = {}
        if chrm not in regions[strain]:
            regions[strain][chrm] = []
        regions[strain][chrm].append(entry)
        line = f.readline()
        i += 1
    f.close()
    return regions

def index_gapped(seq, i):
    '''i is the index in the non-gapped sequence; this function
    returns the index in (gapped) seq that corresponds to i'''
    # if i < 0, return index of first non-gap character
    if i < 0:
        i = 0
    count_non_gaps = 0
    last_non_gap = -1
    x = 0
    while count_non_gaps <= i and x < len(seq):
        if seq[x] != '-':
            count_non_gaps += 1
            last_non_gap = x
        x += 1
    # if i is greater than the number of non-gapped characters, return
    # the index of the last non-gap character
    if x == len(seq):
        return last_non_gap
    return x - 1

def write_region_alignment(block, strain, entry, fn, context = 200):
    '''write relevant portion of alignment block'''

    # entry is for one introgressed region, and will get its own file
    # containing the alignment in the introgressed region plus some
    # context on either side

    relative_start_ungapped = entry['start'] - block['strains'][strain]['start']
    relative_start = index_gapped(block['strains'][strain]['sequence'], \
                                      relative_start_ungapped)
    relative_start_with_context = index_gapped(block['strains'][strain]['sequence'], \
                                                   relative_start_ungapped - context)


    relative_end_ungapped = entry['end'] - block['strains'][strain]['start'] + 1
    relative_end = index_gapped(block['strains'][strain]['sequence'], \
                                    relative_end_ungapped)
    relative_end_with_context = index_gapped(block['strains'][strain]['sequence'], \
                                                 relative_end_ungapped + context)

    # want to be able to look up all the regions in this alignment block later
    if 'regions' not in block:
        block['regions'] = {}
    block['regions'][entry['id']] = { \
        'relative start': relative_start,\
            'relative start with context': relative_start_with_context,\
            'relative end': relative_end,\
            'relative end with context': relative_end_with_context}
    
    f = open(fn, 'w')
    all_strains = block['strains'].keys()
    all_strains.sort(key = lambda x: block['strains'][x]['index'])
    for current_strain in all_strains:
        # TODO put in positions? and make sure to correct based on strain
        f.write('>' + current_strain + '\n')
        # save annotation of introgressed region for later
        #f.write(block['strains'][current_strain]['sequence'][relative_start_with_context:relative_start] + '{')
        #f.write(block['strains'][current_strain]['sequence'][relative_start:relative_end] + '}')
        #f.write(block['strains'][current_strain]['sequence'][relative_end:relative_end_with_context] + '\n')
        f.write(block['strains'][current_strain]['sequence'][relative_start:relative_end].lower() + '\n')

    f.close()
    # return the modified block
    return block

def write_annotated_region_alignment(block, strain, entry, fn, genes, introgressed_genes_strains, strains_introgressed_genes):

    all_region_inds = []
    for region_id in block['regions']:
        inds = block['regions'][region_id]
        all_region_inds.append((inds['relative start'], inds['relative end']))
    all_region_inds.sort(key=lambda x: x[0], reverse=True)

    inds = block['regions'][entry['id']]

    # find relative gene_indices
    block_start = block['strains'][strain]['start']
    block_end = block_start + block['strains'][strain]['length']
    seq = block['strains'][strain]['sequence']
    all_gene_inds = []
    for gene_name in genes:
        gene_start, gene_end = genes[gene_name]
        # add 1 because gene end (gene[1]) is INCLUSIVE
        gene_end += 1
        # all genes that fall within this block
        if gene_start >= block_start and gene_end <= block_end:
            s = index_gapped(seq, gene_start - block_start)
            if gene_end >= block_start and gene_end <= block_end:
                e = index_gapped(seq, gene_end - block_start)
                all_gene_inds.append((s, e))
            else:
                all_gene_inds.append((s, len(seq)))
        elif gene_end >= block_start and gene_end <= block_end:
            e = index_gapped(seq, gene_end - block_start)
            all_gene_inds.append((0, e))

        # keep track of genes that fall in introgressed region (entry)
        if (gene_start >= entry['start'] and gene_start <= entry['end']) or \
                (gene_end - 1 >= entry['start'] and gene_end - 1 <= entry['end']):
            if gene_name not in introgressed_genes_strains:
                introgressed_genes_strains[gene_name] = []
            introgressed_genes_strains[gene_name].append(strain)
            if strain not in strains_introgressed_genes:
                strains_introgressed_genes[strain] = []
            strains_introgressed_genes[strain].append(gene_name)
            if 'genes' not in entry:
                entry['genes'] = []
            entry['genes'].append(gene_name)

    # annotation time
    f = open(fn, 'w')
    all_strains = block['strains'].keys()
    all_strains.sort(key = lambda x: block['strains'][x]['index'])
    for current_strain in all_strains:
        # TODO put in positions? and make sure to correct based on strain
        f.write('>' + current_strain + '\n')
        seq = block['strains'][current_strain]['sequence']
        seq = seq.lower()
        
        # first, capitalize all bp that fall within a gene
        for gene_start, gene_end in all_gene_inds:
            seq = seq[:gene_start] + \
                seq[gene_start:gene_end].upper() + \
                seq[gene_end:]

        # then add in open and close {} for all introgressed regions
        # that fall within the block we're writing
        lower_bound = inds['relative start with context']
        upper_bound = inds['relative end with context']
        new_seq = ''
        p = lower_bound
        for relative_start, relative_end in all_region_inds:
            if relative_start >= p and relative_start <= upper_bound:
                left = seq[p:relative_start]
                new_seq += left + '<'
                p += len(left)
            if relative_end >= p and relative_end <= upper_bound:
                middle = seq[p:relative_end]
                new_seq += middle + '>'
                p += len(middle)
        new_seq += seq[p:upper_bound]
            
        # finally write the sequence
        f.write(new_seq + '\n')

    f.close()

    return entry

#####
# read in introgressed regions
#####

gp_dir = '../'
fn = gp_dir + gp.analysis_out_dir + 'introgressed.txt'
# introgressed regions keyed by strain and then chromosome
regions = read_regions(fn)

##### 
# For each introgressed region, extract relevant part of alignment to a separate file
#####
alignment_blocks = {}
for strain in regions.keys():
    alignment_blocks[strain] = {}
    print strain
    for chrm in regions[strain]:
        print '-', chrm
        # file for alignment of this region
        fn_align = gp.alignments_dir + \
            gp.cer_ref_strain + '_' + gp.par_ref_strain + '_' + \
            strain + '_chr' + chrm + gp.alignment_suffix
        alignment_blocks[strain][chrm] = read_maf.read_maf(fn_align)

        for ei in range(len(regions[strain][chrm])):
            entry = regions[strain][chrm][ei]

            for label in alignment_blocks[strain][chrm]:
                # TODO change this to check label instead of start
                # NOTE: mugsy indexes from zero!
                block = alignment_blocks[strain][chrm][label]
                if strain in block['strains'] and \
                        block['strains'][strain]['start'] == entry['block_start']:

                    entry['block_label'] = label

                    fn_region = gp.regions_out_dir + \
                        gp.cer_ref_strain + '_' + gp.par_ref_strain + '_' + \
                        strain + '_chr' + chrm + '_' + \
                        str(entry['start']) + '-' + str(entry['end']) + \
                        gp.alignment_suffix

                    block_modified = write_region_alignment(block, strain, entry, fn_region)
                    alignment_blocks[strain][chrm][label] = block_modified
                    regions[strain][chrm][ei] = entry

#####
# Identify genes (if any) that overlap each alignment block and
# annotate the alignment file for each region to note where the
# introgressed boundaries and genes are
#####

introgressed_genes_strains = {} 
strains_introgressed_genes = {}

f = open(gp.gb_all, 'r')

# start by finding first species x chromosome
eof = False
line = f.readline()
while not line.strip().startswith('DEFINITION'):
    if line == '':
        eof = True
        break
    line = f.readline()

while not eof:
    # starting with a new species x chromosome
    assert line.strip().startswith('DEFINITION'), line
    done_with_chrm = False
    print line

    # TODO make this more general?
    #m = re.search('Saccharomyces cerevisiae (?P<strain>[a-zA-Z0-9]+) chromosome (?P<chrm>[IVXM]+)', line)
    m = re.search(' (?P<strain>[a-zA-Z0-9]+) chromosome (?P<chrm>[IVXM]+)', line)
    strain = m.group('strain').lower()
    chrm = m.group('chrm')

    # collect all the genes for this species x chromosome before moving on to next
    genes = {}
    
    # find next gene and process it, repeatedly until we hit a new
    # species x chromosome or the end of the file
    while not eof and not done_with_chrm:

        # start by finding first gene
        line = f.readline()
        while not line.strip().startswith('gene'): 
            # these two cases just here in case there are no genes
            if line == '':
                eof = True
                break
            if line.strip().startswith('DEFINITION'):
                done_with_chrm = True
                break
            line = f.readline()
        if eof or done_with_chrm:
            break

        # starting with new gene
        assert line.strip().startswith('gene'), line
        skip_this_gene = False

        # regex for finding coordinates
        m = re.search(r'[><]?(?P<start>[0-9]+)[.><,0-9]*\.\.[><]?(?P<end>[0-9]+)', line)
        # subtract one to index from zero TODO is this correct?
        start = int(m.group('start')) - 1
        end = int(m.group('end')) - 1

        # look for the name of the gene in the lines following the
        # start of the entry
        line = f.readline()
        while not line.strip().startswith('/gene'):
            # sometimes we never run into a gene name for whatever
            # reason, and in that case, we'll just skip over this
            # entry that we found coordinates for
            if line == '':
                eof = True
                break
            if line.strip().startswith('gene'):
                skip_this_gene = True
                break
            if line.strip().startswith('DEFINITION'):
                done_with_chrm = True
                break
            line = f.readline()

        if not skip_this_gene:
            gene_name = line[line.find('/gene="')+7:-2]
            genes[gene_name] = (start, end)

    # now that we have all the genes, proceed with annotations (if
    # there's anything to annotate)
    if strain not in regions or chrm not in regions[strain]:
        continue

    for ei in range(len(regions[strain][chrm])):
        entry = regions[strain][chrm][ei]
        fn_region_annotated = gp.regions_out_dir + \
            gp.cer_ref_strain + '_' + gp.par_ref_strain + '_' + \
            strain + '_chr' + chrm + '_' + \
            str(entry['start']) + '-' + str(entry['end']) + \
            '_annotated' + \
            gp.alignment_suffix

        entry = write_annotated_region_alignment(alignment_blocks[strain][chrm][entry['block_label']], \
                                                     strain, entry, fn_region_annotated, genes, \
                                                     introgressed_genes_strains, strains_introgressed_genes)
        regions[strain][chrm][ei] = entry
#####
# gene x strain and strain x gene files
#####

f = open(gp.analysis_out_dir + '/introgressed_genes.txt', 'w')
genes_sorted = []
for gene in introgressed_genes_strains:
    genes_sorted.append((gene, len(introgressed_genes_strains[gene])))
genes_sorted.sort(key=lambda x: x[1])
for gene, n in genes_sorted:
    f.write(gene + ' ' + str(n))
    for strain in introgressed_genes_strains[gene]:
        f.write(' ' + strain)
    f.write('\n')
f.close()


f = open(gp.analysis_out_dir + '/introgressed_strains.txt', 'w')
strains_sorted = []
for strain in strains_introgressed_genes:
    strains_sorted.append((strain, len(strains_introgressed_genes[strain])))
strains_sorted.sort(key=lambda x: x[1])
for strain, n in strains_sorted:
    f.write(strain + ' ' + str(n))
    for gene in strains_introgressed_genes[strain]:
        f.write(' ' + gene)
    f.write('\n')
f.close()
    

#####
# percentage of introgressed regions due to gaps in reference
#####

f = open(gp.analysis_out_dir + '/introgressed_annotated.txt', 'w')
for strain in regions:
    for chrm in regions[strain]:
        for entry in regions[strain][chrm]:
            block = alignment_blocks[strain][chrm][entry['block_label']]
            relative_start = block['regions'][entry['id']]['relative start']
            relative_end = block['regions'][entry['id']]['relative end']
            # sequence for the reference that is "non-introgressed"
            gap_fraction = -1 # if that reference not aligned in this region
            if gp.cer_ref_strain in block['strains']:
                ref_seq = block['strains'][gp.cer_ref_strain]['sequence']
                gap_count = ref_seq[relative_start:relative_end+1].count('-')
                gap_fraction = float(gap_count) / (relative_end - relative_start + 1)
            f.write(strain + ',' + chrm + ',' + entry['strand'] + ',')
            f.write(str(entry['start']) + ',' + str(entry['end']) + ',' + str(entry['end'] - entry['start'] + 1) + ',')
            f.write(str(gap_fraction) + ',')
            #for gene in entry['genes']:
            #    f.write(gene + ' ')
            f.write('\n')
f.close()




