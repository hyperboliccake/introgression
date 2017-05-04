import re
import sys
import copy
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../sim/')
import sim_analyze_hmm_bw as sim
sys.path.insert(0, '../misc/')
import read_maf

def read_regions(fn):
    regions = {} 
    f = open(fn, 'r')
    f.readline() # header
    line = f.readline()
    i = 0
    while line != '':
        line = line.strip().split('\t')
        if len(line) < 8:
            break
        strain = line[0]
        chrm = line[1]
        # alignment block label, strand, predicted reference, region
        # start, region end, number non gap sites
        alignment_block_label, strand, predicted_reference, \
            region_start, region_end, number_non_gap_sites = line[2:]
        entry = {}
        entry['block_label'] = alignment_block_label
        entry['strand'] = strand
        entry['predicted_reference'] = predicted_reference
        entry['region_start'] = int(region_start)
        entry['region_end'] = int(region_end)
        entry['number_non_gap'] = int(number_non_gap_sites)
        entry['region_id'] = 'region_' + str(i)
        #entry = alignment_block_label, strand, predicted_reference, \
        #    int(region_start), int(region_end), int(number_non_gap_sites)
        if strain not in regions:
            regions[strain] = {}
        if chrm not in regions[strain]:
            regions[strain][chrm] = []
        regions[strain][chrm].append(entry)
        line = f.readline()
        i += 1
    f.close()
    return regions

def read_regions_with_genes(fn):
    regions = {} 
    f = open(fn, 'r')
    f.readline() # header
    line = f.readline()
    i = 0
    while line != '':
        line = line.strip().split('\t')
        if len(line) < 8:
            break
        strain = line[0]
        chrm = line[1]
        # alignment block label, strand, predicted reference, region
        # start, region end, number non gap sites
        entry = {}
        if len(line) == 9:
            genes = line[-1]
            entry['genes'] = genes.split(' ')
            assert '' not in entry['genes']
            line = line[:-1]
        else:
            entry['genes'] = []
        alignment_block_label, strand, predicted_reference, \
            region_start, region_end, number_non_gap_sites = line[2:]
        entry['block_label'] = alignment_block_label
        entry['strand'] = strand
        entry['predicted_reference'] = predicted_reference
        entry['region_start'] = int(region_start)
        entry['region_end'] = int(region_end)
        entry['number_non_gap'] = int(number_non_gap_sites)
        #entry = alignment_block_label, strand, predicted_reference, \
        #    int(region_start), int(region_end), int(number_non_gap_sites)
        if strain not in regions:
            regions[strain] = {}
        if chrm not in regions[strain]:
            regions[strain][chrm] = []
        regions[strain][chrm].append(entry)
        line = f.readline()
    f.close()
    return regions

def index_ignoring_gaps(s, i, s_start):
    '''returns the index of the ith non-gap character in s, given that
    the start index of s is s_start (instead of needing to be 0)'''

    #if mode == 'error' and (i < s_start or i >= s_start + s_length):
    #    return -1
    #if i < s_start:
    #    return 0
    #if i >= s_start + s_length:
    #    return s_start + s_length
    
    x = 0
    non_gap_count = 0
    i -= s_start
    if i < 0:
        return -1
    while x < len(s):
        if s[x] != gp.gap_symbol and non_gap_count >= i:
            return x
        if s[x] != gp.gap_symbol:
            non_gap_count += 1
        x += 1
    return x

def mark_gene(block, start, end, genes):
    seq = block['sequence']
    block_start = block['start']
    block_end = block_start + block['length'] - 1 # subtract 1 so index of end
    a = [''] * len(seq)
    len_a = len(a)

    # just annotate whole alignment block and then pull out portion at
    # end; TODO make this more efficient
    for gene_name in genes:
        gene_start, gene_end = genes[gene_name]
        gene_relative_start = index_ignoring_gaps(seq, gene_start, \
                                                      block_start)
        gene_relative_end = index_ignoring_gaps(seq, gene_end, \
                                                    block_start)

        # if gene start and end are both before the block, they will
        # both be -1; if gene start and end are both after the block,
        # they will both be len_a
        if (gene_relative_start == -1 and gene_relative_end == -1) or \
                (gene_relative_start == len_a and gene_relative_end == len_a):
            continue
        if gene_relative_start == -1:
            gene_relative_start = 0
        if gene_relative_end == -1:
            gene_relative_end = 0
        if gene_relative_start == len_a:
            gene_relative_start = len_a - 1
        if gene_relative_end == len_a:
            gene_relative_end = len_a - 1

        l = gene_relative_end - gene_relative_start + 1
        a[gene_relative_start:gene_relative_end + 1] = [gene_name] * l
        a = a[:len_a] # because the above can actually append

        assert len(a) == len_a, str(len(a)) + ' ' + str(len_a) + \
            ' ' + str(gene_relative_start) + ' ' + str(gene_relative_end)

    return a[start:end + 1]

def write_region_alignment(block, entry, genes, \
                               strain, master_ref, refs, \
                               fn, fn_annotated, context = 0):
    '''write relevant portion of alignment block'''

    # entry is for one introgressed region, and will get its own file
    # containing the alignment in the introgressed region, as well as
    # a similar file annotated with gene info

    seq_master_ref = block['strains'][master_ref]['sequence']
    length_master_ref = block['strains'][master_ref]['length']

    # start index of alignment block
    alignment_block_start = block['strains'][master_ref]['start']

    # relative indices into that block to pull out correct region:
    # - with context before
    i_start_with_context = index_ignoring_gaps(seq_master_ref, \
                                                   entry['region_start'] - context, \
                                                   alignment_block_start)
    # - without context before
    i_start = index_ignoring_gaps(seq_master_ref, \
                                      entry['region_start'], \
                                      alignment_block_start)
    # - without context after
    i_end = index_ignoring_gaps(seq_master_ref, \
                                    entry['region_end'], \
                                    alignment_block_start)
    # - with context after
    i_end_with_context = index_ignoring_gaps(seq_master_ref, \
                                                 entry['region_end'] + context, \
                                                 alignment_block_start)

    # TODO maybe make this less stupid
    if i_start_with_context == -1:
        i_start_with_context = 0
    if i_end_with_context == len(seq_master_ref):
        i_end_with_context = len(seq_master_ref) - 1

    # non-annotated file (and no context)
    f = open(fn, 'w')
    # start with references
    for ref in gp.alignment_ref_order:
        if ref in refs:
            f.write('>' + ref + '\n')
            seq = block['strains'][ref]['sequence'][i_start:i_end+1]
            f.write(seq.lower() + '\n')
    # then the current strain
    f.write('>' + strain + '\n')
    seq = block['strains'][strain]['sequence'][i_start:i_end+1]
    f.write(seq.lower() + '\n')
    f.close()

    # annotated file (with context)
    # mark each site as being in gene or not
    annotation = mark_gene(block['strains'][master_ref],
                           i_start_with_context, 
                           i_end_with_context,
                           genes)

    start_offset = i_start - i_start_with_context
    end_offset = i_end - i_start_with_context

    gene_set = [] # needs to be a list to preserve ordering
    gene_set_region = [] # not including context
    gene_introgressed_bases_count = {}
    gene_lengths = {}
    # loop through all positions and keep track of all unique genes
    # seen within the introgressed region, and also separately within
    # the region including surrounding context
    for i in range(len(annotation)):
        g = annotation[i]
        if g != '':
            if g not in gene_set:
                gene_set.append(g)
            if i >= start_offset and i <= end_offset:
                if g not in gene_set_region:
                    gene_set_region.append(g)
                    gene_lengths[g] = genes[g][1] - genes[g][0] + 1
                    gene_introgressed_bases_count[g] = 0
                # keep track of fraction of each gene that is
                # introgressed; this isn't perfect, but numerator is
                # number of sites within gene in introgressed block
                # (relative to master reference) and denominator is gene
                # length
                gene_introgressed_bases_count[g] += 1
        
    gene_set = filter(lambda x: x != '', gene_set)

    # now write to file
    f = open(fn_annotated, 'w')
    # first list the genes that are contained
    f.write('# ' + ' '.join(gene_set) + '\n')
    # start with references
    region_start_symbol = '<'
    region_end_symbol = '>'
    for ref in gp.alignment_ref_order:
        if ref in refs:
            f.write('>' + ref + '\n')
            seq = block['strains'][ref]['sequence']\
                [i_start_with_context:i_end_with_context+1].lower()
            assert len(seq) == len(annotation), str(len(seq)) + ' ' + \
                str(len(annotation)) + ' ' + str(i_start_with_context) + \
                ' ' + str(i_start) + ' ' + str(i_end) + ' ' + str(i_end_with_context)
            for i in range(len(seq)):
                if i == start_offset:
                    f.write(region_start_symbol)
                if annotation[i] != '':
                    f.write(seq[i].upper())
                else:
                    f.write(seq[i].lower())
                if i == end_offset:
                    f.write(region_end_symbol)
            f.write('\n')

    # then the current strain
    f.write('>' + strain + '\n')
    seq = block['strains'][strain]['sequence']\
        [i_start_with_context:i_end_with_context+1].lower()
    assert len(seq) == len(annotation), str(len(seq)) + ' ' + str(len(annotation)) + \
        ' ' + str(i_start_with_context) + \
        ' ' + str(i_start) + ' ' + str(i_end) + ' ' + str(i_end_with_context)
    for i in range(len(seq)):
        if i == start_offset:
            f.write(region_start_symbol)
        if annotation[i] != '':
            f.write(seq[i].upper())
        else:
            f.write(seq[i].lower())
        if i == end_offset:
            f.write(region_end_symbol)
    f.write('\n')
    f.close()    

    # store the genes we found in this region (NOT including stuff
    # only in context)
    entry['genes'] = gene_set_region

    # introgressed fractions of genes
    x = {}
    for gene in gene_lengths:
        x[gene] = float(gene_introgressed_bases_count[gene]) / gene_lengths[gene]
    entry['genes_introgressed_fractions'] = x

def write_region_alignment_old(block, entry, genes, \
                               strain, master_ref, refs, \
                               fn, fn_annotated, context = 0):
    '''write relevant portion of alignment block'''

    # entry is for one introgressed region, and will get its own file
    # containing the alignment in the introgressed region, as well as
    # a similar file annotated with gene info

    current_ind_master_ref = block['strains'][master_ref]['start']
    target_ind_master_ref = entry[3] - context
    num_non_gap_sites = entry[5]
    seq_master_ref = block['strains'][master_ref]['sequence']
    current_ind_block = 0
    # get up to start of region
    try:
        while True:
            if seq_master_ref[current_ind_block] != gp.gap_symbol and \
                    current_ind_master_ref < target_ind_master_ref:
                break
            current_ind_block += 1
            if seq_master_ref[current_ind_block] != gp.gap_symbol:
                current_ind_master_ref += 1
    except:
        print seq_master_ref, entry, block['strains'][master_ref]['start']
    block_start = current_ind_block
    # keep going until the end of the region
    target_ind_master_ref = entry[4]
    while True:
        if seq_master_ref[current_ind_block] != gp.gap_symbol and \
                current_ind_master_ref < target_ind_master_ref:
            break
        current_ind_block += 1
        if seq_master_ref[current_ind_block] != gp.gap_symbol:
            current_ind_master_ref += 1
    block_end = current_ind_block # inclusive end

    # write to file, first references
    f = open(fn, 'w')
    for ref in gp.alignment_ref_order:
        if ref in refs:
            f.write('>' + ref + '\n')
            seq = block['strains'][master_ref]['sequence'][block_start:block_end+1]
            f.write(seq.lower() + '\n')
    # then the current strain
    f.write('>' + strain + '\n')
    seq = block['strains'][strain]['sequence'][block_start:block_end+1]
    f.write(seq.lower() + '\n')

    f.close()

def read_genes(f):

    genes = {}
    line = f.readline()
    eof = False
    while line != '' and 'ORIGIN' not in line:
        # why yes i am relying on this probably arbitrary number of 5
        # spaces because fuck this file format (random other text can
        # start with 'gene')
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
            if gene_name != '':
                genes[gene_name] = (start, end)
            else:
                print 'gene name not found: ' + line

    return genes

def read_one_strain_chrm(f, line):

    eof = False
    while not line.strip().startswith('DEFINITION'):
        if line == '':
            eof = True
            break
        line = f.readline()

    if eof:
        return None

    done_with_chrm = False

    m = re.search(' (?P<strain>[a-zA-Z0-9]+) chromosome (?P<chrm>[IVXM]+)', line)
    strain = m.group('strain').lower()
    chrm = m.group('chrm')

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
        # subtract one to index from zero TODO is this correct? end is
        # inclusive
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

    return genes, strain, chrm, line
