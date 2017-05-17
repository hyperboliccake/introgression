import re
import sys
import os
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

    # pull out just region of interest from sequence
    #region_relative_start = index_ignoring_gaps(seq, start, block_start)
    #region_relative_end = index_ignoring_gaps(seq, end, block_start)
    #seq_region = seq[region_relative_start:region_relative_end+1]
    seq_region = seq[start:end+1]
    seq_region_start_ind = block_start + start

    a = [''] * len(seq_region)
    len_a = len(a)

    # just annotate whole alignment block and then pull out portion at
    # end; TODO make this more efficient
    for gene_name in genes:
        gene_start, gene_end = genes[gene_name]

        gene_relative_start = index_ignoring_gaps(seq_region, gene_start, \
                                                      seq_region_start_ind)
        gene_relative_end = index_ignoring_gaps(seq_region, gene_end, \
                                                    seq_region_start_ind)

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

    return a

def write_nonannotated(fn, refs, strain, block, i_start, i_end):

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

def write_annotated(block, refs, master_ref, strain, \
                        ref_to, ref_from, genes, \
                        ref_to_code, \
                        i_start_with_context, i_end_with_context, \
                        i_start, i_end, \
                        fn_annotated):

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
    # corresponding sequence (just so we know where gaps are)
    annotation_seq = block['strains'][master_ref]['sequence'][i_start_with_context:i_end_with_context+1]
    # loop through all positions and keep track of all unique genes
    # seen within the introgressed region, and also separately within
    # the region including surrounding context
    for i in range(len(annotation)):
        g = annotation[i]
        if g != '':
            if g not in gene_set:
                gene_set.append(g)
            # within introgressed region, rather than surrounding context
            if i >= start_offset and i <= end_offset:
                if g not in gene_set_region:
                    gene_set_region.append(g)
                    gene_lengths[g] = genes[g][1] - genes[g][0] + 1
                    gene_introgressed_bases_count[g] = 0
                # keep track of fraction of each gene that is
                # introgressed; this isn't perfect, but numerator is
                # number of (non-gap) sites within gene in
                # introgressed block (relative to master reference)
                # and denominator is gene length
                if annotation_seq[i] != gp.gap_symbol:
                    gene_introgressed_bases_count[g] += 1
        
    gene_set = filter(lambda x: x != '', gene_set)

    # now write to file
    f = open(fn_annotated, 'w')
    # first list the genes that are contained
    f.write('# ' + ' '.join(gene_set) + '\n\n')
    # start with references
    region_start_symbol = '<'
    region_end_symbol = '>'

    row_width = 80
    x = 0
    n = i_end_with_context - i_start_with_context + 1
    all_strains = refs + [strain]

    # first write order of strains
    for current_strain in all_strains:
        f.write(current_strain + '\n')
    f.write('match\n')
    f.write('introgressed\n\n')
    
    # get seqs
    seqs = {}
    for current_strain in all_strains:
        seqs[current_strain] = block['strains'][current_strain]['sequence']\
            [i_start_with_context:i_end_with_context+1]

    match_both = ' '
    match_neither = 'x'

    introgressed_code = ref_to_code[ref_from]

    introgressed = ' ' * (i_start - i_start_with_context) + \
        introgressed_code * (i_end - i_start + 1) + \
        ' ' * (i_end_with_context - i_end)

    # number of sites that match only species from and not species to
    # within the introgressed region (indication of how well-supported
    # the region is)
    ref_from_count = 0

    while x <= n:
        m = min(x + row_width, n)
        for strain in all_strains:
            seq = seqs[strain]
            for i in range(x, m):
                if annotation[i] != '':
                    f.write(seq[i].upper())
                else:
                    f.write(seq[i].lower())
            f.write('\n')

        # row for which reference strain matches (only accounts for 2
        # references - the reference for the species the strain is in
        # and the reference that the region is predicted to be
        # introgressed from)
        for i in range(x, m):
            base_x = seqs[strain][i]
            base_from = seqs[ref_from][i]
            base_to = seqs[ref_to][i]
            if base_x == gp.gap_symbol or \
                    base_from == gp.gap_symbol or \
                    base_to == gp.gap_symbol:
                f.write(' ')
            elif base_x == base_to:
                if base_x == base_from:
                    # matches both refs
                    f.write(match_both)
                else:
                    # matches ref 0 only 
                    f.write(ref_to_code[ref_to])
            else:
                if base_x == base_from:
                    # matches ref 1 only
                    f.write(ref_to_code[ref_from])
                    # increment number of matches to reference from,
                    # but only if we're in the introgressed region and
                    # not just surrounding context
                    if introgressed[i] == introgressed_code:
                        ref_from_count += 1
                else:
                    # matches neither ref
                    f.write(match_neither)
        f.write('\n')

        # row for which sites called introgressed
        for i in range(x, m):
            f.write(introgressed[i])
        f.write('\n\n')        

        # next set of rows
        x += row_width

    return gene_set_region, gene_introgressed_bases_count, gene_lengths, ref_from_count

def write_region_alignment(block, entry, genes, \
                               strain, master_ref, refs, \
                               ref_to_code, \
                               ref_to, ref_from, \
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
        

    # write alignment just for region without context
    write_nonannotated(fn, refs, strain, block, i_start, i_end)

    # write annotated alignment
    gene_set_region, gene_introgressed_bases_count, gene_lengths, ref_from_count = \
        write_annotated(block, refs, master_ref, strain, \
                            ref_to, ref_from, genes, \
                            ref_to_code, \
                            i_start_with_context, i_end_with_context, \
                            i_start, i_end, \
                            fn_annotated)
        

    # store the genes we found in this region (NOT including stuff
    # only in context)
    entry['genes'] = gene_set_region

    # introgressed fractions of genes
    x = {}
    for gene in gene_lengths:
        x[gene] = float(gene_introgressed_bases_count[gene]) / gene_lengths[gene]
    entry['genes_introgressed_fractions'] = x

    entry['ref_from_count'] = ref_from_count

def read_gene_file(fn):
    f = open(fn, 'r')
    genes = {}
    line = f.readline()
    while line != '':
        line = line.split('\t')
        genes[line[0]] = (int(line[1]), int(line[2]))
        line = f.readline()
    f.close()
    return genes

def write_gene_file(genes, fn):
    f = open(fn, 'w')
    for gene in genes:
        start, end = genes[gene]
        f.write(gene + '\t' + str(start) + '\t' + str(end) + '\n')
    f.close()

def read_genes(fn, fn_genes):

    if os.path.isfile(fn_genes):
        return read_gene_file(fn_genes)

    f = open(fn, 'r')
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
    f.close()
    write_gene_file(genes, fn_genes)

    return genes

def summarize_gene_info(fn_all, introgressed_genes, tag, threshold=0):
    
    f_all = open(fn_all, 'w')
    f_all.write('gene\tnumber_strains\taverage_introgressed_fraction\taverage_number_non_gap\taverage_ref_from_count\n')

    f_gene_heading = 'region_id\tstrain\tintrogressed_fraction\tnumber_non_gap\tref_from_count\n'

    for gene in introgressed_genes:
        # keyed by strain, because gene can be broken across multiple
        # alignment blocks/regions for the same strain
        sum_introgressed_fraction = {}
        sum_number_non_gap = {}
        sum_ref_from_count = {}
        fn_gene = gp.analysis_out_dir_absolute + tag + '/genes/' + gene + '.txt'
        if not os.path.exists(os.path.dirname(fn_gene)):
            os.makedirs(os.path.dirname(fn_gene))
        f_gene = open(fn_gene, 'w')
        f_gene.write(f_gene_heading)
        for entry in introgressed_genes[gene]:
            region_id, strain, introgressed_fraction, number_non_gap, ref_from_count = entry
            if ref_from_count >= threshold:
                if strain not in sum_introgressed_fraction:
                    sum_introgressed_fraction[strain] = 0
                    sum_number_non_gap[strain] = 0
                    sum_ref_from_count[strain] = 0
                sum_introgressed_fraction[strain] += introgressed_fraction
                sum_number_non_gap[strain] += number_non_gap
                sum_ref_from_count[strain] += ref_from_count
            f_gene.write(region_id + '\t' + strain + '\t' + \
                             str(introgressed_fraction) + '\t' + str(number_non_gap) + '\t' + \
                             str(ref_from_count) + '\n')

        f_gene.close()

        # now do averaging over strains
        num_strains = len(sum_introgressed_fraction)
        if num_strains == 0:
            continue

        avg_introgressed_fraction = sum(sum_introgressed_fraction.values()) / \
            float(num_strains)
        avg_number_non_gap = sum(sum_number_non_gap.values()) / \
            float(num_strains)
        avg_ref_from_count = sum(sum_ref_from_count.values()) / \
            float(num_strains)

        f_all.write(gene + '\t' + str(num_strains) + '\t' + \
                        str(avg_introgressed_fraction) + '\t' + \
                        str(avg_number_non_gap) + '\t' + \
                        str(avg_ref_from_count) + '\n')

    f_all.close()
