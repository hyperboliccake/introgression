import re
import sys
import os
import copy
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../sim/')
import sim_analyze_hmm_bw as sim
sys.path.insert(0, '../misc/')
import write_fasta

def index_ignoring_gaps(s, i, s_start):
    '''returns the index of the ith (starting at 0) non-gap character in
    s, given that the start index of s is s_start (instead of needing
    to be 0); for example, s='-A-AA-A', i=2, s_start=0 => returns 4
    '''

    assert s_start >= 0

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

def get_ref_match_by_site(seqs, labels):

    # for master: matches _only_ that ref
    # for other refs: matches that ref but not master

    nrefs = len(seqs) - 1
    nsites = len(seqs[0])
    ref_match_by_site = [[' ']*nsites for r in range(nrefs)]

    for i in range(nsites):
        gap = False
        for seq in seqs:
            if seq[i] == gp.gap_symbol:
                gap = True
                break
        if gap:
            continue

        if seqs[0][i] == seqs[-1][i]:
            ref_match_by_site[0][i] = labels[0][0]
            
        for r in range(1, nrefs):
            if seqs[r][i] == seqs[-1][i]:
                # matches this ref and master ref -> both blank
                if seqs[0][i] == seqs[-1][i]:
                    ref_match_by_site[0][i] = ' '
                # matches this ref and not master ref -> label
                else:
                    ref_match_by_site[r][i] = labels[r][0]

            else:
                # matches master ref but not this ref -> already good
                if seqs[0][i] == seqs[-1][i]:
                    pass
                # matches neither ref -> . for both
                else:
                    ref_match_by_site[r][i] = '.'
                    ref_match_by_site[0][i] = '.'
                

    return [''.join(s) for s in ref_match_by_site]
    

def get_ref_match_by_site_2(seqs, labels):

    # this is for getting a string for matching/mismatching each
    # reference individually

    nrefs = len(seqs) - 1
    nsites = len(seqs[0])
    ref_match_by_site = [[] for r in range(nrefs)]

    for i in range(nsites):
        add = []
        gap = (seqs[-1][i] == gp.gap_symbol)
        for r in range(nrefs):
            if gap or seqs[r][i] == gp.gap_symbol:
                add = ['.'] * nrefs
                break
            elif seqs[r][i] == seqs[-1][i]:
                add.append(labels[r][0])
            else:
                add.append(' ')
        for r in range(nrefs):
            ref_match_by_site[r].append(add[r])

    return [''.join(s) for s in ref_match_by_site]

def get_genes_by_site(genes, seq):

    genes_by_site = [None for site in seq]
    for gene_name in genes:
        gene = genes[gene_name]
        gene_start_ind = index_ignoring_gaps(seq, gene[0], 0)
        gene_end_ind = index_ignoring_gaps(seq, gene[1], 0)
        for i in range(gene_start_ind, gene_end_ind+1):
            genes_by_site[i] = gene_name
    return genes_by_site

def get_introgressed_by_site(regions, seq):

    introgressed_by_site = [' ' for site in seq]
    for entry in regions:
        start_ind = index_ignoring_gaps(seq, entry[0], 0)
        end_ind = index_ignoring_gaps(seq, entry[1], 0)
        for i in range(start_ind, end_ind+1):
            introgressed_by_site[i] = 'i'
    return ''.join(introgressed_by_site)
    

def write_region_alignment(headers, seqs, fn, start, end, master_ind):
    
    relative_start = max(0, index_ignoring_gaps(seqs[master_ind], start, 0))
    relative_end = index_ignoring_gaps(seqs[master_ind], end, 0)
    
    region_seqs = [seq[relative_start:relative_end+1] for seq in seqs]

    write_fasta.write_fasta(headers, region_seqs, fn)

def get_genes_in_region(start, end, genes):
    
    region_genes = []
    for gene_name in genes:
        gene_start, gene_end = genes[gene_name]
        if (gene_start > start and gene_start <= end) or \
           (gene_end > start and gene_end <= end):
            region_genes.append((gene_name, gene_start, gene_end))

    region_genes.sort(key=lambda x: x[1])
    return region_genes

def write_region_alignment_annotated(labels, seqs, fn, start, end, \
                                     master_ind, genes, ref_match_by_site, \
                                     genes_by_site, \
                                     introgressed_by_site, context):

    relative_start_with_context = \
        max(0, index_ignoring_gaps(seqs[master_ind], start-context, 0))
    relative_start = max(0, index_ignoring_gaps(seqs[master_ind], start, 0))
    relative_end = index_ignoring_gaps(seqs[master_ind], end, 0)
    relative_end_with_context = index_ignoring_gaps(seqs[master_ind], end+context, 0)
    
    region_seqs = [seq[relative_start_with_context:relative_end_with_context+1] \
                   for seq in seqs]

    # for reference matching lines
    ref_match_strings = []
    for r in ref_match_by_site:
        ref_match_strings.append(\
            r[relative_start_with_context:relative_end_with_context+1])

    # for gene line
    region_genes = \
        genes_by_site[relative_start_with_context:relative_end_with_context+1]
    region_genes_set = list(set(region_genes))
    try:
        region_genes_set.remove(None)
    except:
        pass
    region_genes_set.sort(key=lambda x: genes[x][1])
    gene_string = ''.join([' ' if entry == None else '=' for entry in region_genes])

    # for introgression line
    introgressed_string = \
        introgressed_by_site[relative_start_with_context:relative_start] + \
        (relative_end + 1 - relative_start) * 'I' + \
        introgressed_by_site[relative_end+1:relative_end_with_context+1]
    assert len(introgressed_string) == len(region_seqs[0]), \
        str(len(introgressed_string)) + ' ' + str(len(region_seqs[0]))

    # line labels
    f = open(fn, 'w')
    for label in labels:
        f.write(label + '\n')
    # assume master ref comes first
    f.write('matches only ' + labels[0] + '\n')
    # and assume ref seqs come before predict seq
    for label in labels[1:-1]: 
        f.write('matches ' + label + ' and mismatches ' + labels[0] + '\n')
    f.write('genes: ' + ' '.join(region_genes_set) + '\n')
    f.write('introgressed\n\n')

    # ...finally write the lines, width characters at a time
    width = 80
    i = 0
    while i < len(region_seqs[0]):
        # sequences
        for s in range(len(region_seqs)):
            f.write(region_seqs[s][i:i+width] + '\n')
        # matches references?
        for s in range(len(ref_match_strings)):
            f.write(ref_match_strings[s][i:i+width] + '\n')
        # genes
        f.write(gene_string[i:i+width] + '\n')
        # introgressed
        f.write(introgressed_string[i:i+width] + '\n\n')

        i += width

    return relative_start, relative_end

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

def write_region_summary_header(refs, f):
    f.write('region_id\tstrain\tchromosome\tpredicted_species\tstart\tend\t' + \
            'number_non_gap\t')
    f.write('\t'.join(['number_match_' + ref for ref in refs]) + '\t')
    f.write('\t'.join(['number_match_only_' + ref for ref in refs]) + '\t')
    f.write('number_mismatch_all_refs\n')

def write_region_summary_line(region, strain, chrm, predicted_species, seqs, labels, 
                              start, end, f):

    # region_id [strain chromosome predicted_species start end number_non_gap]
    # number_match_ref1 number_match_ref2 number_match_only_ref1
    # number_match_ref2_not_ref1 number_mismatch_all_ref

    sep = '\t'

    f.write(region[3] + sep + strain + sep + chrm + sep + predicted_species + \
            sep + str(region[0]) + sep + str(region[1]) + sep + \
            str(region[2]) + sep)

    ids = [0] * (len(seqs) - 1)
    # matches only master ref, or ref but not master ref
    unique_ids = [0] * (len(seqs) - 1)
    mismatch_all = 0
    for i in range(start, end+1):
        gap = False
        for s in seqs:
            if s[i] == gp.gap_symbol:
                gap = True
                break
        if gap:
            continue

        match_refs = []
        for r in range(0, len(seqs)-1):
            if seqs[r][i] == seqs[-1][i]:
                match_refs.append(True)
                ids[r] += 1
            else:
                match_refs.append(False)
        if match_refs[0]:
            # false = 0, true = 1
            if sum(match_refs[1:]) == 0:
                unique_ids[0] += 1
        else:
            if sum(match_refs[1:]) == 0:
                mismatch_all += 1
                continue
            for r in range(1, len(seqs) - 1):
                unique_ids[r] += match_refs[r]
                
    f.write(sep.join([str(x) for x in ids]) + sep)
    f.write(sep.join([str(x) for x in unique_ids]) + sep)
    f.write(str(mismatch_all) + '\n')
    f.flush()

def write_genes_for_each_region_summary_line(region_id, genes_by_site, gene_summary, \
                                             start, end, seq, f):
    
    # region_id num_genes gene frac_intd gene frac_intd
    genes = genes_by_site[start:end+1]
    genes_set = list(set(genes))
    try:
        genes_set.remove(None)
    except:
        pass
    seq_region = seq[start:end+1]
    gene_site_counts = dict(zip(genes_set, [0]*len(genes_set)))
    for i in range(len(seq_region)):
        if seq_region[i] != gp.gap_symbol and genes[i] != None:
            gene_site_counts[genes[i]] += 1
    frac_intd = {}
    for gene in genes_set:
        gene_length = gene_summary[gene][1] - gene_summary[gene][0] + 1
        frac_intd[gene] = float(gene_site_counts[gene]) / gene_length
    
    sep = '\t'
    f.write(region_id + sep)
    f.write(str(len(genes_set)))
    for gene in genes_set:
        f.write(sep + gene + sep + str(frac_intd[gene]))
    f.write('\n')
    f.flush()

    return frac_intd

def write_regions_for_each_strain(regions, f):

    # strain num_regions region length region length
    sep = '\t'
    for strain in regions:
        f.write(strain + sep)
        num_regions = sum([len(regions[strain][chrm]) for chrm in gp.chrms])
        f.write(str(num_regions))
        for chrm in gp.chrms:
            for region in regions[strain][chrm]:
                region_length = region[1] - region[0] + 1
                f.write(sep + region[3] + sep + str(region_length))
        f.write('\n')
    f.flush()

def write_genes_for_each_strain(strain_genes_dic, f):

    # strain num_genes gene frac_intd gene frac_intd
    sep = '\t'
    for strain in strain_genes_dic:
        f.write(strain + sep + str(len(strain_genes_dic[strain])))
        for gene in strain_genes_dic[strain]:
            f.write(sep + gene + sep + str(strain_genes_dic[strain][gene]))
        f.write('\n')
    f.flush()

def write_strains_for_each_gene_lines(gene_strains_dic, f):

    # (this is actually the same as above function, but it's confusing
    # to have the wrong names)

    # gene num_strains strain frac_intd strain frac_intd
    sep = '\t'
    for gene in gene_strains_dic:
        f.write(gene + sep + str(len(gene_strains_dic[gene])))
        for strain in gene_strains_dic[gene]:
            f.write(sep + strain + sep + str(gene_strains_dic[gene][strain]))
        f.write('\n')
    f.flush()

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

"""
def summarize_gene_info(fn_all, fn_strains, fn_strains_g, \
                            introgressed_genes, gene_info, tag, threshold=0):
    
    f_all = open(fn_all, 'w')
    f_all.write('gene\tchromosome\tstart\tend\tnumber_strains\taverage_introgressed_fraction\taverage_number_non_gap\taverage_ref_from_count\n')

    f_gene_heading = 'region_id\tstrain\tstart\tend\tintrogressed_fraction\tnumber_non_gap\tref_from_count\n'

    strain_genes = {}

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
            strain = entry['strain']
            if entry['ref_from_count'] >= threshold:
                if strain not in sum_introgressed_fraction:
                    sum_introgressed_fraction[strain] = 0
                    sum_number_non_gap[strain] = 0
                    sum_ref_from_count[strain] = 0
                sum_introgressed_fraction[strain] += entry['introgressed_fraction']
                sum_number_non_gap[strain] += entry['number_non_gap']
                sum_ref_from_count[strain] += entry['ref_from_count']
                if strain not in strain_genes:
                    strain_genes[strain] = []
                strain_genes[strain].append(gene)

            f_gene.write(entry['region_id'] + '\t' + \
                             strain + '\t' + \
                             str(entry['start']) + '\t' + \
                             str(entry['end']) + '\t' + \
                             str(entry['introgressed_fraction']) + '\t' + \
                             str(entry['number_non_gap']) + '\t' + \
                             str(entry['ref_from_count']) + '\n')

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

        f_all.write(gene + '\t' + \
                        gene_info[gene][0] + '\t' + \
                        str(gene_info[gene][1]) + '\t' + \
                        str(gene_info[gene][2]) + '\t' + \
                        str(num_strains) + '\t' + \
                        str(avg_introgressed_fraction) + '\t' + \
                        str(avg_number_non_gap) + '\t' + \
                        str(avg_ref_from_count) + '\n')

    f_all.close()

    f_strains = open(fn_strains, 'w')
    f_strains_g = open(fn_strains_g, 'w')
    for strain in strain_genes:
        line = strain + '\t' + str(len(strain_genes[strain]))
        f_strains.write(line + '\n')
        line += '\t' + '\t'.join(strain_genes[strain])
        f_strains_g.write(line + '\n')
    f_strains.close()
    f_strains_g.close()
"""
