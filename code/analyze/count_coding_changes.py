import sys
import os
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc/')
import seq_functions
import read_fasta

def get_aligned_genes(fn, strains):
    headers, seqs = read_fasta.read_fasta(fn)
    d = {}
    for i in range(len(headers)):
        strain = headers[i][1:].split()[0]
        if strain in strains:
            d[strain] = seqs[i]
    n = len(d.values()[0])
    remove_columns = []
    for i in range(n):
        all_gap = True
        for strain in d.keys():
            if d[strain][i] != gp.gap_symbol:
                all_gap = False
                break
        if all_gap:
            remove_columns.append(i)
    for i in remove_columns[::-1]:
        for strain in d.keys():
            d[strain] = d[strain][:i] + d[strain][i+1:]
    return d


def ambiguous(gene, ref_start, ref_end, coords, orfs):
    # this function is for determining whether there's unambiguous
    # presence of gene in another strain, in this case by there being
    # an orf at the corresponding coordinates (with no gaps or
    # differences in length)

    start = coords[ref_start]
    end = coords[ref_end]
    if start == round(start) and end == round(end) and \
       (start, end) in orfs and end - start == ref_end - ref_start:
        return False
    return True


def count_coding(seq_master, seq_ref, seq_strain, start, end):
    
    if not seq_master.startswith('ATG'):
        seq_master = seq_functions.reverse_complement(seq_master)
        assert seq_master.startswith('ATG'), seq_master
    if not seq_ref.startswith('ATG'):
        seq_ref = seq_functions.reverse_complement(seq_ref)
        assert seq_ref.startswith('ATG')
    if not seq_strain.startswith('ATG'):
        seq_strain = seq_functions.reverse_complement(seq_strain)
        assert seq_strain.startswith('ATG')

    if not start % 3 == 0:
        if start % 3 == 1:
            start += 2
        elif start % 3 == 2:
            start += 1
    if not end % 3 == 2:
        if end % 3 == 0:
            end -= 1
        elif end % 3 == 1:
            end -= 2

    s = list(seq_master)
    for i in range(len(s)):
        if seq_strain[i] == seq_ref[i]:
            s[i] = seq_ref[i]
    s = ''.join(s)
    t = 0
    for i in range(start, end + 1, 3):
        if s[i:i+3] != seq_master[i:i+3]:
            t += 1
    t_syn = 0
    t_non = 0
    a = seq_functions.translate(s)
    a_master = seq_functions.translate(seq_master)
    for i in range(start/3, (end + 1)/3):
        if a[i] != a_master[i]:
            t_non += 1
    t_syn = t - t_non
    return t_syn, t_non


def count_coding_with_gaps(seq_master, seq_ref, seq_strain, start, end):

    print seq_master
    print seq_ref
    print seq_strain
    print start, end

    seq_master = seq_master.upper()
    seq_ref = seq_ref.upper()
    seq_strain = seq_strain.upper()
    
    ind_master = 0
    ind_ref = 0
    ind_strain = 0

    if not start % 3 == 0:
        if start % 3 == 1:
            start += 2
        elif start % 3 == 2:
            start += 1
    if not end % 3 == 2:
        if end % 3 == 0:
            end -= 1
        elif end % 3 == 1:
            end -= 2

    # count changed codons (relative to master) that are synonymous
    t_syn = 0
    # count changed codons (relative to master) that are nonsynonymous
    t_non = 0

    # count changed codons (relative to master) that match ref that
    # are synonymous
    t_syn_ref = 0
    # count changed codons (relative to master) that match ref that
    # are nonsynonymous
    t_non_ref = 0

    t_insert = 0
    t_delete = 0
    t_insert_ref = 0
    t_delete_ref = 0

    # only count codons within introgressed region
    in_region = False

    # things are more confusing if there's a frameshift
    frameshift = False
    frameshift_count = 0

    gene_delete = 0
    gene_delete_ref = 0

    if seq_strain.replace(gp.gap_symbol, '') == '':
        gene_delete = 1
        if seq_ref.replace(gp.gap_symbol, '') == '':
            gene_delete_ref = 1
        return t_syn, t_non, t_syn_ref, t_non_ref, \
            t_insert/3.0, t_delete/3.0, t_insert_ref/3.0, t_delete_ref/3.0, \
            gene_delete, gene_delete_ref, frameshift_count

    for i in range(0, len(seq_master), 3):
        if ind_master >= start and ind_master <= end:
            in_region = True
        elif in_region:
            break

        codon_master = seq_master[i:i+3]
        codon_ref = seq_ref[i:i+3]
        codon_strain = seq_strain[i:i+3]

        gaps_master = codon_master.count(gp.gap_symbol)
        gaps_ref = codon_ref.count(gp.gap_symbol)
        gaps_strain = codon_strain.count(gp.gap_symbol)

        ind_master += (3 - gaps_master)
        ind_ref += (3 - gaps_ref)
        ind_strain += (3 - gaps_strain)

        if not in_region:
            continue

        if ind_master % 3 == 0 and ind_ref % 3 == 0 and ind_strain % 3 == 0:
            frameshift = False
        else:
            frameshift = True
            frameshift_count = 1

        if codon_strain != codon_master:

            aa_master = seq_functions.codon_table.get(codon_master)
            aa_ref = seq_functions.codon_table.get(codon_ref)
            aa_strain = seq_functions.codon_table.get(codon_strain)

            if aa_master == None or aa_strain == None:
                if gaps_master > gaps_strain:
                    t_insert += gaps_master - gaps_strain
                else:
                    t_delete += gaps_strain - gaps_master
                if codon_strain == codon_ref:
                    if gaps_master > gaps_strain:
                        t_insert_ref += gaps_master - gaps_strain
                    else:
                        t_delete_ref += gaps_strain - gaps_master

            if frameshift:
                continue

            if gaps_master > 0 or gaps_strain > 0:
                continue

            if aa_strain == aa_master:
                t_syn += 1
            else:
                t_non += 1

            if gaps_ref > 0:
                continue

            if codon_strain == codon_ref:

                if aa_strain == aa_master:
                    t_syn_ref += 1
                else:
                    t_non_ref += 1

    print t_syn, t_non, t_syn_ref, t_non_ref
    print t_insert, t_delete, t_insert_ref, t_delete_ref
    print frameshift
    return t_syn, t_non, t_syn_ref, t_non_ref, \
        t_insert/3.0, t_delete/3.0, t_insert_ref/3.0, t_delete_ref/3.0, \
        gene_delete, gene_delete_ref, frameshift_count

