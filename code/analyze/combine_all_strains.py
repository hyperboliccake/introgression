import re
import sys
import os
import math
import Bio.SeqIO
import copy
import gene_predictions
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../sim/')
import sim_analyze_hmm_bw as sim
sys.path.insert(0, '../misc/')
import seq_functions
import read_table
import read_fasta
import write_fasta
import mystats
import overlap

def get_range_seq(start, end, seq_fn):

    chrm_seq = read_fasta.read_fasta(seq_fn)[1][0]
    range_seq = chrm_seq[start:end+1]
    return range_seq

def get_ref_gene_seq(gene, gene_coords_fn, seq_fn):

    d1, labels = read_table.read_table_rows(gene_coords_fn, '\t', \
                                            header=False, key_ind=0)
    d = {}
    for g in d1:
        if d1[g][0] == '""':
            d[g] = d1[g][1:]
        else:
            d[d1[g][0]] = d1[g][1:]        

    gene_start = int(d[gene][2]) - 1
    gene_end = int(d[gene][3]) - 1
    chrm_seq = read_fasta.read_fasta(seq_fn)[1][0]
    gene_seq = chrm_seq[gene_start:gene_end+1]
    strand = d[gene][1]
    if strand == '-1':
        gene_seq = seq_functions.reverse_complement(gene_seq)
    assert gene_seq.startswith('atg') or gene_seq.startswith('ATG')
    assert gene_start < gene_end
    return gene_seq, gene_start, gene_end, strand

def get_inds_from_alignment(fn, flip_ref, rind=0, sind=1):
    headers, seqs = read_fasta.read_fasta(fn)
    n = len(seqs[0])
    ri = -1
    si = -1
    pr = []
    ps = []
    if flip_ref:
        rind = 1
        sind = 0
    for i in range(n):
        if seqs[sind][i] != gp.gap_symbol:
            si += 1
        if seqs[rind][i] != gp.gap_symbol:
            ri += 1
            pr.append(str(ri))
            ps.append(str(si))
    if flip_ref:
        return {'ps_ref':ps, 'ps_strain':pr}
    return {'ps_ref':pr, 'ps_strain':ps}


# by taking part of sequence aligned with reference coordinates
def get_range_seqs(strains, chrm, start, end, tag, gp_dir = '../'):
    # TODO this shouldn't actually be dependent on tag

    strain_range_seqs = {}
    for strain, d in strains:
        print strain
        fn = d + strain + '_chr' + chrm + gp.fasta_suffix
        chrm_seq = read_fasta.read_fasta(fn)[1][0]

        t = None
        try:
            t, labels = read_table.read_table_columns(gp.analysis_out_dir_absolute + \
                                                      tag + '/' + \
                                                      'site_summaries/predictions_' + \
                                                      strain + \
                                                      '_chr' + chrm + \
                                                      '_site_summary.txt.gz', '\t')
        except:
            # for par reference which doesn't have site summary file
            align_fn = gp_dir + gp.alignments_dir + \
                       '_'.join(gp.alignment_ref_order) + '_chr' + chrm + \
                       '_mafft' + gp.alignment_suffix
            t = get_inds_from_alignment(align_fn, True)

        ref_ind_to_strain_ind = dict(zip(t['ps_ref'], t['ps_strain']))

        start_strain = int(math.ceil(float(ref_ind_to_strain_ind[str(start)])))
        end_strain = int(math.floor(float(ref_ind_to_strain_ind[str(end)])))

        
        strain_range_seqs[strain] = (chrm_seq[start_strain:end_strain+1], \
                                      start_strain, end_strain)
    return strain_range_seqs


def choose_best_hit_rev(hits, query_fn, ref_chrm_fn, orf_headers, orf_seqs, start, end):
    # choosing best hit by reciprocal blast -> not reliable tho
    if len(hits) == 1:
        return hits[0][0]
    outfmt = '"6 sseqid slen evalue bitscore sstart send"'
    greatest_overlap = -1
    best_hit = hits[0]
    out_fn = str(start) + str(end) + '.out'
    orf_query_fn = str(start) + str(end) + '.txt'
    for hit in hits:
        header = None
        seq = None
        for i in range(len(orf_headers)):
            header = orf_headers[i]
            if header[1:].startswith(hit[0]):
                seq = orf_seqs[i]
                break
        f = open(orf_query_fn, 'w')
        f.write(seq + '\n')
        f.close()
        cmd_string = gp.blast_install_path + 'blastn' + \
                     ' -db ' + ref_chrm_fn + \
                     ' -query ' + orf_query_fn + \
                     ' -out ' + out_fn + \
                     ' -outfmt ' + outfmt
        os.system(cmd_string)
        f = open(out_fn, 'r')
        nhits = [line[:-1].split('\t') for line in f.readlines()]
        f.close()
        nstart = int(nhits[0][-2])
        nend = int(nhits[0][-1])
        # this division is hacky and unprincipled
        o = overlap.overlap(start, end, nstart, nend) / float(hit[1]) 
        if o > greatest_overlap:
            greatest_overlap = o
            best_hit = hit[0]
    os.remove(orf_query_fn)
    os.remove(out_fn)
    return best_hit

def choose_best_hit(hits, start, end, tag, strain, chrm, headers, seqs,\
                    strain_ind_to_ref_ind, gp_dir='../'):

    greatest_overlap = 0 # don't want to take overlaps of 0
    best_hit = None
    x = None
    seq = None
    orf_start = None
    orf_end = None
    strand = None
    x_max = None
    seq_max = None
    orf_start_max = None
    orf_end_max = None
    strand_max = None

    for hit in hits:
        orf_start = -1
        orf_end = -1
        for i in range(len(headers)):
            header = headers[i]
            if header[1:].startswith(hit[0]):
                chunk2 = header.split()[1]
                x = chunk2[:chunk2.find('_')]
                c1 = chunk2.find(':')
                c2 = chunk2.find(':', c1+1)
                seq = seqs[i]
                orf_start = int(chunk2[c1+1:c2])
                orf_end = int(chunk2[c2+1:])                
                strand = '1'
                if orf_start > orf_end:
                    temp = orf_end
                    orf_end = orf_start
                    orf_start = temp
                    strand = '-1'
                break
        current_start = int(math.ceil(float(strain_ind_to_ref_ind[str(orf_start)])))
        current_end = int(math.floor(float(strain_ind_to_ref_ind[str(orf_end)])))
        o = overlap.overlap(start, end, current_start, current_end)
        if o > greatest_overlap:
            greatest_overlap = o
            best_hit = hit[0]
            x_max = x
            seq_max = seq
            orf_start_max = orf_start
            orf_end_max = orf_end
            strand_max = strand
            seq_max = seq # don't need to reverse complement (blast does this)

    print greatest_overlap
    return best_hit, x_max, seq_max, orf_start_max, orf_end_max, strand_max

# by blasting ORFs
def get_gene_seqs(query_fn, strains, chrm, ref_chrm_fn, start, end, strand, tag,
                  strain_ind_to_ref_ind):
    
    #outfmt = '"6 qseqid sseqid slen qstart qend length mismatch gapopen gaps sseq"'
    outfmt = '"6 sseqid slen evalue bitscore"'

    strain_gene_seqs = {}
    out_fn = 'blast_chr' + chrm + '.out'
    for strain, d in strains:
        if strain != 'yjm1332':
            continue

        print '-', strain
        sys.stdout.flush()
        fn = d + 'orfs/' + strain + '_chr' + chrm + '_orfs' +  gp.fasta_suffix
        cmd_string = gp.blast_install_path + 'blastn' + \
                     ' -db ' + fn + \
                     ' -query ' + query_fn + \
                     ' -out ' + out_fn + \
                     ' -outfmt ' + outfmt
        #print cmd_string
        os.system(cmd_string)
        hits = [line[:-1].split('\t') for line in open(out_fn, 'r').readlines()]
        num_hits = len(hits)
        if len(hits) == 0:
            strain_gene_seqs[strain] = ('nohit', '', -1, -1, '')
            continue
        #best_orf_id = hits[0][0]
        headers, seqs = read_fasta.read_fasta(fn)
        best_orf_id, x, seq, orf_start, orf_end, orf_strand = \
            choose_best_hit(hits, start, end, tag, strain, chrm, headers, seqs, \
                            strain_ind_to_ref_ind[strain])
        print hits
        print best_orf_id
        print orf_strand, strand
        sys.exit()

        if best_orf_id == None or orf_strand != strand:
            strain_gene_seqs[strain] = ('nohit', '', -1, -1, '')
            continue
        strain_gene_seqs[strain] = (x, seq, orf_start, orf_end, orf_strand)
    os.remove(out_fn)
    return strain_gene_seqs
        

# can't actually count on annotations
def get_gene_seqs_gb(fn, gene, chrm):
    # get gene sequence for each strain
    gb_records = Bio.SeqIO.parse(fn, 'genbank')
    strain_gene_seqs = {}
    strains = set([])
    for strain_chrm_record in gb_records:
        desc = strain_chrm_record.description
        m = re.search(' (?P<strain>[a-zA-Z0-9]+) chromosome (?P<chrm>[IVXM]+)', \
                      desc)
        chrm_current = m.group('chrm')
        strain = m.group('strain').lower()
        strains.add(strain)
        #if len(strain_gene_seqs) > 82:
        #    break
        print strain, chrm_current
        if chrm_current != chrm:
            continue
        for feature in strain_chrm_record.features:
            if feature.type == 'CDS' and feature.qualifiers.has_key('gene') and \
               feature.qualifiers['gene'][0] == gene:
                desc = strain_chrm_record.description
                m = re.search(\
                        ' (?P<strain>[a-zA-Z0-9]+) chromosome (?P<chrm>[IVXM]+)', \
                        desc)
                seq = str(feature.extract(strain_chrm_record.seq).lower())
                start = str(feature.location.start)
                end = str(feature.location.end)
                strand = str(feature.location.strand)
                locus_tag = feature.qualifiers['locus_tag'][0]
                strain_gene_seqs[strain] = {'seq':seq, \
                                            'chrm':chrm, \
                                            'start':start, \
                                            'end':end, \
                                            'strand':strand,\
                                            'locus_tag':locus_tag}
            
                print '- found gene in', strain
    return strain_gene_seqs, list(strains)

# because don't have gb file for paradoxus...
def get_gene_seqs_fsa(fn, gene, chrm):
    f = open(fn, 'r')
    line = f.readline()
    m = 'SCER:' + gene
    while line != '':
        if m in line:
            seq = ''
            line = f.readline()
            while not line.startswith('>'):
                seq += line[:-1]
                line = f.readline()
            f.close()

            seqfa = open(gp.ref_dir['CBS432'] + 'CBS432_chr' + chrm + '.fa', 'r').read()
            seqfa = seqfa.replace('\n', '')
            if seq in seqfa:
                print 'found paradoxus seq'
            else:
                print 'did not find paradoxus seq'
            fg = open('a.txt', 'w')
            fg.write(seq + '\n')
            fg.write(seqfa + '\n')
            fg.close()
            return seq.lower()

        line = f.readline()
        
