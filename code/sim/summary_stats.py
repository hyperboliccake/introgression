# calculate sequence identity between all pairs of cerevisiae in 100
# genomes set, and also with paradoxus reference

# do this from the alignments? just the three-way portions?

# also seq id for simulations

import sys
import os
import re
sys.path.insert(0, '../misc/')
import read_fasta


def seq_id(a, b, l = -1, use_gaps = False):
    assert len(a) == len(b)
    ndiff = 0
    ntotal = 0
    for i in xrange(len(a)):
        # use sequence a as denominator
        if a[i] != '-':
            if b[i] != '-' or use_gaps:
                ntotal += 1
                if a[i] != b[i]:
                    ndiff += 1
    if l == -1:
        if ntotal == 0:
            return -1, 0
        return 1 - float(ndiff) / ntotal, ntotal
    return 1 - float(ndiff) / l

def mean(a):
    return float(sum(a))/len(a)


def maf_id(fn, ref1 = 'S288c', ref2 = 'CBS432'):
    headers, seqs = read_fasta.read_fasta(fn)
    id1, den1 = seq_id(seqs[2], seqs[0])
    id2, den2 = seq_id(seqs[2], seqs[1])
    return id1, id2, den1, den2

# for muscle output
def maf_id_old(fn, only_threeway, ref1 = 'S288c', ref2 = 'CBS432'):
    assert only_threeway, 'reading non-threeway parts of the alignment not yet implemented'
    f = open(fn, 'r')
    line = f.readline()
    id1 = []
    id2 = []
    len1 = []
    len2 = []
    while line != '':
        if line[0] == 'a':
            s1 = f.readline()
            s2 = f.readline()
            if (s2[0] == 's') and (ref1 in s1) and (ref2 in s2):
                s3 = f.readline()
                if s3[0] == 's':
                    s1 = s1.split()
                    s2 = s2.split()
                    s3 = s3.split()
                    r1 = seq_id(s3[-1], s1[-1])
                    r2 = seq_id(s3[-1], s2[-1])
                    # denominator should be the same in both cases
                    # but only if using columns with gaps in reference
                    # assert r1[1] == r2[1]
                    id1.append(r1[0])
                    id2.append(r2[0])
                    len1.append(r1[1])
                    len2.append(r2[1])
        line = f.readline()
    f.close()
    num1 = sum([id1[i] * len1[i] for i in xrange(len(id1))])
    num2 = sum([id2[i] * len2[i] for i in xrange(len(id2))])
    den1 = float(sum(len1))
    den2 = float(sum(len2))
    return num1 / den1, num2 / den2, den1, den2
    


if sys.argv[1] == '100':

    # strategy here is to use threeway alignments to calculate
    # sequence identity, i.e. assume genomes are mostly aligned so
    # this is valid

    chrms = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XIV']

    # get all strain names
    alignment_dir = '../../alignments/genbank/'
    prefix = 'S288c_CBS432_'
    fns = os.listdir(alignment_dir)
    fns = filter(lambda fn: fn.endswith('mafft.maf'), fns)
    strains = set()
    for fn in fns:
        try:
            m = re.match(prefix + '(?P<strain>[a-zA-Z0-9]+)_chr', fn)
            strains.add(m.group('strain'))
        except:
            pass
    # process alignments for each strain and chromosome
    id_cer = []
    id_par = []
    len_cer = []
    len_par = []
    for strain in strains:
        print strain
        id_strain_cer = []
        id_strain_par = []
        len_strain_cer = []
        len_strain_par = []
        for chrm in chrms:
            id_chrm_cer, id_chrm_par, len_chrm_cer, len_chrm_par = maf_id(alignment_dir + prefix + strain + '_chr' + chrm + '_mafft.maf')
            id_strain_cer.append(id_chrm_cer)
            id_strain_par.append(id_chrm_par)
            len_strain_cer.append(len_chrm_cer)
            len_strain_par.append(len_chrm_par)
        num_strain_cer = sum([id_strain_cer[i] * len_strain_cer[i] for i in xrange(len(chrms))])
        num_strain_par = sum([id_strain_par[i] * len_strain_par[i] for i in xrange(len(chrms))])
        den_strain_cer = float(sum(len_strain_cer))
        den_strain_par = float(sum(len_strain_par))
        id_cer.append(num_strain_cer/den_strain_cer)
        id_par.append(num_strain_par/den_strain_par)
        len_cer.append(den_strain_cer)
        len_par.append(den_strain_par)
        print id_cer, id_par
    num_cer = sum([id_cer[i] * len_cer[i] for i in xrange(len(id_cer))])
    num_par = sum([id_par[i] * len_par[i] for i in xrange(len(id_par))])
    den_cer = float(sum(len_cer))
    den_par = float(sum(len_par))
    print 'average id between cer and cer ref:', num_cer/den_cer
    print '(all strains:', id_cer, len_cer, ')'
    print 'average id between cer and par ref:', num_par/den_par
    print '(all strains:', id_par, len_par, ')'    

elif sys.argv[1] == 'sim':
    # for simulations, specifically with one par and one cer population
    out_dir = '../../results/sim/'
    out_prefix = 'sim_out_'
    out_param_fn = '../sim/sim_multi_model_args.txt'
    params = [line.split() for line in open(out_param_fn, 'r').readlines()]
    for p in params:
        num_par = int(p[3])
        num_cer = int(p[4])
        par_ref_ind = 0
        cer_ref_ind = num_par
        seq_len = int(p[7])            
        d = {'par_ref':{'par':[], 'cer':[]}, 'cer_ref':{'par':[], 'cer':[]}}

        f = open(out_dir + out_prefix + p[0] + '.txt', 'r')
        line = f.readline()
        while line != '':
            m = re.match('segsites: (?P<nseg>[0-9]+)', line) 
            if m != None:
                # positions:
                f.readline()
                # begin sequences
                line = f.readline()
                seqs = []
                while line != '' and line != '\n':
                    seqs.append(line[:-1])
                    line = f.readline()
                d_temp = {'par_ref':{'par':[], 'cer':[]}, 'cer_ref':{'par':[], 'cer':[]}}
                for i in range(0, num_par):
                    if i != par_ref_ind:
                        d_temp['par_ref']['par'].append(seq_id(seqs[par_ref_ind], seqs[i], l=seq_len))
                    d_temp['cer_ref']['par'].append(seq_id(seqs[cer_ref_ind], seqs[i], l=seq_len))
                for i in range(num_par, num_par + num_cer):
                    if i != cer_ref_ind:
                        d_temp['cer_ref']['cer'].append(seq_id(seqs[cer_ref_ind], seqs[i], l=seq_len))
                    d_temp['par_ref']['cer'].append(seq_id(seqs[par_ref_ind], seqs[i], l=seq_len))
                d['par_ref']['par'].append(mean(d_temp['par_ref']['par']))
                d['par_ref']['cer'].append(mean(d_temp['par_ref']['cer']))
                d['cer_ref']['par'].append(mean(d_temp['cer_ref']['par']))
                d['cer_ref']['cer'].append(mean(d_temp['cer_ref']['cer']))
            line = f.readline()
        f.close()

        print '-----'
        print out_prefix + p[0]
        print p
        print 'average id within par:', mean(d['par_ref']['par'])
        print 'average id within cer:', mean(d['cer_ref']['cer'])
        print 'average id between cer ref and par:', mean(d['cer_ref']['par'])
        print 'average id between par ref and cer:', mean(d['par_ref']['cer'])

