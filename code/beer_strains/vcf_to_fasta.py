import sys

def read_vcf(fn):
    f = open(fn, 'r')

    v = []
    line = f.readline()
    while line != '':
        chrm, ps, x1, ref, alt, x2, x3, x4, x5, x6 = line.split('\t')
        if chrm not in v:
            v[chrm] = {}
        v[chrm][ps] = (ref, alt)
        line = f.readline()
    f.close()
    return v

def vcf_to_fasta(v, fn_ref, fn_out):
    f_ref = open(fn_ref, 'r')
    f_out = open(fn_out, 'w')

    line = f_ref.readline()
    while line != '':
        
        line = f_ref.readline()
    
v = read_vcf(sys.argv[1])
vcf_to_fasta(v, sys.argv[2])
