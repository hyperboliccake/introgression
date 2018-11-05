import sys
import os
sys.path.insert(0, '../misc/')
import read_fasta
import write_fasta

def concatenate_fastas(fns, fn_out, remove_gaps):
    strains = read_fasta.read_fasta(fns[0])[0]
    concat_seqs = dict(zip(strains, ['' for s in strains]))
    for fn in fns:
        headers, seqs = read_fasta.read_fasta(fn)
        for i in range(len(seqs)):
            concat_seqs[headers[i]] += seqs[i]
    f = open(fn_out, 'w')
    for strain in strains:
        f.write(strain + '\n')
        f.write(concat_seqs[strain] + '\n')
    f.close()
            
