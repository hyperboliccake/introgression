import sys
import os
from convert_coordinates import *
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc/')
import read_fasta

gp_dir = '../'
fns = os.listdir(gp_dir + gp.alignments_dir)
fns = filter(lambda fn: fn.endswith(gp.alignment_suffix), fns)

for fn in fns:
    print fn

    x = fn.split('_')
    chrm = x[-2]
    strain_names = x[0:-2]
    headers, seqs = read_fasta.read_fasta(gp_dir + gp.alignments_dir + fn)
    
    # for each index in cer reference, get index in other strain
    # (either par reference for 2-way alignment or cer strain for
    # 3-way)
    coord_fn = gp.analysis_out_dir_absolute + 'coordinates/' + \
               strain_names[0] + '_to_' + strain_names[-1] + \
               '_' + chrm + '.txt.gz'
    write_coordinates(convert(seqs[0], seqs[-1]), coord_fn)

    # for each index in other strain, get index in cer reference
    coord_fn = gp.analysis_out_dir_absolute + 'coordinates/' + \
               strain_names[-1] + '_to_' + strain_names[0] + \
               '_' + chrm + '.txt.gz'
    write_coordinates(convert(seqs[-1], seqs[0]), coord_fn)
