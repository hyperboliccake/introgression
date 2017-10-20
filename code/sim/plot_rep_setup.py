import sys
import os
import sim_process
import sim_predict
sys.path.append('..')
import global_params as gp

def write_coding_table(seqs_coded, files, rep, header):
    
    for i in range(len(seqs_coded)):
        f = files[i]
        if header:
            f.write('site\tcode\trep\n')
        for j in range(len(seqs_coded[i])):
            f.write(str(j) + '\t' + seqs_coded[i][j] + '\t' + str(rep) + '\n')

def write_blocks_table(blocks_dic, files, rep, header):
    for ind in blocks_dic.keys():
        f = files[ind]
        if header:
            f.write('type\tspecies\tstart\tend\trep\n')
        for block_type in blocks_dic[ind].keys():
            for species in blocks_dic[ind][block_types].keys():
                for block in blocks_dic[ind][block_type][species]:
                    f.write(block_type + '\t' + species + '\t' + \
                            str(block[0]) + '\t' + str(block[1]) + '\t' + \
                            str(rep) + '\n')
            
