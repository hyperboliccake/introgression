import sys
import os
import sim_process
import sim_predict
sys.path.append('..')
import global_params as gp

def write_combined_file(f, rep, codings, introgressed, \
                        actual_block_type, introgressed_ref, probs, header):
    #  rep site coding predicted actual actual_ref prob_cer prob_par
    #    0    0     ++       cer    cer        cer      .99       .1
    #    .
    #    .
    #    .

    categories = sorted(introgressed.keys())
    prob_types = sorted(probs.keys())

    if header:
        f.write('rep\tsite\tcoding')
        for category in categories:
            f.write('\t' + category)
        f.write('\t' + 'ref')
        for prob_type in prob_types:
            for state in probs[prob_type].keys():
                f.write('\t' + 'prob_' + prob_type + '_' + state)
        f.write('\n')

    for i in range(len(codings)):
        f.write(str(rep) + '\t' + str(i) + '\t' + codings[i])
        for category in categories:
            f.write('\t' + introgressed[category][i])
        f.write('\t' + introgressed_ref[actual_block_type][i])
        for prob_type in prob_types:
            for state in probs[prob_type].keys():
                try:
                    f.write('\t' + str(probs[prob_type][state][i]))
                except:
                    print prob_type, state, i, len(probs[prob_type][state])
                    sys.exit()
        f.write('\n')

def write_combined_files(files, inds, rep, codings, introgressed, probs, \
                         actual_block_type, ref_ind, header):
    for i in inds:
        write_combined_file(files[i], rep, codings[i], introgressed[i], \
                            actual_block_type, introgressed[ref_ind], probs[i], header)

def write_coding_table(seqs_coded, files, rep, header):
    
    for i in files.keys():
        f = files[i]
        if header:
            f.write('site\tcode\trep\n')
        for j in range(len(seqs_coded[i])):
            f.write(str(j) + '\t' + seqs_coded[i][j] + '\t' + str(rep) + '\n')

def write_one_block_set(blocks_dic, ind, f, rep, suffix = ''):

    for block_type in blocks_dic[ind].keys():
        for species in blocks_dic[ind][block_type].keys():
            for block in blocks_dic[ind][block_type][species]:
                f.write(block_type + suffix + '\t' + species + '\t' + \
                        str(block[0]) + '\t' + str(block[1]) + '\t' + \
                        str(rep) + '\n')

def write_blocks_table(blocks_dic, files, rep, ref_ind, header):

    for ind in blocks_dic.keys():
        f = files[int(ind)]
        if header:
            f.write('type\tspecies\tstart\tend\trep\n')
        # write blocks for the current individual
        write_one_block_set(blocks_dic, ind, f, rep)
        # and also the reference individual
        write_one_block_set(blocks_dic, ref_ind, f, rep, suffix='_ref')
        
