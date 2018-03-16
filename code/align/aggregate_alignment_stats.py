import os
import sys
sys.path.insert(0, '..')
import global_params as gp

gp_dir = '../'
stats_files = [gp_dir + gp.alignments_dir + x for x in filter(\
        lambda x: 'stats' in x and 'summary' not in x, os.listdir(gp_dir + gp.alignments_dir))]

# goal is to generate file for R:
# chromosome strain frac_S288c_S288c frac_S288c_CBS432 frac_S288c_x frac_CBS432_S288c frac_CBS432_CBS432 frac_CBS432_x frac_x_S288c frac_x_CBS432 frac_x_x  aligned_length_S288c aligned_length_CBS432 aligned_length_x num_align_columns_0 num_align_columns_1 num_align_columns_2 num_align_columns_3

f = open(gp_dir + gp.alignments_dir + 'mafft_stats_summary.txt', 'w')

f.write('chromosome\tstrain\tfrac_S288c_S288c\tfrac_S288c_CBS432\tfrac_S288c_x\tfrac_CBS432_S288c\tfrac_CBS432_CBS432\tfrac_CBS432_x\tfrac_x_S288c\tfrac_x_CBS432\tfrac_x_x\t\taligned_length_S288c\taligned_length_CBS432\taligned_length_x\tnum_align_columns_0\tnum_align_columns_1\tnum_align_columns_2\tnum_align_columns_3\n')

# one line for each of these files
for fn in stats_files:
    print fn
    
    lines = [line.strip() for line in open(fn, 'r').readlines()]

    sc = lines[7].split(',')[0]
    sp = lines[8].split(',')[0]
    sx = lines[9].split(',')[0]

    c0 = float(lines[1].split(',')[1])
    c1 = float(lines[2].split(',')[1])
    c2 = float(lines[3].split(',')[1])
    c3 = float(lines[4].split(',')[1])

    lc = float(lines[7].split(',')[1])
    lp = float(lines[8].split(',')[1])
    lx = float(lines[9].split(',')[1])

    fc = [float(x) for x in lines[12].split(',')[1:]]
    fp = [float(x) for x in lines[13].split(',')[1:]]
    fx = [float(x) for x in lines[14].split(',')[1:]]

    i = fn.find('chr')
    m = fn[i+3:]
    i = m.find('_')
    chrm = m[:i]

    f.write(chrm + '\t' + sx + \
            '\t' + str(fc[0]) + '\t' + str(fc[1]) + '\t' + str(fc[2]) + \
            '\t' + str(fp[0]) + '\t' + str(fp[1]) + '\t' + str(fp[2]) + \
            '\t' + str(fx[0]) + '\t' + str(fx[1]) + '\t' + str(fx[2]) + \
            '\t' + str(lc) + '\t' + str(lp) + '\t' + str(lx) + 
            '\t' + str(c0) + '\t' + str(c1) + '\t' + str(c2) + '\t' + str(c3) + '\n')
f.close()
