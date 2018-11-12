import os
import sys
sys.path.insert(0, '..')
import global_params as gp

gp_dir = '../'
stats_files = [gp_dir + gp.alignments_dir + x for x in filter(\
        lambda x: 'stats' in x and 'summary' not in x, os.listdir(gp_dir + gp.alignments_dir))]

# goal is to generate file for R (e.g. for two references and test strain):
# chromosome strain frac_S288c_S288c frac_S288c_CBS432 frac_S288c_x frac_CBS432_S288c frac_CBS432_CBS432 frac_CBS432_x frac_x_S288c frac_x_CBS432 frac_x_x  aligned_length_S288c aligned_length_CBS432 aligned_length_x num_align_columns_0 num_align_columns_1 num_align_columns_2 num_align_columns_3

f = open(gp_dir + gp.alignments_dir + 'mafft_stats_summary.txt', 'w')

f.write('chromosome\tstrain')

gp.alignment_ref_order

for ref1 in gp.alignment_ref_order + ['x']:
    for ref2 in gp.alignment_ref_order + ['x']:
        f.write('\t' + 'frac_' + ref1 + '_' + ref2)

for ref in gp.alignment_ref_order + ['x']:
    f.write('\t' + 'aligned_length_' + ref)

for i in range(0, len(gp.alignment_ref_order) + 2):
    f.write('\t' + 'num_align_columns_' + str(i))
    
f.write('\n')

all_strains = gp.alignment_ref_order + ['x']

# one line for each of these files
for fn in stats_files:
    print fn
    
    lines = [line.strip() for line in open(fn, 'r').readlines()]

    # histogram of number of number of strains aligned
    c = []
    offset = 1
    for i in range(len(all_strains) + 1):
        c.append(float(lines[i + offset].split(',')[1]))

    # aligned lengths
    l = []
    offset += len(all_strains) + 1 + 2
    for i in range(len(all_strains)):
        l.append(float(lines[i + offset].split(',')[1]))

    sx = lines[offset + len(all_strains) - 1].split(',')[0]

    # frac aligned to reference
    fr = []
    offset += len(all_strains) + 2
    for i in range(len(all_strains)):
        fr.append([float(x) for x in lines[i + offset].split(',')[1:]])

    i = fn.find('chr')
    m = fn[i+3:]
    i = m.find('_')
    chrm = m[:i]

    f.write(chrm + '\t' + sx)

    for i in range(len(all_strains)):
        for j in range(len(all_strains)):
            f.write('\t' + str(fr[i][j]))
    for i in range(len(all_strains)):
        f.write('\t' + str(l[i]))
    for i in range(len(all_strains) + 1):
        f.write('\t' + str(c[i]))
    f.write('\n')

f.close()
