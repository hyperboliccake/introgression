import os
import sys
sys.path.insert(0, '..')
import global_params as gp
import numpy

gp_dir = '../'
stats_files = [gp_dir + gp.alignments_dir + x for x in filter(\
        lambda x: 'stats' in x, os.listdir(gp_dir + gp.alignments_dir))]

#avg_frac_aligned_by_chrm = dict(zip(gp.chrms, [0]*len(gp.chrms)))
avg_frac_aligned_p = 0
avg_frac_aligned_x = 0
total_p = 0
total_x = 0

a = []

for fn in stats_files:
    lines = [line.strip() for line in open(fn, 'r').readlines()]
    lc = float(lines[7].split(',')[1])
    lp = float(lines[8].split(',')[1])
    lx = float(lines[9].split(',')[1])
    fp, fx = [float(x) for x in lines[12].split(',')[2:]]
    avg_frac_aligned_p += fp * lp
    avg_frac_aligned_x += fx * lx
    total_p += lp
    total_x += lx
    #print fn[fn.find('chr')-8:], fx, lx, lc
    a.append(fx)


avg_frac_aligned_p /= total_p
avg_frac_aligned_x /= total_x

print len(stats_files)
print avg_frac_aligned_p
print avg_frac_aligned_x

hist, edges = numpy.histogram(a, bins=30)
print hist
print edges
print sum(hist[:-1])
