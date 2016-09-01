import os

outdir = '../../results/sim/'
fns = os.listdir(outdir)
newoutdir = '../../results/sim/fixed/'
for fn in fns:
    if 'summary' in fn:
        f = open(outdir + fn, 'r')
        line = f.readline()
        ncols = len(line.split('\t'))
        fnew = open(newoutdir + fn + '.txt', 'w')
        fnew.write(line)
        line = f.readline()
        line = line.split('\t')
        for i in range(len(line)):
            fnew.write(line[i])
            if i != 0 and i%(ncols-1) == 0:
                fnew.write('\n')
            else:
                fnew.write('\t')
        fnew.close()
        f.close()
