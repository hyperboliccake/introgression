import os

outdir = '../../results/sim/run_3/'
fns = os.listdir(outdir)
newoutdir = '../../results/sim/run_3/fixed/'
for fn in fns:
    if 'summary' in fn and 'txt' not in fn:
        f = open(outdir + fn, 'r')
        line = f.readline()
        ncols = len(line.strip().split('\t'))
        fnew = open(newoutdir + fn + '.txt', 'w')
        fnew.write(line)
        line = f.readline()
        line = line.strip().split('\t')
        for i in range(len(line)):
            fnew.write(line[i])
            if (i+1)%ncols == 0:
                fnew.write('\n')
            else:
                fnew.write('\t')
        fnew.close()
        f.close()
