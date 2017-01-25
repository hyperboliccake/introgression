import math
import numpy.random
import sys
sys.path.insert(0, '../misc/')
import mystats

#TODO new introgressed file format, order of entries etc

# distribution of gaps between adjacent introgressed sequeneces -->
# want to get an idea of whether we're not stitching nearby ones
# together
def get_gap_hist(lines):

    d = {}
    for line in lines:
        # line[0] example: yjm1244_chrI.yjm1244
        strainchrm, strain = line[0].split('.')
        chrm = strainchrm.split('_')[1]
        if strain not in d:
            d[strain] = {}
        if chrm not in d[strain]:
            d[strain][chrm] = []
        start = line[2]
        end = line[3]
        # TODO deal w reverse strand better?
        d[strain][chrm].append((int(start), int(end)))

    gap_hist = []

    for strain in d:
        for chrm in d[strain]:
            entries = sorted(list(set(d[strain][chrm])), key = lambda x: x[0])
            for i in range(1, len(entries)):
                # TODO WHY???
                if entries[i][0] - entries[i-1][1] < 0:
                    print '********'
                    print entries[i]
                    print entries[i-1]
                else:
                    gap_hist.append(entries[i][0] - entries[i-1][1])

    f_out = open('../../results/summary_stats_gaps.txt', 'w')

    for g in gap_hist:
        f_out.write(str(g) + ' ')
    f_out.write('\n')

    f_out.write(str(mystats.mean(gap_hist)) + ' ' + str(mystats.std_err(gap_hist)) + ' ')
    bs = mystats.bootstrap(gap_hist)
    f_out.write(str(bs[0]) + ' ' + str(bs[1]) + '\n')

    f_out.close()

# lengths of introgressed regions
def get_length_hist(lines):

    length_hist = []
    for line in lines:
        # line[0] example: yjm1244_chrI.yjm1244
        start, end = line[2].split('-')
        length_hist.append(int(end) - int(start) + 1)
            
    f_out = open('../../results/summary_stats_lengths.txt', 'w')

    for l in length_hist:
        f_out.write(str(l) + ' ')
    f_out.write('\n')

    f_out.write(str(mystats.mean(length_hist)) + ' ' + str(mystats.std_err(length_hist)) + ' ')
    bs = mystats.bootstrap(length_hist)
    f_out.write(str(bs[0]) + ' ' + str(bs[1]) + '\n')

    f_out.close()


# average length of introgressed sequence in different strains
def get_lengths_by_strain(lines):

    d = {}
    for line in lines:
        # line[0] example: yjm1244_chrI.yjm1244
        strainchrm, strain = line[0].split('.')
        if strain not in d:
            d[strain] = []
        start, end = line[2].split('-')
        d[strain].append(int(end) - int(start) + 1)

    f_out = open('../../results/summary_stats_lengths_by_strain.txt', 'w')
    for strain in d:
        f_out.write(strain + ' ')
        f_out.write(str(mystats.mean(d[strain])) + ' ' + str(mystats.std_err(d[strain])) + ' ')
        bs = mystats.bootstrap(d[strain])
        f_out.write(str(bs[0]) + ' ' + str(bs[1]))
        #for r in d[strain]:
        #    f_out.write(' ' + str(r))
        f_out.write('\n')

    f_out.close()

# total amount of introgressed sequence in different strains
def get_total_length_by_strain(lines):

    d = {}
    for line in lines:
        # line[0] example: yjm1244_chrI.yjm1244
        strainchrm, strain = line[0].split('.')
        if strain not in d:
            d[strain] = []
        start, end = line[2].split('-')
        d[strain].append(int(end) - int(start) + 1)

    f_out = open('../../results/summary_stats_total_length_by_strain.txt', 'w')
    for strain in d:
        f_out.write(strain + ' ')
        f_out.write(str(sum(d[strain])))
        f_out.write('\n')

    f_out.close()

# get number of strains that each gene is introgressed in (given that
# it's introgressed in at least one strain)
def get_strains_by_gene(lines):
    f_out = open('../../results/summary_stats_strains_by_gene.txt', 'w')
    for line in lines:
        l = len(line) - 1 # subtract gene name column
        f_out.write(line[0] + ' ' + str(l) + '\n')
    f_out.close()

# get histogram of number of strains that each gene is introgressed in
# (given that it's introgressed in at least one strain)
def get_strains_by_gene_hist(lines):
    d = {}
    for line in lines:
        l = len(line) - 1 # subtract gene name column
        if l not in d:
            d[l] = 0
        d[l] += 1
    f_out = open('../../results/summary_stats_strains_by_gene.txt', 'w')
    for i in range(1, max(d) + 1):
        f_out.write(str(i) + ' ')
        if i in d:
            f_out.write(str(d[i]))
        else:
            f_out.write('0')
        f_out.write('\n')
    f_out.close()

# amount of introgressed sequence genic/intergenic

# something about how much introgressed sequence is due to gaps

lines = [x.strip().split(',') for x in open('../../results/introgressed.txt', 'r').readlines()]

get_gap_hist(lines)
sys.exit()
get_length_hist(lines)
get_lengths_by_strain(lines)
get_total_length_by_strain(lines)

lines = [x.strip().split(' ') for x in open('../../results/introgressed_id_genes.txt', 'r').readlines()]

get_strains_by_gene(lines)
