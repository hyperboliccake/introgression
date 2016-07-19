# given a list of introgressed regions in this format:
# Sigma1278b.chrX, + strand, 9910-10909
#
# generate a set of annotations in the introgressed/ folder:
#
# for each gene that is called introgressed in at least one strain,
# created folder introgressed/gene/ containing the alignment
# of cerevisiae and paradoxus references to all the introgressed
# versions (gene_introgressed.fasta), and also to all of the versions
# (gene_all.fasta).
#
# for each region called introgressed, generate a file in introgressed/regions/
# that is the alignment for the two references and the strain with the
# introgressed region, S288c_CBS432_strain_chr_start-end.fasta. In this file, the aligned
# bases within coding sequence are lower case. In addition, there is a corresponding
# file S288c_CBS432_strain_chr_start-end.genes.txt listing the genes that overlap
# with this region, and the indices of the bases they overlap, in this format:
# gene_name, 0-149, 25236-25385
# gene_name, 200-600, ....


import re

def reverse_complement(s):
    t = ''
    for i in s:
        if i == 'A':
            t += 'T'
        elif i == 'T':
            t += 'A'
        elif i == 'G':
            t += 'C'
        elif i == 'C':
            t += 'G'
        else:
            t += i
    return t[::-1]

translate = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
       "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
       "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
       "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
       "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
       "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
       "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
       "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
       "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
       "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
       "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

ref_cer = 'S288c'
ref_par = 'CBS432'


##### 
# For each introgressed region, extract relevant part of alignment to a separate file
#####

suffix = 'id'

# read in introgressed regions
lines = [x.split(',') for x in open('introgressed_' + suffix + '.txt', 'r').readlines()]
regions = {} # introgressed regions keyed by strain and then chromosome
for line in lines:
    strain = line[0][:line[0].find('_')].lower()
    chrm = re.search(r'chr(?P<chrm>[IVXM]+)', line[0]).group('chrm')
    i = line[2].find('-')
    # +/-, start, end
    entry = [line[1][1], int(line[2][:i]), int(line[2][i+1:]), int(line[3]), int(line[4])]
    if strain in regions:
        if chrm in regions[strain]:
            regions[strain][chrm].append(entry)
        else:
            regions[strain][chrm] = [entry]
    else:
        regions[strain] = {}
        regions[strain][chrm] = [entry]

# pull out alignment blocks
for strain in regions.keys():
    for chrm in regions[strain]:
        for entry in regions[strain][chrm]:
            print strain, chrm, entry
            f = open('alignments_gb/S288c_CBS432_' + strain + '_chr' + chrm + '.maf', 'r')
            line = f.readline()        
            block = []
            relative_start = -1
            relative_end = -1
            while line != '':
                if line[0] == 'a':
                    block = [line]
                else:
                    block.append(line)
                # matches from beginning of line
         
                if re.match('s ' + strain + '_chr' + chrm + '.', line) != None:
                    line = line.split()
                    # subtract one to index by zero; also add in some extra context
                    #context = 500
                    start = int(line[2]) - 1
                    #start_with_context = max(0, start - context)
                    # for end, need to ignore gaps in the
                    # non-reference sequence, and subtract one to make
                    # inclusive; note that using the given length in
                    # the alignment file only counts non-gaps (i.e. is
                    # the actual sequence coordinates)
                    end = start + int(line[3]) - 1
                    #end_with_context = min(len(block[0].split()[6]), end + context)
                    # okay, so alignment blocks can be overlapping
                    # which is unexpected, except I _think_ this only
                    # happens when they're on opposite strands? still
                    # weird though...anyway we'll just directly pull the
                    # same alignment block we've notated
                    #print entry[3]-1, entry[4]-1, entry[0], start, int(line[3])-1
                    if entry[3]-1 == start and entry[4]-1 == int(line[3])-1 and entry[0] == line[4]:
                    #if entry[1] >= start and entry[1] <= end:
                        #assert entry[2] >= start and entry[2] <= end, strain + ' ' + chrm + ' ' + str(entry) + ' ' + str(start) + ' ' + str(end) + ', ' + str(entry[1]) + ' ' + str(entry[2]) + '\n' + str(line)
                        assert entry[1] >= start and entry[1] <= end and entry[2] >= start and entry[2] <= end
                        # relative in alignment block, factoring in gaps
                        relative_start = entry[1] - start
                        i = 0
                        a = 0
                        while i + a < len(line[6]) and i < relative_start:
                            if line[6][i+a] == '-':
                                a += 1
                            else:
                                i += 1
                        relative_start += a

                        relative_end = entry[2] - start
                        i = 0
                        a = 0
                        while i + a < len(line[6]) and i < relative_end:
                            if line[6][i+a] == '-':
                                a += 1
                            else:
                                i += 1
                        relative_end += a

                        break
                line = f.readline()
            f.close()

            assert relative_start != -1, strain + ' ' + chrm + ' ' + str(entry)
            fout = open('introgressed/regions/S288c_CBS432_' + strain + '_chr' + chrm + '_' + str(entry[1]) + '-' + str(entry[2]) + '.maf', 'w')
            assert block[0][0] == 'a'
            block = block[1:]
            # use the same coordinates for all of them because we just
            # want to pull out the parts of the references that align
            # to the same part of the current strain

            # add some extra positions onto either end for context
            # might end up being less than 500 if the alignment block isn't that long
            context = 500
            #relative_start_with_context = max(0, relative_start - context)
            #relative_end_with_context = min(len(block[0].split()[6]), relative_end + context)
            relative_start_with_context = relative_start
            r1 = 0
            sx = block[-1].split()[6]
            while r1 < context and relative_start_with_context > 0:
                if sx[relative_start_with_context] != '-':
                    r1 += 1
                relative_start_with_context -= 1

            relative_end_with_context = relative_end
            r2 = 0
            while r2 < context and relative_end_with_context < len(sx):
                if sx[relative_end_with_context] != '-':
                    r2 += 1
                relative_end_with_context += 1
            print '*****', relative_start_with_context, relative_start, relative_end, relative_end_with_context
            if len(block) == 3:
                b = block[0].split()
                s = b[6]
                fout.write('>' + ref_cer + ' ' + chrm + ' ' + str(int(b[2]) + relative_start_with_context - s[:relative_start_with_context].count('-')) + ' ' + str(int(b[2]) + relative_end_with_context - s[:relative_end_with_context].count('-')) + '\n')# + ' ' + str(r1) + ' ' + str(r2) + '\n')
                fout.write(s[relative_start_with_context:relative_start] + '|')
                fout.write(s[relative_start:relative_end+1])
                fout.write('|' + s[relative_end+1:relative_end_with_context] + '\n')
                block = block[1:]

            b = block[0].split()
            s = b[6]
            fout.write('>' + ref_par + ' ' + chrm + ' ' + str(int(b[2]) + relative_start_with_context - s[:relative_start_with_context].count('-')) + ' ' + str(int(b[2]) + relative_end_with_context - s[:relative_end_with_context].count('-')) + '\n')
            fout.write(s[relative_start_with_context:relative_start] + '|')
            fout.write(s[relative_start:relative_end+1])
            fout.write('|' + s[relative_end+1:relative_end_with_context] + '\n')
            block = block[1:]

            b = block[0].split()
            s = b[6]
            fout.write('>' + strain + ' ' + chrm + ' ' + str(int(b[2]) + relative_start_with_context - s[:relative_start_with_context].count('-')) + ' ' + str(int(b[2]) + relative_end_with_context - s[:relative_end_with_context].count('-')) + '\n')
            fout.write(s[relative_start_with_context:relative_start] + '|')
            fout.write(s[relative_start:relative_end+1])
            fout.write('|' + s[relative_end+1:relative_end_with_context] + '\n')
            fout.close()

#####
# Identify genes (if any) that overlap each alignment block
#####
print regions
f = open('/net/akey/vol2/aclark4/nobackup/100_genomes/sequence.gb', 'r')
line = f.readline()
strain = ''
chrm = ''
in_gene = False
gene_inds = []
all_genes = {}
strand = '+'
while line != '':
    # starting with a new species x chromosome
    if line[:10] == 'DEFINITION':
        print line
        m = re.search('Saccharomyces cerevisiae (?P<strain>[a-zA-Z0-9]+) chromosome (?P<chrm>[IVXM]+)', line)
        strain = m.group('strain').lower()
        chrm = m.group('chrm')
    # found a gene block
    elif len(line.split()) > 1 and line.split()[0] == 'gene':
        # if we never found a gene name for the previous one, forget about it
        in_gene = False
        gene_inds = []
        # extract start and end coordinates, and check whether they
        # fall within any of the introgressed regions
        m = re.search(r'[><]?(?P<start>[0-9]+)[.><,0-9]*\.\.[><]?(?P<end>[0-9]+)', line)
        # subtract one to index from zero
        start = int(m.group('start')) - 1
        end = int(m.group('end')) - 1
        if strain in regions and chrm in regions[strain]:
            for i in range(len(regions[strain][chrm])):
                l = regions[strain][chrm][i]
                if (l[1] > start and l[1] < end) or (l[2] > start and l[2] < end):
                    in_gene = True
                    # this is a list just in case the gene happens to
                    # overlap multiple introgressed regions
                    gene_inds.append(i)
                    if 'complement' in line:
                        strand = '-'
                    else:
                        strand = '+'
    # found name line for the gene
    elif in_gene and '/gene' in line:
        gene_name = line[line.find('/gene="')+7:-2]        
        if gene_name in all_genes:
            all_genes[gene_name].append(strain)
        else:
            all_genes[gene_name] = [strain]
        for i in gene_inds:
            regions[strain][chrm][i].append([strand, start, end, gene_name])
        in_gene = False
        gene_inds = []
    # put gene sequences into entry
    elif line[:6] == 'ORIGIN' and strain in regions and chrm in regions[strain]:
        seq = ''
        line = f.readline()
        while line != '//\n':
            for x in line.split()[1:]:
                seq += x
            line = f.readline()
        # for this strain and chrm, given region r
        for r in range(len(regions[strain][chrm])):
            # each gene overlapping that region (possibly none)
            for i in range(5, len(regions[strain][chrm][r])):
                regions[strain][chrm][r][i].append(seq[regions[strain][chrm][r][i][1]:regions[strain][chrm][r][i][2]+1])

    line = f.readline()
f.close()

#print regions
#print all_genes


#####
# Create a modified version of each introgressed alignment file that
# has lowercase letters for positions within genes; also create a file
# for each introgressed region that lists the genes within it
#####

startf = ['atg']
stopf = ['taa', 'tag', 'tga']
startr = ['cat']
stopr = ['tta', 'cta', 'tca']

all_genes_fns = {}

for strain in regions.keys():
    for chrm in regions[strain]:
        print regions[strain][chrm]
        for entry in regions[strain][chrm]:
            
            f = open('introgressed/regions/S288c_CBS432_' + strain + '_chr' + chrm + '_' +  str(entry[1]) + '-' + str(entry[2]) + '.maf', 'r')
            lines = f.readlines()
            headers = lines[::2]
            h = headers[0].split()
            #context_before = int(h[-2])
            #context_after = int(h[-1])
            #relative_start_with_context = int(h[-4])
            #relative_start = int(h[-3])
            #relative_end = int(h[-2])
            #relative_end_with_context = int(h[-1])
            seqs_with_context = [s[:-1].upper() for s in lines[1::2]]
            context_before_ind = seqs_with_context[0].find('|')
            context_before = context_before_ind - seqs_with_context[-1][:context_before_ind].count('-')
            context_after_ind = seqs_with_context[0].rfind('|')
            context_after = len(seqs_with_context[-1]) - 1 - context_after_ind - seqs_with_context[-1][context_after_ind:].count('-')
            seqs = [x.replace('|', '') for x in seqs_with_context]
            #seqs_context_before = []
            #seqs_context_after = []
            #for s in seqs_with_context:
            #    seqs.append(s[s.find('|')+1:s.rfind('|')])
            #    seqs_context_before.append(s[:s.find('|')])
            #    seqs_context_after.append(s[s.rfind('|')+1:])

            f.close()
            fn = 'introgressed/regions/S288c_CBS432_' + strain + '_chr' + chrm + '_' + str(entry[1]) + '-' + str(entry[2]) + '_annotated.maf'
            f_mod = open(fn, 'w')
            f_genes = open('introgressed/regions/S288c_CBS432_' + strain + '_chr' + chrm + '_' + str(entry[1]) + '-' + str(entry[2]) + '_genes.maf', 'w')
            #seqs_annotated = []
            for gene in entry[5:]:
                if gene[3] in all_genes_fns:
                    all_genes_fns[gene[3]].append(fn)
                else:
                    all_genes_fns[gene[3]] = [fn]
                print '============'
                # gene is position of gene in genome
                # entry is position of introgressed region
                # gene[1&2] are 1-indexed, entry[1&2] are 0-indexed
                relative_start = max(0, gene[1] - entry[1] + context_before - 1)
                i = 0
                a = 0
                while i + a < len(seqs[-1]) and i < relative_start:
                    if seqs[-1][i+a] == '-':
                        a += 1
                    else:
                        i += 1
                # advance by number of gaps because those aren't
                # included in position
                relative_start += a
                relative_start = min(relative_start, len(seqs[-1])-1)
                print a
                # inclusive end
                relative_end = min(len(seqs[-1])-1, gene[2] - entry[1] + context_before)
                i = 0
                a = 0
                while i + a < len(seqs[-1]) and i < relative_end:
                    if seqs[-1][i+a] == '-':
                        a += 1
                    else:
                        i += 1
                relative_end += a
                relative_end = min(relative_end, len(seqs[-1])-1)
                print a

                # TODO finish context stuff; especially figure out how to index gene within context, taking gaps into account
                for s in range(len(seqs)):
                    seqs[s] = seqs[s][:relative_start] + seqs[s][relative_start:relative_end].lower() + seqs[s][relative_end:]
                    try:
                        if s == len(seqs) - 1 and relative_start != 0:
                            if gene[0] != entry[0]:
                                assert seqs[s][relative_start:relative_start+3] in stopr, seqs[s] + '\n' + seqs[s][relative_start:relative_start+3] + ' ' + \
                                    str(entry[1]) + ' ' +  str(entry[2]) + ' ' + str(gene[1]) + ' ' + str(gene[2]) + ' ' + str(relative_start) + ' ' + str(relative_end)
                            else:
                                assert seqs[s][relative_start:relative_start+3] in startf, seqs[s] + '\n' +  seqs[s][relative_start:relative_start+3] + ' ' + \
                                    str(entry[1]) + ' ' +  str(entry[2]) + ' ' + str(gene[1]) + ' ' + str(gene[2]) + ' ' + str(relative_start) + ' ' + str(relative_end)
                        if s == len(seqs) - 1 and relative_end != len(seqs[0]) - 1:
                            if gene[0] != entry[0]:
                                assert seqs[s][relative_end-3:relative_end] in startr, seqs[s] + '\n' + seqs[s][relative_end-3:relative_end] + ' ' + \
                                    str(entry[1]) + ' ' +  str(entry[2]) + ' ' + str(gene[1]) + ' ' + str(gene[2]) + ' ' + str(relative_start) + ' ' + str(relative_end)
                            else:
                                assert seqs[s][relative_end-3:relative_end] in stopf, seqs[s] + '\n' + seqs[s][relative_end-3:relative_end] + ' ' + \
                                    str(entry[1]) + ' ' +  str(entry[2]) + ' ' + str(gene[1]) + ' ' + str(gene[2]) + ' ' + str(relative_start) + ' ' + str(relative_end)
                    except Exception, e:
                        print e
                        print strain, chrm
                        print gene
                        print gene[0], entry[0]
                        print gene[-1][:3], gene[-1][-3:]
                        continue
                f_genes.write(gene[3] + ', ')
                # complement?
                if gene[0]:
                    f_genes.write('-')
                else:
                    f_genes.write('+')
                # the gene start and end here are 1-indexed
                f_genes.write(', ' + str(relative_start) + '-' + str(relative_end) + ', ' + str(gene[1]) + '-' + str(gene[2]) + '\n')
            f_genes.close()
            for i in range(len(headers)):
                f_mod.write(headers[i][:-1])
                # name, start
                g = [(x[3], x[1]) for x in entry[5:]]
                g.sort(key = lambda x: x[1])
                for gi in g:
                    f_mod.write(' ' + gi[0])
                f_mod.write('\n')
                f_mod.write(seqs[i][:context_before_ind] + '|' + \
                                seqs[i][context_before_ind:context_after_ind-1] + '|' + \
                                seqs[i][context_after_ind-1:] + '\n')
            f_mod.close()
            
#####
# gene x strain and strain x gene files
#####

f = open('introgressed_id_genes.txt', 'w')
for gene in all_genes:
    f.write(gene)
    for strain in all_genes[gene]:
        f.write(' ' + strain)
    f.write('\n')
f.close()

f = open('introgressed_id_genes_fns.txt', 'w')
for gene in all_genes_fns:
    f.write(gene)
    for fn in all_genes_fns[gene]:
        f.write(' ' + fn)
    f.write('\n')
f.close()

f = open('introgressed_id_strains.txt', 'w')
for strain in regions:
    num_regions = 0
    num_bp = 0
    f.write(strain + ' ')
    for chrm in regions[strain]:
        for entry in regions[strain][chrm]:
            num_regions += 1
            num_bp += entry[2] - entry[1] + 1
            for gene in entry[5:]:
                f.write(gene[3] + ' ')
            
    # total number of introgressed regions
    f.write(str(num_regions) + ' ')
    # total number of introgressed bases
    f.write(str(num_bp) + '\n')
f.close()
    




