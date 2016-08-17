# take fastq files containing reads and quality information, along
# with reference genome, and convert to fasta file ... or vcf file and then fasta??

import os
import sys

fastq_dir = '/net/dunham/vol2/Giang/DunhamBeer/DunhamBeer'

quality_chars = list('!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~')
char_to_score = dict(zip(quality_chars, range(1, len(quality_chars))))

fastq_dir = '/net/dunham/vol2/Giang/DunhamBeer/DunhamBeer/'
ls = os.listdir(fastq_dir)
fns = []
for l in ls:
    if '.1.fastq' in l and 'stats' not in l and l[0] != 'N':
        fns.append(l[:-8])

ref_fasta = '/net/akey/vol2/aclark4/nobackup/100_genomes/genomes/S288c_SGD-R64.fa'

#####
# align reads with bwa
#####

samdir = '/net/akey/vol2/aclark4/nobackup/introgression/data/beer/dunham/sam/'
os.system('module load bwa/latest') # this doesn't actually work because it makes a new shell instance every time - TODO fix this
cmd = 'bwa index ' + ref_fasta 
#print cmd
#os.system(cmd)
for fn in fns:
    cmd = 'bwa mem ' + ref_fasta + ' ' + fastq_dir + fn + '.1.fastq ' + fastq_dir + fn + '.2.fastq' + ' > ' + samdir + fn + '.sam'
    print cmd
    os.system(cmd)

sys.exit()

#####
# run base recalibrator
#####

outdir = '/net/akey/vol2/aclark4/nobackup/introgression/data/beer/dunham/fasta/'
for fn in fns:
    # -knownSites database of previously known polymorphisms
    os.system('java -jar ~/software/GenomeAnalysisTK.jar -T BaseRecalibrator -R ' + ref_fasta + ' -I ' + fastq_dir + fn + ' -o ' + outdir + fn[:-1] + 'a')
    
#####
# run 
#####
    
