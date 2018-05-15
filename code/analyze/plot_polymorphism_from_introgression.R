library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)

a = read.table('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/polymorphism/polymorphism_summary.txt', header=T, stringsAsFactors=F)

b = a[which(a$sites == 'biallelic'),]
ggplot(b, aes(x=chromosome, y=frac, fill=match)) +
    geom_bar(stat='identity', position='dodge') +
    ylab('fraction of biallelic sites in cerevisiae resulting from paradoxus introgression')

ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/polymorphism/polymorphism_by_chrm_biallelic.pdf', height=7, width=12)

b = a[which(a$sites == 'polymorphic'),]
ggplot(b, aes(x=chromosome, y=frac, fill=match)) +
    geom_bar(stat='identity', position='dodge')+
    ylab('fraction of polymorphism in cerevisiae resulting from paradoxus introgression')

ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/polymorphism/polymorphism_by_chrm.pdf', height=7, width=12)


h = read.table('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/count_introgressed.txt', header=T, sep="\t", stringsAsFactors=F)

ggplot(h, aes(x=chromosome, y=at_least_one_frac)) +
    geom_bar(stat='identity') +
    ylab('fraction introgressed in at least one strain')

ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/polymorphism/introgression_by_chrm.pdf', height=7, width=12)


axdsg


a$type = paste(a$number_cer_alleles, a$cer_ref_match_par_ref,
               a$one_match_intd, a$all_match_intd)

b = a[which(a$number_cer_alleles == '2'),]
total_biallelic = sum(b[which(b$chromosome =='all'),]$count)
b$frac = b$count / total_biallelic
print(total_biallelic)
b = b[which(b$cer_ref_match_par_ref == 'False'),]
b = b[which(b$one_match_intd == 'True'),]
d = aggregate(list(frac=b$frac), by=list(chromosome=b$chromosome), FUN=sum)
d$number_cer_alleles = '2'
d$cer_ref_match_par_ref = 'False'
d$one_match_intd = 'True'
d$all_match_intd = 'TrueFalse'
d$count = 'NA'
d$type = paste(d$number_cer_alleles, d$cer_ref_match_par_ref,
               d$one_match_intd, d$all_match_intd)
#b = rbind(b, d)

