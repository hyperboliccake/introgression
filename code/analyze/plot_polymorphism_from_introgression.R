# also for plotting nucleotide diversity!

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)
source('../my_color_palette.R')

tag = 'u3_i.001_tv_l1000_f.01'

a = read.table('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/polymorphism/nucleotide_diversity.txt', header=T, stringsAsFactors=F)

a$chromosome = factor(a$chromosome, levels = c('all', 'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI'))

a1 = a[which(a$chromosome != 'all'),]
average_frac = 1 - a[which(a$chromosome=='all'),]$pi_nonint / a[which(a$chromosome=='all'),]$pi

ggplot(a1, aes(x=chromosome, y=1-pi_nonint/pi, fill='x')) +
    geom_bar(stat='identity', position='dodge') +
    scale_fill_manual(values = c(my_color_palette[['introgressed']])) +
    ylab('Fraction of nucleotide diversity resulting from paradoxus introgression') +
    xlab('Chromosome') +
    geom_hline(yintercept=average_frac, linetype='dashed', color='black') +
    scale_y_continuous(expand = c(0,0), limits=c(0,1.1*max(1-a1$pi_nonint/a1$pi)))+ 
        theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          legend.position="none",
          axis.line=element_line(),
          axis.title.x = element_text(size=16), 
          axis.title.y = element_text(size=15), 
          axis.text.x = element_text(colour="black",size=12), 
          axis.text.y = element_text(colour="black",size=12),
          legend.text=element_text(size=14),
          legend.title=element_blank())

ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/polymorphism/nuc_div_by_chrm.pdf', height=7, width=12)



a = read.table('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/polymorphism/polymorphism_summary.txt', header=T, stringsAsFactors=F)

a$chromosome = factor(a$chromosome, levels = c('all', 'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI'))

a1 = a[(which((a$sites == 'biallelic')&(a$match == 'all')&(a$chromosome != 'all'))),]
average_frac = a[which((a$chromosome=='all')&(a$sites == 'biallelic')&(a$match == 'all')),]$frac

ggplot(a1, aes(x=chromosome, y=frac, fill='x')) +
    geom_bar(stat='identity', position='dodge') +
    scale_fill_manual(values = c(my_color_palette[['introgressed']])) +
    ylab('fraction of biallelic sites resulting from paradoxus introgression') +
    geom_hline(yintercept=average_frac, linetype='dashed', color='black') +
    scale_y_continuous(expand = c(0,0), limits=c(0,1.1*max(a1$frac)))+ 
        theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          legend.position="none",
          axis.line=element_line(),
          axis.title.x = element_text(size=16), 
          axis.title.y = element_text(size=16), 
          axis.text.x = element_text(colour="black",size=12), 
          axis.text.y = element_text(colour="black",size=12),
          legend.text=element_text(size=14),
          legend.title=element_blank())


#b = a[which(a$sites == 'biallelic'),]
#ggplot(b, aes(x=chromosome, y=frac, fill=match)) +
#    geom_bar(stat='identity', position='dodge') +
#    ylab('fraction of biallelic sites in cerevisiae resulting from paradoxus introgression')

ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/polymorphism/polymorphism_by_chrm_biallelic.pdf', height=7, width=12)

stop here


b = a[which(a$sites == 'polymorphic'),]
ggplot(b, aes(x=chromosome, y=frac, fill=match)) +
    geom_bar(stat='identity', position='dodge')+
    ylab('fraction of polymorphism in cerevisiae resulting from paradoxus introgression')

ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/polymorphism/polymorphism_by_chrm.pdf', height=7, width=12)


h = read.table('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/count_introgressed.txt', header=T, sep="\t", stringsAsFactors=F)

h$chromosome = factor(h$chromosome, levels = c('all', 'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI'))

ggplot(h, aes(x=chromosome, y=at_least_one_frac)) +
    geom_bar(stat='identity') +
    ylab('fraction introgressed in at least one strain')

ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/polymorphism/introgression_by_chrm.pdf', height=7, width=12)

# polymorphism plus total biallelic sites in one plot for comparison porpoises

a1 = a[(which((a$sites == 'biallelic')&(a$match == 'all'))),]
d1 = data.frame(chromosome = a1$chromosome, type = 'biallelic due to introgression from paradoxus', frac = a1$frac)
d2 = data.frame(chromosome = h$chromosome, type = 'introgressed in at least one strain', frac = h$at_least_one_frac)
d = rbind(d1, d2)

ggplot(d, aes(x=chromosome, y=frac, fill=type)) +
    geom_bar(stat='identity', position='dodge') +
    ylab('fraction of sites') +
    scale_y_continuous(expand = c(0,0), limits=c(0,1.1*max(d$frac)))+ 
        theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(colour="black",size=12), 
          axis.text.y = element_text(colour="black",size=12),
          legend.text=element_text(size=14),
          legend.title=element_blank())


ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/polymorphism/introgression_polymorphism_by_chrm.pdf', height=6, width=12)


# coding vs nonconcoding, synonymous vs nonsynonymous

a = read.table('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/coding_changes_summary_u3_i.001_tv_l1000_f.01.txt', header=T, stringsAsFactors=F)

b = a[which(a$label_type == 'all'),]

d = data.frame('coding'=c('coding', 'coding', 'coding changes', 'coding changes', 'noncoding', 'noncoding'),
               'category'=c('coding', 'ref_gene_only', 'syn', 'non', 'noncoding', 'strain_orf_only'),
               'count'=c(b[which(b$'total_type'=='coding'),]$total,
                         b[which(b$'total_type'=='ref_gene_only'),]$total,
                         b[which(b$'total_type'=='syn'),]$total,
                         b[which(b$'total_type'=='non'),]$total,
                         b[which(b$'total_type'=='noncoding'),]$total,
                         b[which(b$'total_type'=='strain_orf_only'),]$total))


ggplot(d, aes(x=coding, y=count, fill=category)) +
    geom_bar(stat='identity', position='dodge') +
    ylab('number of sites')

ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/plots/count_coding.pdf', height=7, width=12)


aw = dcast(a, label + label_type ~ total_type, value.var='total')
aw$frac_coding = aw$coding / (aw$coding + aw$noncoding)
aw$frac_non = aw$non / (aw$non + aw$syn)
al = melt(aw, id.vars = c('label', 'label_type'))

d = al[which((al$label_type == 'strain')&(al$variable == 'frac_coding')),]

ggplot(d, aes(x=reorder(label, -value), y=value, fill='x')) + 
    geom_bar(stat='identity',position='dodge') + 
    xlab('strain') + ylab('fraction coding') +
    scale_fill_viridis(discrete=TRUE) +
    guides(fill=FALSE) +
    scale_y_continuous(expand = c(0,0), limits=c(0,max(d$value))) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(angle = 45,vjust = 1,hjust=1,colour="black"), 
          axis.text.y = element_text(colour="black"))
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/frac_coding_by_strain.pdf',sep=''), width = 12, height = 6)


d = al[which((al$label_type == 'gene')&(al$variable == 'frac_non')),]

ggplot(d, aes(x=reorder(label, -value), y=value, fill='x')) + 
    geom_bar(stat='identity',position='dodge') + 
    xlab('gene') + ylab('fraction nonsynonymous') +
    scale_fill_viridis(discrete=TRUE) +
    guides(fill=FALSE) +
    scale_y_continuous(expand = c(0,0), limits=c(0,max(d$value))) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(angle = 45,vjust = 1,hjust=1,colour="black",size=1), 
          axis.text.y = element_text(colour="black"))
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/frac_nonsynonymous_by_gene.pdf',sep=''), width = 12, height = 6)



jgfjhf













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


