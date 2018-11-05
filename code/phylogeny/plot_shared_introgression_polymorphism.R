library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)
source('../my_color_palette.R')

options(stringsAsFactors=FALSE)

a = read.table('shared_introgression_nonsingleton_polymorphism_3strains.txt', header=T, sep='\t')

d = a
d$group = paste(a$in_3strains, a$in_only_3strains)
dtt = d[which(d$group=='TRUE TRUE'),]
dtf = d[which(d$group=='TRUE FALSE'),]
dff = d[which(d$group=='FALSE FALSE'),]

dtt$group = 'in only 3 strains'
ndtt = nrow(dtt)
dtf$group = 'in 3 strains and others'
ndtf = nrow(dtf)
dff$group = 'not in 3 strains'
ndff = nrow(dff)

n = min(ndtt, ndtf, ndff)
dtts = dtt[sample(ndtt, n, replace=F),]
dtfs = dtf[sample(ndtf, n, replace=F),]
dffs = dff[sample(ndff, n, replace=F),]
dg = rbind(dtts, dtfs, dffs)

d = rbind(dtt, dtf, dff)


ggplot(dg, (aes(x=group, y=pi))) +
    geom_boxplot() + geom_jitter(width=.3, alpha=.5) +
    #xlab('group (in 3 strains & in only 3 strains)') +
    xlab('') +
    ylab('nucleotide diversity (pi)') 

ggsave('a.pdf', width = 12, height = 7)

ggplot(d, (aes(x=group, y=pi))) +
    geom_boxplot() + geom_jitter(width=.3, alpha=.5) +
    #xlab('group (in 3 strains & in only 3 strains)') +
    xlab('') +
    ylab('nucleotide diversity (pi)') 

ggsave('a2.pdf', width = 12, height = 7)



ggplot(d, aes(x=pi, fill=group)) + geom_histogram(aes(y =..scaled..), alpha=.3, position='identity')
    #geom_density(alpha=.3)
ggsave('b.pdf', width = 12, height = 7)


asdf

b = a[which(a$in_3strains == TRUE),]

ggplot(b, (aes(x=num_total, y=frac_poly, colour=in_only_3strains))) + geom_point(size=1, alpha=.5)+
    xlab('shared region size') +
    ylab('fraction polymorphic sites') + 
        theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          legend.title = element_blank(),
          legend.text = element_text(size=18), 
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"))

ggsave('b.pdf', width = 12, height = 7)
