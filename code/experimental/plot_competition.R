library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)
source('../my_color_palette.R')

a = read.csv('cindy_SUL1_competition.csv', header=T, stringsAsFactors=F)
a$strain = factor(a$strain, levels=c('cerevisiae S288c', 'paradoxus CBS432', 'cerevisiae YJM320', 'none'))


ggplot(a, aes(x=strain, y=slope, shape=as.factor(rep), colour=as.factor(rep))) +
    geom_point(size=3, alpha=.8) + 
    xlab('Version of SUL1 on plasmid') + ylab('Fitness') +
    #guides(fill=FALSE) +
    scale_colour_manual(values=rep(my_color_palette[['mixed']],4)) + 
    labs(shape='Replicate', colour='Replicate')+
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour='gray92'),
          panel.grid.major=element_line(colour='gray92'),
          axis.ticks=element_line(colour="black"),
          axis.line=element_line(),
          axis.title.x = element_text(size=18, margin=margin(t=15)), 
          axis.title.y = element_text(size=18, margin=margin(r=15)), 
          axis.text.x = element_text(size=12,colour="black"), 
          axis.text.y = element_text(colour="black"),
          legend.key = element_rect(fill = "transparent"))

ggplot(a, aes(x=strain, y=slope, colour='x')) +
    geom_point(size=4, alpha=.8) + 
    xlab('Version of SUL1 on plasmid') + ylab('Fitness') +
    #guides(fill=FALSE) +
    scale_colour_manual(values=rep(my_color_palette[['mixed']],4)) + 
    labs(shape='Replicate', colour='Replicate')+
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour='gray92'),
          panel.grid.major=element_line(colour='gray92'),
          axis.ticks=element_line(colour="black"),
          axis.line=element_line(),
          axis.title.x = element_text(size=18, margin=margin(t=15)), 
          axis.title.y = element_text(size=18, margin=margin(r=15)), 
          axis.text.x = element_text(size=12,colour="black"), 
          axis.text.y = element_text(colour="black"),
          legend.key = element_rect(fill = "transparent"))

ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/experimental/cindy_SUL1', '.pdf',sep=''), width = 9, height = 6)

t.test(a[which(a$strain=='cerevisiae S288c'),]$slope, a[which(a$strain=='cerevisiae YJM320'),]$slope)
t.test(a[which(a$strain=='cerevisiae S288c'),]$slope, a[which(a$strain=='cerevisiae YJM320'),]$slope, alternative='less')
#t.test(a[which(a$strain=='cerevisiae S288c'),]$slope, a[which(a$strain=='cerevisiae YJM320'),]$slope, paired=TRUE)
t.test(a[which(a$strain=='cerevisiae S288c'),]$slope, a[which(a$strain=='cerevisiae YJM320'),]$slope, alternative='less', paired=TRUE)
#t.test(a[which(a$strain=='none'),]$slope, a[which(a$strain=='cerevisiae YJM320'),]$slope, alternative='less', paired=TRUE)
