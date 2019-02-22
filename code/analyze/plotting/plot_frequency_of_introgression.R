library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(dplyr)
source('../my_color_palette.R')

args = commandArgs(trailingOnly=TRUE)
tag = args[1]

intd_refs = c("par")

## bases
a = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/',
                     'analysis/', tag, '/introgressed_frequency.txt', sep=''),
               sep = '\t',
               header = TRUE,
               stringsAsFactors = FALSE)

ggplot(a, aes(x=num_strains, y=count/1000, fill='a')) +
    geom_histogram(stat="identity") +
    xlab('Number of strains') +
    ylab('Introgressed bases (kb)') +
    scale_fill_manual(values =c(my_color_palette[['introgressed']]) ) +
    guides(fill=FALSE) +
    scale_y_continuous(expand = c(0,0)) + 
    scale_x_continuous(expand = c(0,0), limits=c(.5,max(a$num_strains)*1.05)) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.line=element_line(),
          axis.ticks=element_line(colour="black"),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(colour="black",size=12), 
          axis.text.y = element_text(colour="black",size=12))
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/frequency_hist_',tag,'.pdf',sep=''), width = 8.4, height = 7.6)

## genes bar chart
ag = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/',
                      'introgression/results/analysis/', tag,
                      '/genes_strain_hist_', tag,
                      '.txt', sep=''), sep='\t', header=T,
                stringsAsFactors=F)
ggplot(ag, aes(x=reorder(gene, -num_strains), y=num_strains, fill='x')) + 
    geom_bar(stat='identity',position='dodge') + 
    xlab('Gene') +
    ylab('Number of strains introgressed') +
    scale_fill_manual(values =c(my_color_palette[['introgressed']]) ) +
    guides(fill=FALSE) +
    scale_y_continuous(expand = c(0,0), limits=c(0,93)) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.ticks=element_line(colour="black"),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(size=.5,angle = 45,vjust = 1,hjust=1,colour="black"), 
          axis.text.y = element_text(colour="black"))
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/num_strains_per_gene_bar_',tag,'.pdf',sep=''), width = 12, height = 6)

## genes histogram
ag2 = ag %>%
    group_by(num_strains) %>%
    summarise(n=n())
ggplot(ag2, aes(x=num_strains, y=n, fill='x')) + 
    geom_bar(stat='identity',position='dodge') + 
    xlab('Number of strains introgressed') +
    ylab('Number of genes') +
    scale_fill_manual(values = c(my_color_palette[['introgressed']]) ) +
    guides(fill=FALSE) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), limits=c(.5,93.5), breaks=c(seq(5,93,5),93)) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.ticks=element_line(colour="black"),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(size=12,colour="black"), 
          axis.text.y = element_text(size=16,colour="black"))
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/num_strains_per_gene_hist_',tag,'.pdf',sep=''), width = 12, height = 6)







