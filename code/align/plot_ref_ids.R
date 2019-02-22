library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(hexbin)
library(stringr)
library(data.table)
library(dplyr)
source('../my_color_palette.R')

args = commandArgs(trailingOnly=TRUE)
refs = c("S288c", "CBS432", "N_45", "DBVPG6304", "UWOPS91_917_1")

rainbow = c("#E13939", "#DA7921", "#E1A939", "#009E2A", "#007CEB", "#8447EB")

## pair id
a = read.table(
    paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/',
          'analysis/ref_ids_',
          paste(refs, collapse='_'), '.txt', sep=''),
    sep='\t', header=T, stringsAsFactors=F)

## pair chromosome id
a_summary = read.table(
    paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/',
          'analysis/ref_ids_summary_',
          paste(refs, collapse='_'), '.txt', sep=''),
    sep='\t', header=T, stringsAsFactors=F)


## =========
## summary divergence info
##=========

k = c(paste(refs[1], refs[2], sep=','),
      paste(refs[1], refs[3], sep=','),
      paste(refs[1], refs[4], sep=','),
      paste(refs[1], refs[5], sep=','),
      paste(refs[2], refs[3], sep=','), 
      paste(refs[4], refs[5], sep=','),
      paste(refs[2], refs[4], sep=','),
      paste(refs[2], refs[5], sep=','),
      paste(refs[3], refs[4], sep=','),
      paste(refs[3], refs[5], sep=',')
      )
cols = c("#E13939", "#E1A939", "#007CEB", "#009E2A", "#DA7921", "#00A89C",
         rep('gray75', 4))

a$pair = factor(a$pair, levels = k)

# violin plot for each pair

ggplot(a, aes(x = pair, fill = pair, y=id)) +
    geom_violin(adjust=3, draw_quantiles = c(0.25, 0.5, 0.75)) + 
    xlab('strain pair') +
    ylab('identity (in 100 bp windows)') + 
    #scale_colour_manual(values=cols) + 
    scale_fill_manual(values=cols) +
    #scale_alpha_manual(values=alphas) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"),
          panel.grid.major=element_line(colour="gray92"),
          legend.position = "none",
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(size=6, colour="black"),
          axis.text.y = element_text(colour="black"))

ggsave(
    paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',
          'plots/', 'ref_ids_violin_quantile_', paste(refs, collapse='_'),'.png', sep=''),
       width = 14, height = 8, dpi=300)

ggplot(a, aes(x = pair, fill = pair, y=id)) +
    geom_violin(adjust=3) + 
    geom_boxplot(width=.1, outlier.colour = NA, position = 'dodge') +
    xlab('strain pair') +
    ylab('identity (in 100 bp windows)') + 
    #scale_colour_manual(values=cols) + 
    scale_fill_manual(values=cols) +
    #scale_alpha_manual(values=alphas) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"),
          panel.grid.major=element_line(colour="gray92"),
          legend.position = "none",
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(size=6, colour="black"),
          axis.text.y = element_text(colour="black"))

ggsave(
    paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',
          'plots/', 'ref_ids_violin_box_', paste(refs, collapse='_'),'.png', sep=''),
       width = 14, height = 8, dpi=300)



# divergence
a$div = 1 - a$id

transparent = '#00000000'

ggplot(a, aes(x = pair, fill = pair, y=div)) +
    geom_violin(adjust=3) + 
    geom_boxplot(width=.1, outlier.colour = NA, position = 'dodge') +
    xlab('strain pair') +
    ylab('divergence (in 100 bp windows)') + 
    scale_fill_manual(values=cols) +
    scale_y_continuous(breaks = seq(0, 1, .1),
                       minor_breaks = seq(0, 1, .01),
                       expand = c(0,0)) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray90"),
          panel.grid.major.y=element_line(colour="gray80"),
          panel.grid.major.x=element_line(colour=transparent),
          legend.position = "none",
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(size=6, colour="black"),
          axis.text.y = element_text(colour="black"))

ggsave(
    paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',
          'plots/', 'ref_ids_violin_div_box_', paste(refs, collapse='_'),'.png', sep=''),
       width = 14, height = 8, dpi=300)

## =========
## histogram and density plot for each pair, overlaid
##=========

## exclude pairs with S288c
a = a[which(!grepl(refs[1], a$pair)),]

k = c(paste(refs[2], refs[3], sep=','), 
      paste(refs[4], refs[5], sep=','),
      paste(refs[2], refs[4], sep=','),
      paste(refs[2], refs[5], sep=','),
      paste(refs[3], refs[4], sep=','),
      paste(refs[3], refs[5], sep=',')
      )
cols = c("#DA7921", "#00A89C", rep('gray75', 4))

a$pair = factor(a$pair, levels = k)

alphas = c(.4, .4, rep(.2, 4))

ggplot(a, aes(x = id, fill = pair, colour = pair, alpha=pair)) +
    geom_density(adjust=3) + 
    xlab('pairwise identity (window size = 100bp)') +
    ylab('frequency') + 
    scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
    scale_y_continuous(limits=c(0,37), expand=c(0,0)) +
    scale_colour_manual(values=cols) + 
    scale_fill_manual(values=cols) +
    scale_alpha_manual(values=alphas) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"),
          panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"))

ggsave(
    paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',
          'plots/', 'ref_ids_density_', paste(refs, collapse='_'),'.png', sep=''),
       width = 10, height = 8, dpi=300)

ggplot(a, aes(x = id, fill = pair, colour = pair, alpha=pair)) +
    geom_density(adjust=3) + 
    xlab('pairwise identity (window size = 100bp)') +
    ylab('frequency') + 
    scale_x_continuous(limits=c(.85,1), expand=c(0,0)) +
    scale_y_continuous(limits=c(0,37), expand=c(0,0)) +
    scale_colour_manual(values=cols) + 
    scale_fill_manual(values=cols) +
    scale_alpha_manual(values=alphas) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"),
          panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"))

ggsave(
    paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',
          'plots/', 'ref_ids_density_zoomed_', paste(refs, collapse='_'),'.png', sep=''),
       width = 10, height = 8, dpi=300)

cols = c("#DA7921FF", "#00A89CFF", rep('#A1A1A115', 4))

ggplot(a, aes(x = id, fill = pair, colour = pair, alpha = pair)) +
    geom_histogram(binwidth = .01, position = 'identity') + 
    xlab('pairwise identity (window size = 100bp)') +
    ylab('count') + 
    #scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
    #scale_y_continuous(limits=c(0,7000), expand=c(0,0)) +
    scale_colour_manual(values=cols) +
    scale_fill_manual(values=cols) + 
    scale_alpha_manual(values=alphas) + 
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"),
          panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"))

ggsave(
    paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',
          'plots/', 'ref_ids_histogram_', paste(refs, collapse='_'),'.png', sep=''),
       width = 10, height = 8, dpi=300)

b = a %>%
    group_by(pair) %>%
    summarize(pair_total = n())
d = a %>%
    group_by(.dots=c("pair", "id")) %>%
    summarize(id_count = n())
g = merge(b, d)
g$freq = g$id_count / g$pair_total

ggplot(g, aes(x = id, y = freq, fill = pair, colour = pair, alpha = pair)) +
    geom_bar(stat='identity', position='identity') + 
    xlab('pairwise identity (window size = 100bp)') +
    ylab('count') + 
    #scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
    #scale_y_continuous(limits=c(0,7000), expand=c(0,0)) +
    scale_colour_manual(values=cols) +
    scale_fill_manual(values=cols) + 
    scale_alpha_manual(values=alphas) + 
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"),
          panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"))

ggsave(
    paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',
          'plots/', 'ref_ids_histogram_freq_', paste(refs, collapse='_'),'.png', sep=''),
       width = 10, height = 8, dpi=300)


