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
tag = args[1]
refs = c("CBS432", "N_45", "DBVPG6304", "UWOPS91_917_1")

rainbow = c("#E13939", "#DA7921", "#E1A939", "#009E2A", "#007CEB", "#8447EB")

a = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/',
                     'analysis/', tag, '/filter_2_thresholds_', tag, '.txt', sep=''),
               sep='\t', header=T, stringsAsFactors=F)

## threshold predicted_state alternative_states count

## =========
## plot number assigned to 1, 2, 3, 4 states vs threshold, where
## threshold is the fraction of identity with the assigned reference
## that another reference needs to meet in order to be considered an
## alternative---so if threshold is very high, other references will
## meet threshold less often, and we'll assign things to fewer
## references
##=========

a$num_alts = str_count(a$alternative_states, ',') + 1

b = a %>%
    group_by(.dots=c('num_alts', 'threshold')) %>%
    summarize(count = sum(count))
b_totals = b %>%
    group_by(threshold) %>%
    summarize(total=sum(count))
b = merge(b, b_totals)
b$frac = b$count / b$total

ggplot(b, aes(x = threshold, y = frac,
              colour = as.factor(num_alts))) +
    geom_point() + 
    geom_line() +
    ylab('fraction assigned') +
    xlab('threshold') + 
    scale_colour_manual(values=rainbow) + 
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
          tag, '/plots/', 'filter_thresholds_1.png', sep=''),
       width = 10, height = 8, dpi=300)

## =========
## plot number assigned to each reference state only vs threshold,
## as well as the interesting pairs
##=========

k = c(refs[1], paste(refs[1], refs[2], sep=','), refs[2],
      refs[3], paste(refs[3], refs[4], sep=','), refs[4])

b = a[which(a$alternative_states %in% k),]
b = b %>%
    group_by(.dots = c('threshold', 'alternative_states')) %>%
    summarize(count = sum(count))

b_all = b %>%
    group_by(threshold) %>%
    summarize(count=sum(count))
b_all$alternative_states = 'all'
b = bind_rows(b, b_all)

b_totals = a %>%
    group_by(threshold) %>%
    summarize(total=sum(count))
b = merge(b, b_totals)

b$frac = b$count / b$total

k = c(k, 'all')
b$alternative_states = factor(b$alternative_states, levels = k)

cols = c("#E13939", "#DA7921", "#E1A939", "#007CEB", "#00A89C", "#009E2A", 'gray75')

ggplot(b, aes(x = threshold, y = frac,
              colour = as.factor(alternative_states))) +
    geom_line() +
    ylab('fraction assigned') +
    xlab('threshold') + 
    scale_colour_manual(values=cols) + 
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
          tag, '/plots/', 'filter_thresholds_2.png', sep=''),
       width = 10, height = 8, dpi=300)


# exclude doubles

k = c(refs[1], refs[2],
      refs[3], refs[4])

b = a[which(a$alternative_states %in% k),]
b = b %>%
    group_by(.dots = c('threshold', 'alternative_states')) %>%
    summarize(count = sum(count))

b_all = b %>%
    group_by(threshold) %>%
    summarize(count=sum(count))
b_all$alternative_states = 'all'
b = bind_rows(b, b_all)

b_totals = a %>%
    group_by(threshold) %>%
    summarize(total=sum(count))
b = merge(b, b_totals)

b$frac = b$count / b$total

k = c(k, 'all')
b$alternative_states = factor(b$alternative_states, levels = k)

cols = c("#E13939", "#E1A939", "#007CEB", "#009E2A", 'gray75')

ggplot(b, aes(x = threshold, y = frac,
              colour = as.factor(alternative_states))) +
    geom_line() +
    ylab('fraction assigned') +
    xlab('threshold') + 
    scale_colour_manual(values=cols) + 
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
          tag, '/plots/', 'filter_thresholds_3.png', sep=''),
       width = 10, height = 8, dpi=300)


# stacked line plot

k = c(refs[1], paste(refs[1], refs[2], sep=','), refs[2],
      refs[3], paste(refs[3], refs[4], sep=','), refs[4])

b = a[which(a$alternative_states %in% k),]
b = b %>%
    group_by(.dots = c('threshold', 'alternative_states')) %>%
    summarize(count = sum(count))

b_totals = a %>%
    group_by(threshold) %>%
    summarize(overall_total=sum(count))
total = b_totals$overall_total[1]

b_remaining = b %>%
    group_by(threshold) %>%
    summarize(count=total-sum(count))
b_remaining$alternative_states = 'other'
b = bind_rows(b, b_remaining)

b$frac = b$count / total


k = c(k, 'other')
b$alternative_states = factor(b$alternative_states, levels = rev(k))

cols = rev(c("#E13939", "#DA7921", "#E1A939", "#007CEB", "#00A89C", "#009E2A", 'gray75'))

ggplot(b, aes(x = threshold, y = frac,
              fill = as.factor(alternative_states))) +
    geom_area(position='stack') +
    ylab('fraction assigned') +
    xlab('threshold') + 
    scale_x_continuous(limits=c(.5,1), expand=c(0,0)) +
    scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
    scale_fill_manual(values=cols) + 
    theme(panel.background=element_rect(fill=NA),
          #panel.ontop = TRUE,
          legend.title = element_blank(),
          panel.grid.minor=element_line(colour="gray92"),
          panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"))

ggsave(
    paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',
          tag, '/plots/', 'filter_thresholds_4.png', sep=''),
       width = 10, height = 8, dpi=300)
