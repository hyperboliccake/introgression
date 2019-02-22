library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(dplyr)
library(data.table)
library(VennDiagram)
source('../my_color_palette.R')

args = commandArgs(trailingOnly=TRUE)
tag = args[1]
refs = c("CBS432", "DBVPG6304", "N_45", "UWOPS91_917_1")

rainbow = c("#E13939", "#DA7921", "#E1A939", "#009E2A", "#007CEB", "#8447EB")
ref_colors = c(rainbow[1], rainbow[3], rainbow[5], rainbow[4], 'gray70') # east west


## strain label count
## label is comma-separated, values are reference strain names, any
## (indicating any reference strain(s)), and other (indicating other
## set of predictions)
a = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/',
                     'results/analysis/', tag, '/state_counts_comparison.txt', sep=''),
               sep='\t', header=T, stringsAsFactors=F)

## venn diagram with two circles: other and everything not inluding
## other (overlap is other,any)

#b = a[which(!grepl('other', a$label)),] %>%
#    group_by(label) %>%
#    summarize(total = sum(count))
x_r = sum(a[which(!grepl('other', a$label)),]$count)
x_o = sum(a[which(a$label=='other'),]$count)
x_or = sum(a[which(a$label=='other,any'),]$count)
x_r = x_r + x_or
x_o = x_o + x_or

v = draw.pairwise.venn(x_r, x_o, x_or, category = c('four references', 'one reference'))

png(filename = paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/',
                     'results/analysis/', tag,
                     '/plots/compare_predictions_venn_general.png', sep=''),
    width = 12, height = 12, units = 'in', res=300)
grid.draw(v)
dev.off()


## venn diagram with five circles: each of four references plus other
## (ignore entries with any)

b = a[which(!grepl('any', a$label)),]


area1 = sum(b[which(grepl(refs[1], b$label)),]$count)
area2 = sum(b[which(grepl(refs[2], b$label)),]$count)
area3 = sum(b[which(grepl(refs[3], b$label)),]$count)
area4 = sum(b[which(grepl(refs[4], b$label)),]$count)
area5 = sum(b[which(grepl('other', b$label)),]$count)

n12 = sum(b[which(grepl(refs[1], b$label) &
                  grepl(refs[2], b$label)),]$count)
n13 = sum(b[which(grepl(refs[1], b$label) &
                  grepl(refs[3], b$label)),]$count)
n14 = sum(b[which(grepl(refs[1], b$label) &
                  grepl(refs[4], b$label)),]$count)
n15 = sum(b[which(grepl(refs[1], b$label) &
                  grepl('other', b$label)),]$count)

n23 = sum(b[which(grepl(refs[2], b$label) &
                  grepl(refs[3], b$label)),]$count)
n24 = sum(b[which(grepl(refs[2], b$label) &
                  grepl(refs[4], b$label)),]$count)
n25 = sum(b[which(grepl(refs[2], b$label) &
                  grepl('other', b$label)),]$count)

n34 = sum(b[which(grepl(refs[3], b$label) &
                  grepl(refs[4], b$label)),]$count)
n35 = sum(b[which(grepl(refs[3], b$label) &
                  grepl('other', b$label)),]$count)

n45 = sum(b[which(grepl(refs[4], b$label) &
                  grepl('other', b$label)),]$count)

n123 = sum(b[which(grepl(refs[1], b$label) &
                   grepl(refs[2], b$label) &
                   grepl(refs[3], b$label)),]$count)
n124 = sum(b[which(grepl(refs[1], b$label) &
                   grepl(refs[2], b$label) &
                   grepl(refs[4], b$label)),]$count)
n125 = sum(b[which(grepl(refs[1], b$label) &
                   grepl(refs[2], b$label) &
                   grepl('other', b$label)),]$count)

n134 = sum(b[which(grepl(refs[1], b$label) &
                   grepl(refs[3], b$label) &
                   grepl(refs[4], b$label)),]$count)
n135 = sum(b[which(grepl(refs[1], b$label) &
                   grepl(refs[3], b$label) &
                   grepl('other', b$label)),]$count)

n145 = sum(b[which(grepl(refs[1], b$label) &
                   grepl(refs[4], b$label) &
                   grepl('other', b$label)),]$count)

n234 = sum(b[which(grepl(refs[2], b$label) &
                   grepl(refs[3], b$label) &
                   grepl(refs[4], b$label)),]$count)
n235 = sum(b[which(grepl(refs[2], b$label) &
                   grepl(refs[3], b$label) &
                   grepl('other', b$label)),]$count)

n245 = sum(b[which(grepl(refs[2], b$label) &
                   grepl(refs[4], b$label) &
                   grepl('other', b$label)),]$count)

n345 = sum(b[which(grepl(refs[3], b$label) &
                   grepl(refs[4], b$label) &
                   grepl('other', b$label)),]$count)

n1234 = sum(b[which(grepl(refs[1], b$label) &
                    grepl(refs[2], b$label) &
                    grepl(refs[3], b$label) &
                    grepl(refs[4], b$label)),]$count)
n1235 = sum(b[which(grepl(refs[1], b$label) &
                    grepl(refs[2], b$label) &
                    grepl(refs[3], b$label) &
                    grepl('other', b$label)),]$count)
n1245 = sum(b[which(grepl(refs[1], b$label) &
                    grepl(refs[2], b$label) &
                    grepl(refs[4], b$label) &
                    grepl('other', b$label)),]$count)
n1345 = sum(b[which(grepl(refs[1], b$label) &
                    grepl(refs[3], b$label) &
                    grepl(refs[4], b$label) &
                    grepl('other', b$label)),]$count)

n2345 = sum(b[which(grepl(refs[2], b$label) &
                    grepl(refs[3], b$label) &
                    grepl(refs[4], b$label) &
                    grepl('other', b$label)),]$count)

n12345 = sum(b[which(grepl(refs[1], b$label) &
                     grepl(refs[2], b$label) &
                     grepl(refs[3], b$label) &
                     grepl(refs[4], b$label) &
                     grepl('other', b$label)),]$count)

v = draw.quintuple.venn(area1, area2, area3, area4, area5,
                        n12, n13, n14, n15, n23, n24, n25, n34, n35, n45,
                        n123, n124, n125, n134, n135, n145, n234, n235, n245, n345,
                        n1234, n1235, n1245, n1345, n2345,
                        n12345,
                        fill=ref_colors,
                        cat.col=ref_colors,
                        category = c(refs, 'one paradoxus reference'))

png(filename = paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/',
                     'results/analysis/', tag,
                     '/plots/compare_predictions_venn_specific.png', sep=''),
    width = 12, height = 12, units = 'in', res=300)
par(mar = c(4, 4, 4, 4))
grid.draw(v)
dev.off()
