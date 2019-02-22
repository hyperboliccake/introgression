library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(dplyr)
library(data.table)
library(VennDiagram)
source('../my_color_palette.R')

tag = 'u3_i.001_tv_l1000_f.01'

strains = c('yjm1252', 'yjm1078', 'yjm248')
# green orange purple
colors = c("#009E2A", "#DA7921", "#8447EB")

## group(comma-separated) count
a = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/',
                     'results/analysis/', tag,
                     '/base_overlap_yjm1252_yjm1078_yjm248.txt', sep=''),
               sep='\t', header=T, stringsAsFactors=F)

area1 = sum(a[which(grepl(strains[1], a$group)),]$count)
area2 = sum(a[which(grepl(strains[2], a$group)),]$count)
area3 = sum(a[which(grepl(strains[3], a$group)),]$count)

n12 = sum(a[which(grepl(strains[1], a$group) &
                  grepl(strains[2], a$group)),]$count)
n23 = sum(a[which(grepl(strains[2], a$group) &
                  grepl(strains[3], a$group)),]$count)
n13 = sum(a[which(grepl(strains[1], a$group) &
                  grepl(strains[3], a$group)),]$count)

n123 = sum(a[which(grepl(strains[1], a$group) &
                   grepl(strains[2], a$group) &
                   grepl(strains[3], a$group)),]$count)

v = draw.triple.venn(area1, area2, area3, n12, n23, n13, n123,
                     category = strains,
                     cat.cex = 3,
                     cex = 3,
                     cat.col = colors,
                     fill = colors,
                     print.mode="percent",
                     alpha = .5)

png(filename = paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/',
                     'results/analysis/', tag,
                     '/plots/base_overlap_3strains.png', sep=''),
    width = 12, height = 12, units = 'in', res=300)
grid.draw(v)
dev.off()

