library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(hexbin)
library(data.table)
library(VennDiagram)
source('../my_color_palette.R')

args = commandArgs(trailingOnly=TRUE)
tag = args[1]
#refs = c("S288c", "CBS432", "N_45", "DBVPG6304", "UWOPS91_917_1")
refs = c("CBS432", "DBVPG6304", "N_45", "UWOPS91_917_1")

a = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/state_counts_by_strain.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)

##=====
# venn diagram showing introgression from all the different sources
##=====

## filtering 2i

a1 = sum(a[,which((names(a) %like% 'num_bases') & (names(a) %like% 'CBS432') &
                  (names(a) %like% 'filtered2') & !(names(a) %like% 'inclusive'))])
a2 = sum(a[,which((names(a) %like% 'num_bases') & (names(a) %like% 'DBVPG6304') &
                  (names(a) %like% 'filtered2') & !(names(a) %like% 'inclusive'))])
a3 = sum(a[,which((names(a) %like% 'num_bases') & (names(a) %like% 'N_45') &
                  (names(a) %like% 'filtered2') & !(names(a) %like% 'inclusive'))])
a4 = sum(a[,which((names(a) %like% 'num_bases') & (names(a) %like% 'UWOPS91_917_1') &
                  (names(a) %like% 'filtered2') & !(names(a) %like% 'inclusive'))])

a12 = sum(a[,which((names(a) %like% 'num_bases') & (names(a) %like% 'CBS432') &
                   (names(a) %like% 'DBVPG6304') &
                   (names(a) %like% 'filtered2') & !(names(a) %like% 'inclusive'))])
a13 = sum(a[,which((names(a) %like% 'num_bases') & (names(a) %like% 'CBS432') &
                   (names(a) %like% 'N_45') &
                   (names(a) %like% 'filtered2') & !(names(a) %like% 'inclusive'))])
a14 = sum(a[,which((names(a) %like% 'num_bases') & (names(a) %like% 'CBS432') &
                   (names(a) %like% 'UWOPS91_917_1') &
                   (names(a) %like% 'filtered2') & !(names(a) %like% 'inclusive'))])
a23 = sum(a[,which((names(a) %like% 'num_bases') & 
                   (names(a) %like% 'DBVPG6304') & (names(a) %like% 'N_45') &
                   (names(a) %like% 'filtered2') & !(names(a) %like% 'inclusive'))])
a24 = sum(a[,which((names(a) %like% 'num_bases') & 
                   (names(a) %like% 'DBVPG6304') & (names(a) %like% 'UWOPS91_917_1') &
                   (names(a) %like% 'filtered2') & !(names(a) %like% 'inclusive'))])
a34 = sum(a[,which((names(a) %like% 'num_bases') & 
                   (names(a) %like% 'N_45') & (names(a) %like% 'UWOPS91_917_1') &
                   (names(a) %like% 'filtered2') & !(names(a) %like% 'inclusive'))])

a123 = sum(a[,which((names(a) %like% 'num_bases') & (names(a) %like% 'CBS432') &
                    (names(a) %like% 'DBVPG6304') & (names(a) %like% 'N_45') &
                    (names(a) %like% 'filtered2') & !(names(a) %like% 'inclusive'))])

a124 = sum(a[,which((names(a) %like% 'num_bases') & (names(a) %like% 'CBS432') &
                    (names(a) %like% 'DBVPG6304') & (names(a) %like% 'UWOPS91_917_1') &
                    (names(a) %like% 'filtered2') & !(names(a) %like% 'inclusive'))])
a134 = sum(a[,which((names(a) %like% 'num_bases') & (names(a) %like% 'CBS432') &
                    (names(a) %like% 'N_45') & (names(a) %like% 'UWOPS91_917_1') &
                    (names(a) %like% 'filtered2') & !(names(a) %like% 'inclusive'))])
a234 = sum(a[,which((names(a) %like% 'num_bases') & (names(a) %like% 'N_45') &
                    (names(a) %like% 'DBVPG6304') & (names(a) %like% 'UWOPS91_917_1') &
                    (names(a) %like% 'filtered2') & !(names(a) %like% 'inclusive'))])

a1234 = sum(a$num_bases_CBS432_or_DBVPG6304_or_N_45_or_UWOPS91_917_1_filtered2i)

venn.plot = draw.quad.venn(a1, a2, a3, a4, a12, a13, a14, a23, a24, a34, a123, a124, a134, a234, a1234,
                           category = refs,
                           fill=c("#E13939", "#007CEB", "#E1A939", "#009E2A"),
                           cat.col=c("#E13939", "#007CEB", "#E1A939", "#009E2A")
                           )

png(filename = paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/plots/venn.png', sep=''))
grid.draw(venn.plot)
dev.off()
