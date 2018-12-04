# this is going to make a bunch of plots to help evaluate which
# regions look real and which should be filtered out

# in particular, it should be useful to look at regions overlapping
# and not overlapping genes; the ones overlapping genes should on
# average be more accurate (or at least less often due to poor
# alignment)

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(hexbin)
source('../my_color_palette.R')

args = commandArgs(trailingOnly=TRUE)
tag = args[1]

regions_S288c = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/introgressed_blocks_S288c_', tag, '_quality.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)

regions_C = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/introgressed_blocks_CBS432_', tag, '_quality.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)

regions_N = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/introgressed_blocks_N_45_', tag, '_quality.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)

regions_D = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/introgressed_blocks_DBVPG6304_', tag, '_quality.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)

regions_U = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/introgressed_blocks_UWOPS91_917_1_', tag, '_quality.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)

# no quality file for unknown regions
regions_unk = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/introgressed_blocks_unknown_', tag, '.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)

