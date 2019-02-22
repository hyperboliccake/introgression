# pretty much the same as plot_genome, but instead of separate rows
# for all strains, plots them on top of each other; darker regions
# have more introgression
# also, plot all chromosomes at once

library(ggplot2)
library(viridis)
require(grDevices)
source('../my_color_palette.R')

unknown = FALSE

chrms = c('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI')
num_chrms = length(chrms)

# centromere and telomere locations for each chromosome 
cen_starts = c(151465, 238207, 114385, 449711, 151987, 148510, 496920, 105586, 355629, 436307, 440129, 150828, 268031, 628758, 326584, 555957) - 1

cen_ends = c(151582,238323,114501,449821,152104,148627,497038,105703,355745,436425,440246,150947,268149,628875,326702,556073) - 1

tel_coords = c(1,801,229411,230218,1,6608,812379,813184,1,1098,315783,316620,1,904,1524625,1531933,1,6473,569599,576874,1,5530,269731,270161,1,781,1083635,1090940,1,5505,556105,562643,1,7784,439068,439888,1,7767,744902,745751,1,807,665904,666816,1,12085,1064281,1078177,1,6344,923541,924431,1,7428,783278,784333,1,847,1083922,1091291,1,7223,942396,948010) - 1

tel_left_starts = tel_coords[seq(1, 64, 4)]
tel_left_ends = tel_coords[seq(2, 64, 4)]
tel_right_starts = tel_coords[seq(3, 64, 4)]
tel_right_ends = tel_coords[seq(4, 64, 4)]

# colors
par_color_dark = my_color_palette[['introgressed']]
unk_color_dark = 'gray50'
r=col2rgb(par_color_dark)/255
par_color = rgb(r[1],r[2],r[3],alpha=1)
r=col2rgb(unk_color_dark)/255
unk_color = rgb(r[1],r[2],r[3],alpha=.1)

# create plots directory for this tag (might already exist)
args = commandArgs(trailingOnly=TRUE)
tag = args[1]
suffix = ''
if (length(args) == 2) {
    suffix = args[2]
}

dir.create(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/',sep=''),showWarnings=FALSE)

# read in annotated and filtered regions file (paradoxus/introgressed)
regions = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/introgressed_blocks', suffix, '_par_', tag, '_summary_plus.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)

# read in regions file (cerevisiae/not introgressed)
regions_cer = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/introgressed_blocks_cer_', tag, '.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)

# read in regions file (unknown, if it exists)
regions_unk = data.frame(matrix(ncol=length(names(regions_cer)), nrow=0))
names(regions_unk) = names(regions_cer)
if (unknown){
regions_unk = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/introgressed_blocks_unknown_', tag, '.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)
}

# all strains, in arbitrary order
strains_list = unique(regions$strain)
strains = data.frame(strain=strains_list, index=1:length(strains_list))

# row/spacing parameters
row_height = 2 # chromosome height
padding = 1.5 # between chromosomes
vmargin = 2 # for top and bottom

# plot png file
png(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/plots/all_chrms', suffix, '_', tag, '.png', sep=''), 6400, 3600)

# set margins: bottom, left, top, right
par(mar=c(30, 20, 4, 2))

# coordinates for longest chromosome
seq_start = 0
seq_end = max(tel_right_ends)

# base plot
# type n doesn't produce any points or lines
# (left, right), (bottom, top)
# xaxs and yaxs args to specify axes go to edge of plot region
# xaxt and yaxt args to specify...?
total_height = 2 * vmargin + num_chrms * (row_height + padding) - padding
plot(c(seq_start, seq_end+1),
     c(0, total_height),
     type = "n",
     xlab = "Position (Mb)", cex.lab=13, mgp = c(22, 0, 0),
     ylab = "", main = "", xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=F)

# move x axis label and title closer to axis
# title(xlab = paste("chromosome ", chrm, " ", m, sep = ''), line = 3, cex.lab=1.7)

for (ci in 1:length(chrms))
{
    chrm = chrms[ci]
    print(chrm)

    regions_chrm = regions[which(regions$chromosome == chrm),]
    regions_cer_chrm = regions_cer[which(regions_cer$chromosome == chrm),]
    regions_unk_chrm = regions_unk[which(regions_unk$chromosome == chrm),]

    row_bottom = vmargin + (ci - 1) * (row_height + padding)
    row_top = row_bottom + row_height
    row_middle = (row_bottom + row_top) / 2
    

    # plot par regions
    regions_chrm = merge(strains, regions_chrm, by = 'strain', all.x = T)
    for (r in 1:nrow(regions_chrm)) {
        region_start = regions_chrm[r,]$start
        region_end = regions_chrm[r,]$end
        rect(region_start, 
             row_bottom,
             region_end, 
             row_top,
             col = par_color, border=par_color)
    }

    # plot chromosome
    rect(seq_start,
         row_bottom,
         tel_right_ends[ci],
         row_top,
         col = NULL, border = 'black')
    feature_height = 0
    feature_color = 'black'
    # plot centromere box and point
    rect(cen_starts[ci],
         row_bottom-feature_height,
         cen_ends[ci],
         row_top+feature_height,
         col = feature_color, border = feature_color)
    points((cen_starts[ci]+cen_ends[ci])/2, row_middle, cex=5,col=feature_color)
    # plot telomere boxes and point
    rect(tel_left_starts[ci],
         row_bottom-feature_height,
         tel_left_ends[ci],
         row_top+feature_height,
         col = feature_color, border = feature_color)
    #points(tel_left_ends[ci], row_middle, pch=19)
    rect(tel_right_starts[ci],
         row_bottom-feature_height,
         tel_right_ends[ci],
         row_top+feature_height,
         col = feature_color, border = feature_color)
    #points(tel_right_starts[ci], row_middle, pch=19)

}

# plot chromosome labels
positions = seq(vmargin + row_height/2,
                num_chrms*(padding + row_height) - padding + vmargin,
                padding+row_height)
axis(2, at=positions, labels=chrms, las=1, cex.axis=10, tick=FALSE, line=-.5, mgp=c(0,5,0))

# plot position labels
positions = seq(seq_start, seq_end+1, 100000)
axis(1, at=positions, labels=positions/1000000, las=1, cex.axis=10, line = 0, mgp = c(15,8,0))

dev.off()

