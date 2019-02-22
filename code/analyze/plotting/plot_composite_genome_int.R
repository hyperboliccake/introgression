# for each strain, make one plot showing all chromosomes, with
# population blocks and introgression blocks shownlibrary(ggplot2)

library(viridis)
library(RColorBrewer)
require(grDevices)

###########
## set up stuff: references, colors, analysis tag, plot parameters
###########

args = commandArgs(trailingOnly=TRUE)
tag = args[1]
structure_id = "27"
num_pops = 6

refs = c("CBS432", "N_45", "DBVPG6304", "UWOPS91_917_1") 

chrms = c('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII',
          'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI')
num_chrms = length(chrms)

# centromere and telomere locations for each chromosome 
cen_starts = c(151465, 238207, 114385, 449711,
               151987, 148510, 496920, 105586,
               355629, 436307, 440129, 150828,
               268031, 628758, 326584, 555957) - 1

cen_ends = c(151582, 238323, 114501, 449821,
             152104, 148627, 497038, 105703,
             355745, 436425, 440246, 150947,
             268149, 628875, 326702, 556073) - 1

tel_coords = c(1,801,229411,230218,
               1,6608,812379,813184,
               1,1098,315783,316620,
               1,904,1524625,1531933,
               1,6473,569599,576874,
               1,5530,269731,270161,
               1,781,1083635,1090940,
               1,5505,556105,562643,
               1,7784,439068,439888,
               1,7767,744902,745751,
               1,807,665904,666816,
               1,12085,1064281,1078177,
               1,6344,923541,924431,
               1,7428,783278,784333,
               1,847,1083922,1091291,
               1,7223,942396,948010) - 1

tel_left_starts = tel_coords[seq(1, 64, 4)]
tel_left_ends = tel_coords[seq(2, 64, 4)]
tel_right_starts = tel_coords[seq(3, 64, 4)]
tel_right_ends = tel_coords[seq(4, 64, 4)]

# colors
rainbow = c("#E13939", "#DA7921", "#E1A939", "#009E2A", "#007CEB", "#8447EB")
pop_cols = list()
for (i in 1:num_pops)
    pop_cols[[paste(i)]] = rainbow[i]
# refs east(1&2) west(3&4) other
int_cols = list()
int_cols[[refs[1]]] = rainbow[1]
int_cols[[refs[2]]] = rainbow[3]
int_cols[[refs[3]]] = rainbow[5]
int_cols[[refs[4]]] = rainbow[4]
int_cols[[paste(refs[1], refs[2], sep=',')]] = rainbow[2]
int_cols[[paste(refs[3], refs[4], sep=',')]] = "#00A89C"
other_col = "gray50"

# other plot dimensional parameters
row_height = 2 # chromosome height
padding = 3 # between chromosomes
vmargin = 3 # for top and bottom


###########
## read in all introgressed blocks
###########

blocks = list()
strains = c()
for (ref in refs) {
    
    a = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/',
                         'results/analysis/', tag, '/blocks_', ref,
                         '_', tag, '_filtered2intermediate.txt', sep=''),
                   sep='\t', header=T, stringsAsFactors=F)
    blocks[[ref]] = a
    strains = c(strains, a$strain)
}
strains = unique(strains)

###########
## loop through all strains, adding to same plot for each
###########

## plot png file
png(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/',
          'results/analysis/', tag, '/plots/genome_composite_int_',
          tag, '.png', sep=''), 6400, 3600)

## set margins: bottom, left, top, right
par(mar=c(5, 10, 4, 2))

## coordinates for longest chromosome
seq_start = 0
seq_end = max(tel_right_ends)

## base plot
## type n doesn't produce any points or lines
## (left, right), (bottom, top)
## xaxs and yaxs args to specify axes go to edge of plot region
## xaxt and yaxt args to specify...?
total_height = 2 * vmargin + num_chrms * (row_height + padding) - padding
plot(c(seq_start, seq_end+1),
     c(0, total_height),
     type = "n", xlab = "",
     ylab = "", main = "", xaxt='n', yaxt='n', xaxs='i', yaxs='i',
     mgp=c(2,2,.5), axes=F)

for (ci in 1:length(chrms)) {

    chrm = chrms[ci]
            
    print(chrm) 

    nrows = length(refs)
    row_bottoms = c(vmargin + (ci - 1) * (row_height + padding))
    row_tops = c(row_bottoms[1] + row_height/nrows)
    for (i in 2:nrows) {
        row_bottoms = c(row_bottoms, row_tops[i - 1])
        row_tops = c(row_tops, row_tops[i - 1] + row_height/nrows)
    }
    
    for (strain in strains) {
        
        for (ri in 1:length(refs)) {
            ref = refs[ri]
            
            regions = blocks[[ref]]
        
            regions_chrm = regions[which(regions$chromosome == chrm &
                                         regions$strain == strain),]

            if (nrow(regions_chrm) != 0) {
            
                ## plot introgressed regions (top half of chromosome)
                for (r in 1:nrow(regions_chrm)) {
                    region_start = regions_chrm[r,]$start
                    region_end = regions_chrm[r,]$end

                    col = other_col
                    if (regions_chrm$alternative_states[r] %in% names(int_cols)) {
                        col = int_cols[[regions_chrm$alternative_states[r]]]
                    }
                    rect(region_start,
                         row_bottoms[ri],
                         region_end, 
                         row_tops[ri],
                         col = col, border = col)
                }
            }
        }

        row_bottom = row_bottoms[1]
        row_top = row_tops[nrows]
        row_middle = (row_top + row_bottom) / 2
        
        ## plot chromosome
        rect(seq_start,
             row_bottom,
             tel_right_ends[ci],
             row_top,
             col = NULL, border = 'black')
        feature_height = 0
        feature_color= 'black'
        ## plot centromere box and point
        rect(cen_starts[ci],
             row_bottom-feature_height,
             cen_ends[ci],
             row_top+feature_height,
             col = 'black', border = 'black')
        points((cen_starts[ci]+cen_ends[ci])/2, row_middle, cex=5,col=feature_color)
        ## plot telomere boxes and point
        rect(tel_left_starts[ci],
             row_bottom-feature_height,
             tel_left_ends[ci],
             row_top+feature_height,
             col = 'black', border = 'black')
        ##points(tel_left_ends[ci], row_middle, pch=19)
        rect(tel_right_starts[ci],
             row_bottom-feature_height,
             tel_right_ends[ci],
             row_top+feature_height,
             col = 'black', border = 'black')
        ##points(tel_right_starts[ci], row_middle, pch=19)

    }
    
    
}
## plot chromosome labels
positions = seq(vmargin + row_height/2,
                num_chrms*(padding + row_height) - padding + vmargin,
                padding+row_height)
axis(2, at=positions, labels=chrms, las=1, cex.axis=6, tick=FALSE, line=-.5)

## plot position labels
positions = seq(seq_start, seq_end+1, 100000)
axis(1, at=positions, las=1, mgp=c(3,4,0), cex.axis=6, line = 0)

dev.off()
    



