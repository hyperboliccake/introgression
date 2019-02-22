library(ggplot2)
library(viridis)
require(grDevices)

uncolor_mode = FALSE

unknown = TRUE

chrms = c('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI')

cen_starts = c(151465, 238207, 114385, 449711, 151987, 148510, 496920, 105586, 355629, 436307, 440129, 150828, 268031, 628758, 326584, 555957) - 1


cen_ends = c(151582,238323,114501,449821,152104,148627,497038,105703,355745,436425,440246,150947,268149,628875,326702,556073) - 1

tel_coords = c(1,801,229411,230218,1,6608,812379,813184,1,1098,315783,316620,1,904,1524625,1531933,1,6473,569599,576874,1,5530,269731,270161,1,781,1083635,1090940,1,5505,556105,562643,1,7784,439068,439888,1,7767,744902,745751,1,807,665904,666816,1,12085,1064281,1078177,1,6344,923541,924431,1,7428,783278,784333,1,847,1083922,1091291,1,7223,942396,948010) - 1

tel_left_starts = tel_coords[seq(1, 64, 4)]
tel_left_ends = tel_coords[seq(2, 64, 4)]
tel_right_starts = tel_coords[seq(3, 64, 4)]
tel_right_ends = tel_coords[seq(4, 64, 4)]

vcolors = viridis(7, option = 'plasma')
vcolors100 = viridis(100, option='plasma')

if (uncolor_mode) {
   vcolors = rep('black', 7)
   vcolors100 = rep('black', 100)
}

chrm_color = viridis(5)[2] # at bottom
unk_color = viridis(5)[3]

######

vcolors = viridis(7, option = 'inferno')
chrm_color = 'gray70' # at bottom
par_color = vcolors[5] # orange
unk_color = vcolors[2] # purple

######

args = commandArgs(trailingOnly=TRUE)
tag = args[1]
suffix = ''
if (length(args) == 2) {
    suffix = args[2]
}
dir.create(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/',sep=''),showWarnings=FALSE)


regions = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/introgressed_blocks', suffix, '_par_', tag, '_summary_plus.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)

regions_cer = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/introgressed_blocks_cer_', tag, '.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)

regions_unk = data.frame(matrix(ncol=length(names(regions_cer)), nrow=0))
names(regions_unk) = names(regions_cer)
if (unknown){
regions_unk = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/introgressed_blocks_unknown_', tag, '.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)
}

a_strain = read.table('strains.txt', header=F, stringsAsFactors=F)
strains_list = as.character(a_strain[1,3:ncol(a_strain)])
strains = data.frame(strain=strains_list, index=1:length(strains_list))

row_height = 2
padding = 3
vmargin = 1 # for top and bottom

m=NULL
if (unknown) {
    m = max(regions$number_match_only_par) + 1
} else {
    m = max(regions$number_match_par_not_cer) + 1 }


m = 5000

for (ci in 1:length(chrms)) 
{
    chrm = chrms[ci]
    print(chrm)

    regions_chrm = regions[which(regions$chromosome == chrm),]
    regions_cer_chrm = regions_cer[which(regions_cer$chromosome == chrm),]
    regions_unk_chrm = regions_unk[which(regions_unk$chromosome == chrm),]

    png(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/plots/chrm', chrm, suffix, '_', tag, '.png', sep=''), 6400, 3600)

    # set margins: bottom, left, top, right
    par(mar=c(5, 10, 4, 2))

    seq_start = 0
    seq_end = tel_right_ends[ci]

    # base plot
    # type n doesn't produce any points or lines
    # (left, right), (bottom, top)
    # xaxs and yaxs args to specify axes go to edge of plot region
    # xaxt and yaxt args to specify...?
    total_height = vmargin + row_height + padding + nrow(strains) * (row_height + padding) + vmargin
    plot(c(seq_start, seq_end+1), 
    	 c(0, total_height),
         type = "n", xlab = "",
         ylab = "", main = "", xaxt='n', yaxt='n', xaxs='i', yaxs='i', mgp=c(2,2,.5), axes=F)

    # move x axis label and title closer to axis
    title(xlab = paste("chromosome ", chrm, " ", m, sep = ''), line = 3, cex.lab=1.7)

    # plot chromosome
    rect(seq_start,
         padding,
         seq_end,
         padding + row_height,
         col = chrm_color, border = chrm_color)
    # plot centromere
    rect(cen_starts[ci],
         padding,
         cen_ends[ci],
         padding + row_height,
         col = 'black', border = 'black')
    # plot telomeres
    rect(tel_left_starts[ci],
         padding,
         tel_left_ends[ci],
         padding + row_height,
         col = 'black', border = 'black')
    rect(tel_right_starts[ci],
         padding,
         tel_right_ends[ci],
         padding + row_height,
         col = 'black', border = 'black')


    # plot gray bar for all strains
    for (strain_index in 1:nrow(strains))
    {
        rect(seq_start, 
             (strain_index - 1) * row_height + strain_index * padding + vmargin + row_height, 
             seq_end, 
             (strain_index - 1) * row_height + strain_index * padding + row_height + vmargin + row_height, 
             col = 'gray75', border=F)

    }

    # plot cer regions in white on top of that
    regions_cer_chrm = merge(strains, regions_cer_chrm, by = 'strain', all.x = T)
    for (r in 1:nrow(regions_cer_chrm)) {
        strain_index = regions_cer_chrm[r,]$index
        region_start = regions_cer_chrm[r,]$start
        region_end = regions_cer_chrm[r,]$end
        rect(region_start, 
             (strain_index - 1) * row_height + strain_index * padding + vmargin + row_height, 
             region_end, 
             (strain_index - 1) * row_height + strain_index * padding + row_height + vmargin + row_height, 
             col = 'white', border='white')
    }

    # plot unknown regions (unknown state as opposed to just uncalled)
    regions_unk_chrm = merge(strains, regions_unk_chrm, by = 'strain', all.x = T)
    for (r in 1:nrow(regions_unk_chrm)) {
        strain_index = regions_unk_chrm[r,]$index
        region_start = regions_unk_chrm[r,]$start
        region_end = regions_unk_chrm[r,]$end
        rect(region_start, 
             (strain_index - 1) * row_height + strain_index * padding + vmargin + row_height, 
             region_end, 
             (strain_index - 1) * row_height + strain_index * padding + row_height + vmargin + row_height, 
             col = unk_color, border=unk_color)
    }


    # plot par regions
    regions_chrm = merge(strains, regions_chrm, by = 'strain', all.x = T)
    for (r in 1:nrow(regions_chrm)) {
    	mx = NULL 
	#if (unknown) { 
	#    mx = regions_chrm[r,]$number_match_only_par
	#} else {
	#    mx = regions_chrm[r,]$number_match_par_not_cer }
    	#if (mx > 0)
	#{
	#    	cx = vcolors100[log(mx)/log(m) * 100 + 1]
	#}
	#else
	#{
	#	cx = 'black'
	#}
        #######
        cx = par_color
        #######
        strain_index = regions_chrm[r,]$index
        region_start = regions_chrm[r,]$start
        region_end = regions_chrm[r,]$end
        rect(region_start, 
             (strain_index - 1) * row_height + strain_index * padding + vmargin + row_height, 
             region_end, 
             (strain_index - 1) * row_height + strain_index * padding + row_height + vmargin + row_height, 
             col = cx, border=cx)
    }

    # plot strain labels
    positions = seq((padding + row_height) + vmargin + row_height/2,
	    	    nrow(strains)*(padding + row_height) + row_height + vmargin,
    	            padding+row_height)
    axis(2, at=positions, labels=strains$strain, las=1, cex.axis=1, tick=FALSE, line=-.5)

    # plot position labels
    positions = seq(seq_start, seq_end+1, (seq_end - seq_start + 1) / 10)
    axis(1, at=positions, las=1, cex.axis=1, line = 0)

    # plot color scale along top 
    box_height = row_height
    box_width = (seq_end-seq_start+1)/100
    for (ci in 1:100)
    {
	# xleft, ybottom, xright, ytop
	rect(seq_start + (ci - 1) * box_width, total_height - box_height,
             seq_start + ci * box_width, total_height,
	     col=vcolors100[ci])
	# x, y of middle of text
	val = round(exp(ci/100*log(m)))
	text(seq_start + (ci-1)*box_width + box_width/2, total_height - box_height/2, val)
    }

    dev.off()

}
