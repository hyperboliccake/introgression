library(ggplot2)
library(viridis)
require(grDevices)

unknown = FALSE

vcolors = viridis(7, option = 'plasma')
vcolors100 = viridis(100, option='plasma')
m = 5000

par_color = vcolors[6]
cer_color = vcolors[2]
both_color = vcolors[5]
neither_color = 'gray70'
unk_color = viridis(5)[3]

tag = commandArgs(trailingOnly=TRUE)[1]
region = commandArgs(trailingOnly=TRUE)[2]


regions = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/introgressed_blocks_par_', tag, '_summary.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)

mx = NULL
if (grepl('u', tag)) {
mx = regions[which(regions$region_id==region),]$number_match_only_par

} else {
mx = regions[which(regions$region_id==region),]$number_match_par_not_cer 
}
cx = vcolors100[log(mx)/log(m) * 100 + 1]
intd_color = cx


a = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', region, '/', 'ps_annotations.txt',sep=''), stringsAsFactors=FALSE, sep='\t', header=TRUE)

gene_height = 3
match_height = 5
intd_height = 3
prob_height = 0
padding = 1
vmargin = 1 # for top and bottom

# base plot
png(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/plots/', region, '/', region, '_', tag, '.png', sep=''), 1280, 720)


seq_start = min(a$ps)
seq_end = max(a$ps)

total_height = vmargin * 2 + gene_height + match_height + intd_height + prob_height + padding * 3
    plot(c(seq_start, seq_end+1), 
    	 c(0, total_height),
         type = "n", xlab = "",
         ylab = "", main = "", xaxt='n', yaxt='n', xaxs='i', yaxs='i', mgp=c(2,2,.5), axes=F)

    # move x axis label and title closer to axis
    #title(xlab = paste("chromosome ", chrm, " ", m, sep = ''), line = 3, cex.lab=1.7)


# plot matching


	# plot gene row
	all_genes = unique(a$gene)
	all_genes = all_genes[which(all_genes != '')]
	for (gene in all_genes) {
	    gene_start = min(a[which(a$gene == gene),]$ps)
	    gene_end = max(a[which(a$gene == gene),]$ps)
	    starty = vmargin
	    endy = vmargin + gene_height
            rect(gene_start, 
                 starty,
                 gene_end, 
             	 endy,
                 col = both_color, border='black')
           # x, y of middle of text
	   text((gene_start + gene_end)/2, (starty+endy)/2, gene)
	    
	}

    for (i in 1:nrow(a))
    {
	

        # plot codings (bottom row)
	code = a[i,]$match
	x = a[i,]$ps
	y = vmargin + gene_height + padding + match_height/8
	color = both_color
	if (code == 'S ') {
	   y = y + match_height/4
	   color = cer_color}
	if (code == ' C') {
	   y = y + 2*match_height/4
	   color = par_color}
	if (code == '..') {
	   y = y + 3*match_height/4
	   color = neither_color}
	points(x, y, col=color, pch=19)


	# plot introgression
	if (a[i,]$intd %in% c('i', 'I')) {
	   cx = intd_color
	   if (a[i,]$intd == 'i'){
	        cx = 'gray70'
	   }
               starty = vmargin + gene_height + padding + match_height + padding
	       endy = starty + intd_height
               segments(x, starty, x, endy, col=cx, lend=1, lwd=3)
	    }

	# plot probabilities - manual graph, woooo
	#graph_bottom = vmargin + (row_height + padding) * (num_block_types + 1)
	#for (px in 1:length(prob_labels)) {
	#prob_label = prob_labels[px]
	#p = a[i,prob_label]
	#y = p * prob_height + graph_bottom
	#color = vcolors[2]
	#if (p > threshold) { 
	#    color = vcolors[4] }
	#points(x, y, pch = prob_points[px], col=color)
	#}	
    }
    
    # plot probability threshold
    # segments(seq_start, y, seq_end, y)

    # plot block_type labels
    #positions = seq((padding + row_height) + vmargin + row_height/2,
    #	    num_block_types*(padding + row_height) + row_height + vmargin,
    #	            padding+row_height)
    #axis(2, at=positions, labels=block_labels, las=1, cex.axis=1, tick=FALSE, line=-.5)

    # plot position labels
    np = (seq_end - seq_start) / 5
    positions = round(seq(seq_start, seq_end+1, np))
    axis(1, at=positions, las=1, cex.axis=1, line = 0)

    dev.off()







asdgsag





for (ci in 1:length(chrms)) 
{
    chrm = chrms[ci]
    print(chrm)

    regions_chrm = regions[which(regions$chromosome == chrm),]
    regions_cer_chrm = regions_cer[which(regions_cer$chromosome == chrm),]
    regions_unk_chrm = regions_unk[which(regions_unk$chromosome == chrm),]

    png(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/plots/chrm', chrm, '_', tag, '.png', sep=''), 6400, 3600)

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
	if (unknown) { 
	    mx = regions_chrm[r,]$number_match_only_par
	} else {
	    mx = regions_chrm[r,]$number_match_par_not_cer }
    	if (mx > 0)
	{
	    	cx = vcolors100[log(mx)/log(m) * 100 + 1]
	}
	else
	{
		cx = 'black'
	}
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