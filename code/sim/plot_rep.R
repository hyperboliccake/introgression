library(ggplot2)
library(viridis)
require(grDevices)

args = commandArgs(trailingOnly=TRUE)

sim_id = args[1]
pred_id = args[2]
#num_reps = 500
num_reps = 100
reps = 0:(num_reps-1)
threshold = as.numeric(args[4])

prob_height = 20

row_height = 2
padding = 3
vmargin = 1 # for top and bottom

vcolors = viridis(7, option = 'plasma')

prob_labels = c(paste('prob_predicted_', pred_id, '_par', sep=''), 
	      	paste('prob_predicted_phylohmm_', pred_id, '_par', sep=''))
prob_points = c(20, 8)

for (rep in reps)
{
    print(rep)

    # read in input file for this rep
    a = read.table(paste('../../results/sim/', sim_id, '/', pred_id, '/sim_out_', sim_id, '_', pred_id, '_combined_strain_1_rep', rep, '.txt', sep = ''), sep='\t', header=T)

    # output file location and size
    #pdf(paste('../../results/sim/introgression_plots/', sim_id, '_', pred_id, '_strain_1_rep_', rep, '.pdf', sep=''), width = 16, height = 9)
    png(paste('../../results/sim/introgression_plots/', sim_id, '_', pred_id, '_strain_1_rep_', rep, '.png', sep=''), 1600, 900)

    # set margins: bottom, left, top, right
    par(mar=c(5, 10, 4, 2))

    block_types = c('ref', 'actual', paste('predicted_phylohmm_', pred_id, sep=''), paste('predicted_', pred_id, sep=''))
    block_labels = c('reference', 'actual', 'predicted-phylonet', 'predicted-posterior')
    num_block_types = length(block_types)

    seq_start = min(a$site)
    seq_end = max(a$site)

    # base plot
    # type n doesn't produce any points or lines
    # (left, right), (bottom, top)
    # xaxs and yaxs args to specify axes go to edge of plot region
    # xaxt and yaxt args to specify...?
    plot(c(seq_start, seq_end+1), 
    	 c(0, vmargin + row_height + padding + num_block_types * (row_height + padding) + prob_height + vmargin),
         type = "n", xlab = "",
         ylab = "", main = "", xaxt='n', yaxt='n', xaxs='i', yaxs='i', mgp=c(2,2,.5), axes=F)

    # move x axis label and title closer to axis
    title(xlab = paste("position"), line = 3, cex.lab=1.7)


    for (i in 1:nrow(a))
    {
	
        # plot codings (bottom row)
	code = a[i,]$coding
	x = a[i,]$site
	y = vmargin + row_height/8
	color = 'gray70'
	if (code == '+-') {
	   y = y + row_height/4
	   color = vcolors[2]}
	if (code == '-+') {
	   y = y + 2*row_height/4
	   color = vcolors[4]}
	if (code == '--') {
	   y = y + 3*row_height/4
	   color = vcolors[6]}
	points(x, y, col=color, pch=20, alpha=.5)

	# plot all block type rows, starting from bottom
	for (b in 1:num_block_types)
	{
	    species = a[i,block_types[b]]
	    color = vcolors[2]
	    if (species == 'par') {
	       color = vcolors[4]
               starty = vmargin + (padding + row_height) * b
	       endy = starty + row_height
               segments(x, starty, x, endy, col=color, lend=1, alph=.5)
	    }
	}

	# plot probabilities - manual graph, woooo
	graph_bottom = vmargin + (row_height + padding) * (num_block_types + 1)
	for (px in 1:length(prob_labels)) {
	prob_label = prob_labels[px]
	p = a[i,prob_label]
	y = p * prob_height + graph_bottom
	color = vcolors[2]
	if (p > threshold) { 
	    color = vcolors[4] }
	points(x, y, pch = prob_points[px], col=color)
	}	
    }
    
    # a little hacky - show which probability points go with which prediction method
    if (length(prob_labels > 1) {
    for (px in 1:length(prob_labels)) {
    	x = 3
        y = vmargin + (padding + row_height) * 
	    (num_block_types-px+1) + row_height / 2
    	points(x, y, pch = prob_points[px], col=vcolors[2])
    }
    }

    # plot probability threshold
    # segments(seq_start, y, seq_end, y)

    # plot block_type labels
    positions = seq((padding + row_height) + vmargin + row_height/2,
	    	    num_block_types*(padding + row_height) + row_height + vmargin,
    	            padding+row_height)
    axis(2, at=positions, labels=block_labels, las=1, cex.axis=1, tick=FALSE, line=-.5)

    # plot position labels
    positions = seq(seq_start, seq_end+1, (seq_end - seq_start + 1) / 10)
    axis(1, at=positions, las=1, cex.axis=1, line = 0)

    dev.off()
	 
}


