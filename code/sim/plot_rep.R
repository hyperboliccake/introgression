library(ggplot2)
library(viridis)
require(grDevices)

sim_id = 'p1'
pred_id = 'pred1'
num_reps = 500

prob_height = 20

row_height = 2
padding = 3
vmargin = 5 # for top and bottom

vcolors = viridis(7, option = 'plasma')

for (rep in 0:(num_reps-1))
{
    print(rep)

    # read in input file for this rep
    a = read.table(paste('../../results/sim/sim_out_', sim_id, '_', pred_id, '_combined_strain_1_rep', rep, '.txt', sep = ''), sep='\t', header=T)

    # output file location and size
    pdf(paste('../../results/sim/introgression_plots/', sim_id, '_', pred_id, '_strain_1_rep_', rep, '.pdf', sep=''), width = 16, height = 9)

    # set margins: bottom, left, top, right
    par(mar=c(5, 10, 4, 2))

    block_types = c('ref', 'actual', paste('predicted_', pred_id, sep=''))
    block_labels = c('reference', 'actual', 'predicted')
    num_block_types = length(block_types)

    seq_start = min(a$site)
    seq_end = max(a$site)

    # base plot
    # type n doesn't produce any points or lines
    # (left, right), (bottom, top)
    # xaxs and yaxs args to specify axes go to edge of plot region
    # xaxt and yaxt args to specify no axes drawn
    plot(c(seq_start, seq_end+1), 
    	 c(0, row_height + padding + num_block_types * (row_height + padding) + prob_height),
         type = "n", xlab = "",
         ylab = "", main = "", xaxt='n', yaxt='n', xaxs='i', yaxs='i', mgp=c(2,2,.5))

    # move x axis label and title closer to axis
    title(xlab = paste("position"), line = 3, cex.lab=1.7)

    for (i in 1:nrow(a))
    {
        # plot codings (bottom row)
	code = a[i,]$code
	color <- switch(code,
			'++' = NA,
        		'+-' = vcolors[2],
        		'-+' = vcolors[4],
			'--' = vcolors[6])
	# startx, starty, endx, endy
	x = a[i,]$site
	starty = vmargin
	endy = starty + row_height
	segments(x, starty, x, endy, col=color)

	# plot all block type rows, starting from bottom
	for (b in 1:num_block_types)
	{
	    species = a[i,block_types[i]]
            color <- switch(code,
	    	   	    'cer' = vcolors[2],
        		    'par' = vcolors[4])
	    starty = vmargin + row_height + padding * i
	    endy = starty + row_heigt
            segments(x, starty, x, endy, col=color)

	}

	# plot probabilities - manual graph, woooo
    }

    # plot block_type labels
    positions = seq((padding + row_height) + vmargin + row_height/2,
	    	    num_block_types*(padding + row_height) + row_height + vmargin,
    	            padding+row_height)
    axis(2, at=positions, labels=block_types, las=1, cex.axis=1, tick=FALSE, line=-.5)

    # plot position labels
    positions = seq(seq_start, seq_end+1, (seq_end - seq_start + 1) / 10)
    axis(1, at=positions, las=1, cex.axis=1, line = 0)

    dev.off()
	 
}


