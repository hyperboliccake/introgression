library(ggplot2)
library(viridis)
require(grDevices)

# read in table of introgressed blocks
# type species start end rep
blocks = read.table('../../results/sim/sim_out_s1_introgressed_blocks_strain_1.txt', sep='\t', header=T)

# read in site codings
# site code rep
codings = read.table('../../results/sim/sim_out_s1_site_codings_strain_1.txt', sep='\t', header=T)

row_height = 2
padding = 3
vmargin = 5 # for top and bottom

vcolors = viridis(7, option = 'plasma')

for (rep in unique(blocks$rep))
{

    print(rep)
    pdf(paste('../../results/sim/introgression_plots/s1_strain_1_rep_', rep, '.pdf', sep=''), width = 16, height = 9)

    par(mar=c(5, 10, 4, 2))

    blocks_subset = blocks[which(blocks$rep == rep),]
    block_types = levels(blocks$type)
    num_block_types = length(block_types)

    codings_subset = codings[which(codings$rep == rep),]
    seq_start = min(codings_subset$site)
    seq_end = max(codings_subset$site)

    # base plot
    # type n doesn't produce any points or lines
    # (left, right), (bottom, top)
    plot(c(seq_start, seq_end+1), 
    	 c(0, num_block_types * (row_height + padding) + row_height + 2 * vmargin),
         type = "n", xlab = "",
         ylab = "", main = "", xaxt='n', yaxt='n', xaxs='i', yaxs='i',
mgp=c(2,2,.5)	 )


    # move x axis label and title closer to axis
    title(xlab = paste("position"), line = 3, cex.lab=1.7)

    # plot codings (bottom row)
    for (i in seq_start:seq_end)
    {
	code = as.character(codings_subset[which(codings_subset$site == i),]$code)
	color <- switch(code,
			'++' = 'gray70',
        		'+-' = vcolors[2],
        		'-+' = vcolors[4],
			'--' = vcolors[6])
	if (code != '++') {
	rect(i, vmargin, i+1, vmargin + row_height, col=color, border=color)}

    }

    # plot actual/predicted introgressed rows
    block_row = 1
    for (block_type in block_types)
    {
	block_type_subset = blocks_subset[which(blocks_subset$type == block_type),]
	for (i in 1:nrow(block_type_subset))
	{
	    # left, bottom, right, top
	    print(block_row*(padding + row_height) + padding)
	    print(block_row*(padding + row_height) + padding + row_height)
	    rect(block_type_subset[i,]$start,
	    	 block_row*(padding + row_height) + vmargin,
        	 block_type_subset[i,]$end,
	    	 block_row*(padding + row_height) + row_height + vmargin,
         	 col = vcolors[2], border = vcolors[2])
	}
	block_row = block_row + 1
    }

    # plot block_type labels
    positions = seq((padding + row_height) + vmargin + row_height/2,
	    	    num_block_types*(padding + row_height) + row_height + vmargin,
    	            padding+row_height)
    axis(2, at=positions, labels=block_types, las=1, cex.axis=1, tick=FALSE, line=-.5)

    # plot position labels
    positions = seq(seq_start, seq_end+1, (seq_end - seq_start + 1) / 10)
    axis(1, at=positions, las=1, cex.axis=1, line = 0)

    # overall outline on top
    rect(seq_start,
         0,
         seq_end + 1,
         num_block_types * (row_height + padding) + row_height + 2 * vmargin,
         col = FALSE, border = "black")

    dev.off()
	 
}


