library(ggplot2)
library(viridis)
require(grDevices)

sim_num = '7'

# read in table of introgressed blocks
# type species start end rep
blocks = read.table(paste('../../results/sim/sim_out_s', sim_num, 
       	 		  '_introgressed_blocks_strain_1.txt', sep=''), 
                    sep='\t', header=T)

# read in site codings
# site code rep
codings = read.table(paste('../../results/sim/sim_out_s', sim_num, 
       	 		  '_site_codings_strain_1.txt', sep=''), 
                    sep='\t', header=T)

row_height = 2
padding = 3
vmargin = 3 # for top and bottom

vcolors = viridis(7, option = 'plasma')

block_types_ordered = c('actual_ref', 'actual','predicted','predicted_phylohmm')
block_type_names_ordered = c('reference', 'actual introgression','predicted introgression (my method)','predicted introgression (PhyloNet-HMM)')

for (rep in unique(blocks$rep))
{

    print(rep)

    pdf(paste('../../results/sim/introgression_plots/s', sim_num, 
    	      '_strain_1_rep_', rep, '.pdf', sep=''), width = 16, height = 9)
    par(mar=c(5, 20, 4, 2))

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
    for (block_type in block_types_ordered)
    {
	block_type_subset = blocks_subset[which(blocks_subset$type == block_type),]
	for (i in 1:nrow(block_type_subset))
	{
	    # left, bottom, right, top
	    print(block_row*(padding + row_height) + padding)
	    print(block_row*(padding + row_height) + padding + row_height)

            color <- switch(as.character(block_type_subset[i,]$species),
        		'cer' = vcolors[2],
        		'par' = vcolors[4],
			'bay' = vcolors[6])
	
	    rect(block_type_subset[i,]$start,
	    	 block_row*(padding + row_height) + vmargin,
        	 block_type_subset[i,]$end,
	    	 block_row*(padding + row_height) + row_height + vmargin,
         	 col = color, border = FALSE)

	    # should really do these lines all after or before 
	    # so they consistently go over/under
	    #abline(v = block_type_subset[i,]$start, col = 'gray50')

	}
	block_row = block_row + 1
    }

    # plot block_type labels
    positions = seq((padding + row_height) + vmargin + row_height/2,
	    	    num_block_types*(padding + row_height) + row_height + vmargin,
    	            padding+row_height)
    axis(2, at=positions, labels=block_type_names_ordered, las=1, cex.axis=1, tick=FALSE, line=-.5)

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


