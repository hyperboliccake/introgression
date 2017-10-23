library(ggplot2)
library(viridis)
require(grDevices)

# read in table of introgressed blocks
# type species start end rep
blocks = read.table('sim_out_s1_introgressed_blocks_strain_1.txt', sep='\t', header=T)

# read in site codings
# site code rep
codings = read.table('sim_out_s1_site_codings_strain_1.txt', sep='\t', header=T)

row_height = 5
padding = 1
vmargin = 5 # for top and bottom

vcolors = viridis(7, option = 'plasma')

for (rep in codings$rep)
{
    blocks_subset = blocks[which(blocks$rep == rep),]
    num_block_types = length(levels(blocks$type))

    codings_subset = codings[which(codings$rep == rep),]
    seq_start = min(codings_subset$site)
    seq_end = max(codings_subset$site)

    # base plot
    # type n doesn't produce any points or lines
    # (left, right), (bottom, top)
    plot(c(seq_start, seq_end), 
    	 c(0, num_block_types * (row_height + padding) + row_height + 2 * vmargin),
         type = "n", xlab = "",
         ylab = "", main = "", xaxt='n', yaxt='n', xaxs='i', yaxs='i')

    # move x axis label and title closer to axis
    title(xlab = paste("position"), line = 1.8)

    # plot codings (bottom row)
    for (i in seq_start:seq_end)
    {
	code = blocks_subset[which(blocks_subset$site == i),]$code
	color <- function(type) {
 	      	 switch(code,
			'++' = 'gray50',
        		'+-' = vcolors[2],
        		'-+' = vcolors[4],
			'--' = vcolors[6])
	}


    }

    # plot actual/predicted introgressed rows
    block_row = 0
    for (block_type in levels(blocks_subsest))
    {
	block_type_subset = blocks_subset[which(blocks_subset$types == block_type),]
	for (i in 1:nrow(block_type_subset))
	{
	    # left, bottom, right, top
	    rect(block_type_subset[i,]$start,
	    	 block_row*(padding + row_height) + padding,
        	 block_type_subset[i,]end,
	    	 block_row*(padding + row_height) + padding + row_height,
         	 col = alpha(vcolors[2],1), border = vcolors[2])
	}
    }
    	

    # overall outline on top
    rect(seq_start,
         0,
         seq_end,
         num_block_types * (row_height + padding) + row_height + 2 * vmargin,
         col = FALSE, border = "black")

    dev.off()
	 
}


