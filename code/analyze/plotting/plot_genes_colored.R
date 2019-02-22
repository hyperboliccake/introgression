library(ggplot2)
library(viridis)
require(grDevices)

chrms = c('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI')

#gene_summary = read.table('introgressed_hmm_9_genes_summary_filtered.txt', sep='\t', header=T)
genes_int = read.table('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/genes_strain_hist_u3_i.001_tv_l1000_f.01.txt', sep='\t', header=T)
gene_summary = data.frame()
for (chrm in chrms) {
    g = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/S288c_chr', chrm, '_genes.txt', sep=''), header=F, sep='\t', stringsAsFactors=F)
    names(g) = c("gene", "start", "end")
    g$chromosome = chrm
    gene_summary = rbind(gene_summary, g)
}
gene_summary = gene_summary[which(gene_summary$gene %in% genes_int$gene),]
print(nrow(gene_summary))
print(nrow(genes_int))
region_summary = read.table('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/introgressed_blocks_filtered_par_u3_i.001_tv_l1000_f.01_summary_plus.txt', sep='\t', header=T)
strains = read.table('strain_list.txt', header=F)
strains = data.frame(strain=strains[,1], index=1:nrow(strains))
region_summary = merge(region_summary, strains, all.x=T)

#strain_info = read.table('../../100_genomes_info.txt', sep='\t', header=T)
#strain_info$strain = tolower(strain_info$strain)
#strain_info$is_beer_strain = grepl('(beer|Beer)', strain_info$environmental_origin)
## ginger beer is not beer :P
#strain_info[which(grepl('(ginger|Ginger)',strain_info$environmental_origin)),]$is_beer_s#train = FALSE
#strain_info$is_wine_strain = grepl('(wine|Wine)', strain_info$environmental_origin)
#strain_info$is_cider_strain = grepl('(cider|Cider)', strain_info$environmental_origin)
#strain_info$is_other_alcohol_strain = grepl('(ferment|Ferment|distill|Distill)', strain_#info$environmental_origin)
#strain_env = data.frame(strain=strain_info$strain, 
#	                is_beer_strain=strain_info$is_beer_strain, 
#			is_wine_strain=strain_info$is_wine_strain, 
#			is_cider_strain=strain_info$is_cider_strain, 
#			is_other_alcohol_strain=strain_info$is_other_alcohol_strain)

#strains = merge(strains, strain_env, by = 'strain', all.x = T)

context_length = 200
gene_height = 17
num_strains = nrow(strains)
strain_height = 2
padding = 1
genome_width = 2 
hmargin = 5 # for top and bottom

vcolors = viridis(7, option = 'plasma')

for (i in 1:nrow(gene_summary)) {
   
    gene_start = gene_summary[i,]$start
    gene_end = gene_summary[i,]$end
    gene_length = gene_end - gene_start + 1
    gene_name = gene_summary[i,]$gene

    print(gene_name)
    pdf(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/gene_plots/', gene_name, '.pdf', sep=''))

    # base plot
    # type n doesn't produce any points or lines
    # first strain is at padding
    plot(c(-context_length + gene_start, gene_end + context_length), 
         c(-gene_height + 1 - padding * 2 - hmargin,
           num_strains * (padding + strain_height) + hmargin), 
         type = "n", xlab = "", 
	 ylab = "", main = "", xaxt='n', yaxt='n', xaxs='i', yaxs='i')

    # move x axis label and title closer to axis
    title(xlab = paste("position on chromosome", gene_summary[i,]$chromosome), line = 1.8)
    

    # plot intergenic lines
    rect(-context_length + gene_start, 
         -gene_height/2 - padding - genome_width/2, 
         gene_start, 
         -gene_height/2 - padding + genome_width/2, 
         col = "gray50", border = "gray50")
    rect(gene_end, 
         -gene_height/2 - padding - genome_width/2, 
         gene_end + context_length, 
         -gene_height/2 - padding + genome_width/2, 
         col = "gray50", border = "gray50")

    # plot gene
    rect(gene_start,
         -gene_height - padding,
         gene_end,
         -padding,
         col = alpha(vcolors[4],.8), border = vcolors[4])
    # gene name
    text((gene_start + gene_end) / 2, -gene_height/2-padding, labels = gene_name, adj=c(.5,.5), cex=1.3)

    # move position labels closer
    axis(1, mgp=c(3, .5, 0))

    # plot strain labels
    positions = seq(padding + strain_height/2, (num_strains*(padding + strain_height)),
    	        (padding+strain_height))
    axis(2, at=positions, labels=toupper(strains$strain[order(strains$index)]), las=1, cex.axis=.35, mgp=c(3, .1, 0), tick=FALSE)
    # plot gridlines for strains
    for (p in positions)
    {   
        abline(a=p, b=0, col=alpha(vcolors[2],.3))
    }

    #####
    # plot strain introgressed regions
    #####

    #regions = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/genes/', gene_name, '.txt', sep=''), header=T)
    #regions = merge(strains, regions, by = 'strain', all.x = T)

    chrm = gene_summary[i,]$chromosome
    regions = region_summary[which(region_summary$chromosome == chrm),]
    
    for (r in 1:nrow(regions)) {
	strain_index = regions[r,]$index
	region_start = regions[r,]$start
	region_end = regions[r,]$end
	environment_color = vcolors[6]
	environment_border = vcolors[6]
	# viridis color palettes
	#if(regions[r,]$is_beer_strain) {
	#	environment_color = '#482677FF'
	#	environment_border = '#481567FF'
	#}
	#if(regions[r,]$is_wine_strain) {
	#	environment_color = '#2D708EFF'
	#	environment_border = '#33638DFF'
	#}
	#if(regions[r,]$is_cider_strain) {
	#	environment_color = '#29AF7FFF'
	#	environment_border = '#20A387FF'
	#}
	#if(regions[r,]$is_other_alcohol_strain) {
	#	environment_color = '#B8DE29FF'
	#	environment_border = '#95D840FF'
	#}
    	rect(region_start, 
             (strain_index - 1) * strain_height + strain_index * padding, 
             region_end, 
             (strain_index - 1) * strain_height + strain_index * padding + strain_height, 
	     col = environment_color, border = environment_border)

        x1 = region_start
        x2 = region_end

        if ((x1 >= gene_start & x1 <= gene_end) |
            (x2 >= gene_start & x2 <= gene_end) |
            (x1 < gene_start & x2 > gene_end)) {

            x1 = max(x1, gene_start)
            x2 = min(x2, gene_end)
            
            text((x1 + x2) / 2,
            (strain_index - 1) * strain_height + strain_index * padding + strain_height / 2,
            adjust = c(.5,.5), labels = regions[r,]$region_id, cex=.3)
        }

    }

    # overall outline on top
    rect(-context_length + gene_start, 
         -gene_height + 1 - padding * 2 - hmargin,
         gene_end + context_length, 
	 num_strains * (padding + strain_height) + hmargin,	 
         col = FALSE, border = "black")

    dev.off()
    # plot sites that match introgressed reference

}








