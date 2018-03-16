require(grDevices)

gene ='SUL1'
chrm = 'II'

fn = paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/genes_for_each_region_chr',chrm,'_u3_i.001_tv_l1000_f.01.txt', sep='')
regions_genes = read.table(fn, sep=c('\t', ' '), fill=T, col.names=1:max(count.fields(fn, sep = c("\t", ' '))))

fn = '/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/introgressed_blocks_filtered_par_u3_i.001_tv_l1000_f.01_summary_plus.txt'
region_summary = read.table(fn, stringsAsFactors=F, header=T)

regions_with_gene = regions_genes[which(apply(regions_genes, 1, function(r) any(r %in% c(gene)))),1]

strains = read.table('strains.txt', header=F)
strains = data.frame(strain=t(strains[3:length(strains)]))
strains$index=1:nrow(strains)

fn = paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/S288c_chr',chrm,'_genes.txt', sep='')
gene_coords = read.table(fn)

context_length = 200
gene_height = 10
num_strains = nrow(strains)
strain_height = 2
padding = 1
genome_width = 2

#for (i in 1:nrow(gene_summary)) {

gene_name = gene
i = which(gene_coords[,1]==gene)

gene_start = gene_coords[i,2]#$start
gene_end = gene_coords[i,3]#$end
gene_length = gene_end - gene_start + 1
    #gene_name = gene_summary[i,]$gene


    print(gene_name)
fn = paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/genes/',gene_name,'.pdf', sep='')
    pdf(fn)
    
    
    ## base plot
    ## type n doesn't produce any points or lines
    ## first strain is at padding
    plot(c(-context_length + gene_start, gene_end + context_length),
         c(-gene_height + 1 - padding * 2, num_strains * (padding + strain_height)),
         type = "n", xlab = paste("chromosome", chrm),
         ylab = "", main = gene_name, yaxt='n')

    ## plot genome line
    rect(-context_length + gene_start,
         -gene_height/2 - padding - genome_width/2,
         gene_end + context_length,
         -gene_height/2 - padding + genome_width/2,
         col = "black", border = "black")

    ## plot gene
    rect(gene_start,
         -gene_height - padding,
         gene_end,
         -padding,
         col = "darkorchid1", border = "darkorchid4")


    ## plot strain labels
    positions = seq(padding + strain_height/2, (num_strains*(padding + strain_height)),
    (padding+strain_height))
    axis(2, at=positions, labels=strains$strain, las=1, cex.axis=.35)

    ##=======
    ## plot strain introgressed regions
    ##=======
    #regions = read.table(paste('genes/', gene_name, '.txt', sep=''), header=T)
    #regions = merge(strains, regions, by = 'strain', all.x = T)
    for (r in regions_with_gene) {
        ri = which(r==region_summary$region_id)
        strain_index = strains[which(strains$strain == region_summary[ri,]$strain),]$index
        print(strains)
        print(strain_index)
        region_start = region_summary[ri,]$start
        region_end = region_summary[ri,]$end
        print((strain_index - 1) * strain_height + strain_index * padding)
        rect(region_start,
        (strain_index - 1) * strain_height + strain_index * padding,
        region_end,
        (strain_index - 1) * strain_height + strain_index * padding + strain_height,
        col = 'dodgerblue3')
    }
    
    dev.off()
    ## plot sites that match introgressed reference
#}
