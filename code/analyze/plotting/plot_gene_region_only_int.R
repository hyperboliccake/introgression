# like plot_genes, except we want to plot a whole set of genes together to show a complete region
# also plot gaps and variants colored by whether they match each reference

library(ggplot2)
library(viridis)
require(grDevices)
source('../my_color_palette.R')

options(stringsAsFactors = FALSE) # fuck you R

# r2259 -> SUL1
#region_start = 787000
#region_end = 794000
#chrm = 'II'

## r4560 -> SIR1
region_start = 917571 - 100
region_end = 921647 + 100
chrm = 'IV'

# save plot
fn = paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/plots/region_chr', chrm, '_', region_start, '-', region_end,'.pdf', sep='')
pdf(fn, width = 10, height = 5)

# get strains
strains = read.table('strains.txt', header=F)
strains = t(strains)
row.names(strains) = NULL
strains = data.frame(strain=strains[3:nrow(strains),1])

#strains = data.frame(strain=c('yjm450', 'yjm320', 'yjm1399', 'yjm1355'))#
strains = data.frame(strain=c('yjm1304', 'yjm1202', 'yjm1199', 'yjm681'))#
strains$index=1:nrow(strains)
print(strains)#

## formatting parameters
#context_length = 200
genome_height = 1
num_strains = nrow(strains)
strain_height = 1
padding = .2
genome_width = genome_height/2
hmargin = .2 # for top and bottom

## plot overall outline, line for each strain etc

## base plot
## type n doesn't produce any points or lines
## first strain is at padding
plot(c(region_start, region_end),
     c(-genome_height - padding - hmargin,
       num_strains * (padding + strain_height) + hmargin),
     type = "n", xlab = "", 
     ylab = "", main = "", xaxt='n', yaxt='n', xaxs='i', yaxs='i')

## move x axis label and title closer to axis
title(xlab = paste("Position on chromosome", chrm), line = 1.8, cex.lab=1.5)

## plot overall line for genome at bottom
rect(region_start, 
     -genome_height/2 - padding - genome_width/2, 
     region_end, 
     -genome_height/2 - padding + genome_width/2, 
     col = "gray50", border = "gray50")

## move position labels closer
axis(1, mgp=c(3, .5, 0))

## plot strain labels
positions = seq(padding + strain_height/2,
                (num_strains*(padding + strain_height)),
                (padding+strain_height))
axis(2, at=positions, labels=toupper(strains$strain[order(strains$index)]),
     las=1, cex.axis=1, mgp=c(3, .1, 0), tick=FALSE)
## plot gridlines for strains
for (p in positions) {   
    abline(a=p, b=0, col=alpha(my_color_palette[['nonintrogressed']],0), lwd=2)
}

## plot genes

# get set of all genes

fn = paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/S288c_chr',chrm,'_genes.txt', sep='')
gene_coords = read.table(fn, header=F, stringsAsFactors=F)
names(gene_coords) = c('gene', 'start', 'end')
genes = gene_coords[which((gene_coords$start > region_start &
                           gene_coords$start < region_end) |
                          (gene_coords$end > region_start &
                           gene_coords$end < region_end)),]

for (i in 1:nrow(genes)) {

    gene_name = genes[i,]$gene
    gene_start = max(region_start, genes[i,]$start)
    gene_end = min(region_end, genes[i,]$end)
    gene_length = gene_end - gene_start + 1

    # gene box
    rect(gene_start,
         -genome_height - padding,
         gene_end,
         -padding,
         col = 'white', border = 'white')
    rect(gene_start,
         -genome_height - padding,
         gene_end,
         -padding,
         col = alpha(my_color_palette[['nonintrogressed']],1),
         border = my_color_palette[['nonintrogressed']])
    # gene name
    text((gene_start + gene_end) / 2, -genome_height/2-padding,
         labels = gene_name, adj=c(.5,.5), cex=1.3, col="white")
    
}

## plot variants
fn = 'gene_region_variants.txt'
variants = read.table(fn, header=T, stringsAsFactors=F)

for (i in 1:nrow(variants)) {

    if (i %% 100 == 0) {
        print(i)
    }
    for (j in strains$index) {
        
        strain = strains[which(strains$index == j),]$strain
        x = variants[[strain]][i]
        ps = variants$ps[i]
        if (x == 'p') {
            #points(ps, positions[j], col = alpha(my_color_palette[['introgressed']],.5), pch='|',cex=.3)
            segments(ps, positions[j] - strain_height / 3, y1=positions[j] + strain_height/3, col = alpha(my_color_palette[['introgressed']],.5), lend='butt', lwd=.5)
        }
        else if (x == 'c') {
            #points(ps, positions[j], col = alpha(my_color_palette[['nonintrogressed']],.5),pch='|',cex=.3)
            segments(ps, positions[j] - strain_height / 3, y1=positions[j] + strain_height/3, col = alpha(my_color_palette[['nonintrogressed']],.5), lend='butt', lwd=.5)
        }
        else if (x == 'n') {
            #points(ps, positions[j], col = alpha('gray50',.5),pch='|',cex=.3)
            segments(ps, positions[j] - strain_height / 3, y1=positions[j] + strain_height/3, col = alpha('black',.5), lend='butt', lwd=.5)
        }
        else if (x == '-') {
            segments(ps, positions[j] - strain_height / 3, y1=positions[j] + strain_height/4, col = alpha('black',0), lend='butt', lwd=.5)
            #segments(ps-.5, positions[j], x1=ps+.5, col = 'white',lend='butt')
        }
    }
}

## plot boxes around introgressed regions
fn = '/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/u3_i.001_tv_l1000_f.01/introgressed_blocks_filtered_par_u3_i.001_tv_l1000_f.01_summary_plus.txt'
region_summary = read.table(fn, stringsAsFactors=F, header=T)
for (i in 1:nrow(region_summary)) {
    if (region_summary[i,]$chromosome == chrm) {
    
        r_start = region_summary[i,]$start
        r_end = region_summary[i,]$end

        if ((r_start > region_start & r_start < region_end) | 
            (r_end > region_start & r_end < region_end)) {
        
            r_start = max(r_start, region_start)
            r_end = min(r_end, region_end)

            strain_index = strains[which(strains$strain == region_summary[i,]$strain),]$index
            print(strain_index)
            print(region_summary[i,]$strain)
            print(strains)
            print(which(strains$strain == region_summary[i,]$strain))
            rect(r_start,
                positions[strain_index]-strain_height/2,
                r_end,
                positions[strain_index]+strain_height/2,
                col = alpha(my_color_palette[['introgressed']],0), border=my_color_palette[['introgressed']])
            }
        }
    }

## overall outline on top
rect(region_start, 
     -genome_height - padding - hmargin,
     region_end, 
     num_strains * (padding + strain_height) + hmargin,	 
     col = FALSE, border = "black")

dev.off()





    
