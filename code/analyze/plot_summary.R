library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)
source('../my_color_palette.R')

args = commandArgs(trailingOnly=TRUE)
tag = args[1]
suffix = ''
if (length(args) == 2)
{
    suffix = args[2]
}


regions = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/introgressed_blocks', suffix,'_par_', tag, '_summary_plus.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)
regions$variable = "variable" # make things a lil easier with dcast
regions$overlap_gene = regions$number_genes >= 1
regions$length = regions$end - regions$start + 1


# bar chart: fraction genome introgressed per strain
g = 12071326
a = dcast(regions, strain ~ variable, value.var="length", fun.aggregate=sum)
a$frac = a$variable/g
ggplot(a, aes(x=reorder(strain, -frac), y=frac, fill='x')) + 
    geom_bar(stat='identity',position='dodge') + 
    xlab('Strain') + ylab('Fraction of genome introgressed') +
    scale_fill_manual(values =c(my_color_palette[['introgressed']]) ) +
    guides(fill=FALSE) +
    scale_y_continuous(expand = c(0,0), limits=c(0,max(a$frac))) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.ticks=element_line(colour="black"),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(size=8,angle = 45,vjust = 1,hjust=1,colour="black"), 
          axis.text.y = element_text(colour="black"))
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/frac_bar',suffix,'_',tag,'.pdf',sep=''), width = 12, height = 6)
print(mean(a$frac))
print(median(a$frac))
print(min(a$frac))
print(max(a$frac))


# bar chart: genes introgressed per strain
a = dcast(regions, strain ~ variable, value.var="number_genes", fun.aggregate=sum)
ggplot(a, aes(x=reorder(strain, -variable), y=variable, fill='x')) + 
    geom_bar(stat='identity',position='dodge') + 
    xlab('Strain') + ylab('Number of genes introgressed') +
    scale_fill_manual(values =c(my_color_palette[['introgressed']]) ) +
    guides(fill=FALSE) +
    scale_y_continuous(expand = c(0,0), limits=c(0,max(a$variable))) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.ticks=element_line(colour="black"),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(size=8,angle = 45,vjust = 1,hjust=1,colour="black"), 
          axis.text.y = element_text(colour="black"))
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/num_genes_bar',suffix,'_',tag,'.pdf',sep=''), width = 12, height = 6)


# histogram: region lengths
ggplot(regions, aes(x=length, fill='x')) + geom_histogram(binwidth=100) +
    xlab('Region length') + ylab('Number of regions') +
    scale_fill_manual(values =c(my_color_palette[['introgressed']]) ) +
    guides(fill=FALSE) +
    scale_y_continuous(expand = c(0,0), limits=c(0,1500)) + 
    scale_x_continuous(expand = c(0,0), limits=c(0,max(regions$length)*1.05)) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.line=element_line(),
          axis.ticks=element_line(colour="black"),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(colour="black",size=12), 
          axis.text.y = element_text(colour="black",size=12))
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/length_hist',suffix,'_',tag,'.pdf',sep=''), width = 9, height = 6)
print(mean(regions$length))
print(median(regions$length))
print(min(regions$length))
print(max(regions$length))
print(length(regions$length))

# region lengths zoomed
ggplot(regions, aes(x=length, fill='x')) + geom_histogram(binwidth=100) +
    xlab('Region length') + ylab('Number of regions') +
    scale_fill_manual(values =c(my_color_palette[['introgressed']]) ) +
    guides(fill=FALSE) +
    scale_y_continuous(expand = c(0,0), limits=c(0,42)) + 
    scale_x_continuous(expand = c(0,0), limits=c(0,max(regions$length)*1.05)) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.line=element_line(),
          axis.ticks=element_line(colour="black"),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(colour="black",size=12), 
          axis.text.y = element_text(colour="black",size=12))
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/length_hist_zoomed',suffix,'_',tag,'.pdf',sep=''), width = 9, height = 6)


# histogram: region lengths (log scale)
ggplot(regions, aes(x=length, fill='x')) + geom_histogram(binwidth=50) +
    xlab('Region length') + ylab('Number of regions') +
    scale_fill_manual(values =c(my_color_palette[['introgressed']]) ) +
    guides(fill=FALSE) +
    #scale_y_continuous(expand = c(0,0), limits=c(-1,750)) + 
    scale_y_log10(limits=c(1,1000)) +
    scale_x_continuous(expand = c(0,0), limits=c(0,max(regions$length)*1.05)) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.line=element_line(),
          axis.ticks=element_line(colour="black"),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(colour="black",size=12), 
          axis.text.y = element_text(colour="black",size=12))
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/length_hist_log',suffix,'_',tag,'.pdf',sep=''), width = 9, height = 6)
print(mean(regions$length))
print(median(regions$length))
print(min(regions$length))
print(max(regions$length))
print(length(regions$length))




asdgasdg

predict_args = read.table('predict_args.txt', sep=' ', stringsAsFactors=F)
predict_args = predict_args[c(seq(1,19),seq(21,25),seq(27,36)),]
tags = predict_args[,1]

##======
# plot: lengths of all introgressed regions
##======

print('plotting lengths of introgressed regions')

#----
# for each tag, make histogram - regions vs lengths
#----

for (tag in tags) {
    print(tag)
    dir.create(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/',sep=''),showWarnings=FALSE)


    a = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plot_region_lengths_',tag,'.txt',sep=''), sep='\t')
    names(a) = c('region_id', 'strain', 'chromosome', 'region_length')

    plot = ggplot(a, aes(region_length/1000, fill='a')) + 
    	      geom_histogram(binwidth=1) +
	      xlab('region length (kb)') + ylab('number of regions') +
	      scale_fill_viridis(discrete=TRUE) +
	      guides(fill=FALSE) +
	      scale_x_continuous(expand = c(0,0), limits=c(-1,max(a$region_length)/1000)) +
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))

    plot + scale_y_continuous(expand = c(0,0),
               limits=c(0,max(ggplot_build(plot)$data[[1]]$count)*1.1))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/region_length_hist_',tag,'.pdf',sep=''), width = 12, height = 7)


    a = a[which(a$region_length < 1000),]
    plot = ggplot(a, aes(region_length, fill='a')) + 
    	      geom_histogram(binwidth=1) +
	      xlab('region length (kb)') + ylab('number of regions') +
	      scale_fill_viridis(discrete=TRUE) +
	      guides(fill=FALSE) +
	      #scale_x_continuous(expand = c(0,0), limits=c(-1,max(a$region_length)/1000)) +
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))

    #plot + scale_y_continuous(expand = c(0,0),
    #           limits=c(0,max(ggplot_build(plot)$data[[1]]$count)*1.1))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/region_small_length_hist_',tag,'.pdf',sep=''), width = 12, height = 7)

}

sdagdsg

#----
# plot average region length for all tags (or boxplot?)
#----

a = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/plot_region_lengths.txt',sep=''), sep='\t')
names(a) = c('tag', 'improvement_frac', 'threshold', 'expected_length', 'expected_frac', 'mean', 'bs_lower', 'bs_upper', 'median', 'min', 'max')

#d = dcast(d, threshold + expected_length + expected_frac ~ improvement_frac, value.var='mean')

# group by threshold
ggplot(a, aes(x=tag, y=mean/1000)) + 
    	      geom_bar(aes(fill=as.factor(threshold)), stat='identity',position='dodge') +
	      ylab('average region length (kb)') + xlab('parameter set') +
	      scale_fill_viridis(discrete=TRUE) +
	      labs(fill='improvement fraction') +
	      scale_y_continuous(expand = c(0,0), limits=c(0,max(a$mean)/1000)) +
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(angle = 45,vjust = 1,hjust=1,colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/region_length_t.pdf', width = 12, height = 7)


##======
# plot: number of genes per introgressed region
##======

print('plotting number of genes per introgressed region')

#----
# for each tag, make histogram - number of regions vs number of genes
#----

for (tag in tags) {
    print(tag)
    dir.create(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/',sep=''),showWarnings=FALSE)


    a = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plot_number_genes_by_region_',tag,'.txt',sep=''), sep='\t')
    names(a) = c('region_id', 'strain', 'chromosome', 'num_genes')

    plot = ggplot(a, aes(num_genes, fill='a')) + 
    	      geom_histogram(binwidth=1) +
	      scale_fill_viridis(discrete=TRUE) +
	      guides(fill=FALSE) +
	      scale_x_continuous(expand = c(0,0), limits=c(-1,max(a$num_genes))) +
	      xlab('number of genes per region') + ylab('number of regions') +
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))

    plot + scale_y_continuous(expand = c(0,0),
               limits=c(0,max(ggplot_build(plot)$data[[1]]$count)*1.1))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/number_genes_in_region_hist_',tag,'.pdf',sep=''), width = 12, height = 7)
}

#----
# plot average number of genes per region for all tags (or boxplot?)
#----

a = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/plot_number_genes_by_region.txt',sep=''), sep='\t')
names(a) = c('tag', 'improvement_frac', 'threshold', 'expected_length', 'expected_frac', 'mean', 'bs_lower', 'bs_upper', 'median', 'min', 'max')

ggplot(a, aes(x=tag, y=mean, fill=as.factor(improvement_frac))) + 
    	      geom_bar(stat='identity', position='dodge') +
	      ylab('average number of genes per region') + xlab('parameter set') +
	      scale_fill_viridis(discrete=TRUE) +
	      labs(fill='improvement fraction') +
	      scale_y_continuous(expand = c(0,0), limits=c(0,max(a$mean))) +
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(angle = 45,vjust = 1,hjust=1,colour="black"), 
	      axis.text.y = element_text(colour="black"))

ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/number_genes_in_region.pdf', width = 12, height = 7)


##======
# plot: number of introgressed bases for each strain
##======

print('plotting number of introgressed bases per strain')

	
a = read.table('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/plot_number_introgressed_bases_by_strain.txt', sep='\t')
names(a) = c('tag', 'improvement_frac', 'threshold', 'expected_length', 'expected_frac', 'strain', 'number_bases')

d = data.frame(tag = unique(a$tag), number_bases=NA, improvement_frac=NA)

#----
# for each tag, plot number of bases for each strain
#----

for (tag in tags) {
    print(tag)
    dir.create(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/',sep=''),showWarnings=FALSE)

    b = a[which(a$tag == tag),]
    d[which(d$tag == tag),]$number_bases = mean(b$number_bases)
    d[which(d$tag == tag),]$improvement_frac = b$improvement_frac[[1]]

    ggplot(b, aes(x=reorder(strain, -number_bases), y=number_bases, fill='x')) + 
    	      geom_bar(stat='identity',position='dodge') + 
	      xlab('strain') + ylab('number of bases introgressed') +
	      scale_fill_viridis(discrete=TRUE) +
	      guides(fill=FALSE) +
	      scale_y_continuous(expand = c(0,0), limits=c(0,max(b$number_bases))) +
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(angle = 45,vjust = 1,hjust=1,colour="black"), 
	      axis.text.y = element_text(colour="black"))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/number_bases_by_strain_',tag,'.pdf',sep=''), width = 12, height = 7)

}


#----
# plot average number of bases for all tags
#----

ggplot(d, aes(x=tag, y=number_bases, fill=as.factor(improvement_frac))) + 
    	      geom_bar(stat='identity',position='dodge') + 
	      xlab('tag') + ylab('number of bases introgressed') +
	      scale_y_continuous(expand = c(0,0), limits=c(0,max(d$number_bases))) +
	      scale_fill_viridis(discrete=TRUE) +
	      labs(fill='improvement fraction') +
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(angle = 45,vjust = 1,hjust=1,colour="black"), 
	      axis.text.y = element_text(colour="black"))
    
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/average_bases_per_strain.pdf', width = 12, height = 7)


##======
# plot: number of introgressed genes for each strain
##======

print('plotting number of introgressed genes per strain')
	
a = read.table('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/plot_number_introgressed_genes_by_strain.txt', sep='\t')
names(a) = c('tag', 'improvement_frac', 'threshold', 'expected_length', 'expected_frac', 'strain', 'number_genes')

d = data.frame(tag = unique(a$tag), number_genes=NA, improvement_frac=NA)

#----
# for each tag, plot number of genes for each strain
#----

for (tag in tags) {
    print(tag)
    dir.create(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/',sep=''),showWarnings=FALSE)

    b = a[which(a$tag == tag),]
    d[which(d$tag == tag),]$number_genes = mean(b$number_genes)
    d[which(d$tag == tag),]$improvement_frac = b$improvement_frac[[1]]
    ggplot(b, aes(x=reorder(strain, -number_genes), y=number_genes, fill='a')) + 
    	      geom_bar(stat='identity',position='dodge') + 
	      xlab('strain') + ylab('number of genes introgressed') +
	      scale_y_continuous(expand = c(0,0), limits=c(0,max(b$number_genes))) +
	      scale_fill_viridis(discrete=TRUE) +
	      guides(fill=FALSE) +
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(angle = 45,vjust = 1,hjust=1,colour="black"), 
	      axis.text.y = element_text(colour="black"))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/number_genes_by_strain_',tag,'.pdf',sep=''), width = 12, height = 7)

}


#----
# plot average number of genes for all tags
#----

ggplot(d, aes(x=tag, y=number_genes, fill=as.factor(d$improvement_frac))) + 
    	      geom_bar(stat='identity',position='dodge') + 
	      xlab('tag') + ylab('average number of bases introgressed per strain') +
	      scale_y_continuous(expand = c(0,0), limits=c(0,max(d$number_genes))) +
	      scale_fill_viridis(discrete=TRUE) +
	      labs(fill='improvement fraction') +
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(angle = 45,vjust = 1,hjust=1,colour="black"), 
	      axis.text.y = element_text(colour="black"))
    
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/average_genes_per_strain.pdf', width = 12, height = 7)


##======
# plot: number of strains each gene introgressed in 
##======

print('plotting number of strains per introgressed gene')

#----
# for each tag, plot number of strains for each gene
#----

for (tag in tags) {
    print(tag)
    dir.create(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/',sep=''),showWarnings=FALSE)

    a = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plot_number_strains_by_genes_',tag,'.txt', sep=''), sep='\t')
    names(a) = c('gene', 'number_strains')

    ggplot(a, aes(x=reorder(gene, -number_strains), y=number_strains, fill='a')) + 
    	      geom_bar(stat='identity',position='dodge') + 
	      xlab('gene') + ylab('number of strains introgressed in') +
	      scale_y_continuous(expand = c(0,0), limits=c(0,max(a$number_strains))) +
	      scale_fill_viridis(discrete=TRUE) +
	      guides(fill=FALSE) +
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(size=.5, angle = 45,vjust = 1,hjust=1,colour="black"), 
	      axis.text.y = element_text(colour="black"))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/number_strains_by_gene_',tag,'.pdf',sep=''), width = 12, height = 7)

}


#----
# plot average number of strains for all tags
#----

a = read.table('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/plot_number_strains_by_genes.txt', sep='\t')
names(a) = c('tag', 'improvement_frac', 'threshold', 'expected_length', 'expected_frac', 'mean', 'lower', 'upper', 'median', 'min', 'max')

ggplot(a, aes(x=tag, y=mean, fill=as.factor(a$improvement_frac))) + 
    	      geom_bar(stat='identity',position='dodge') + 
	      scale_y_continuous(expand = c(0,0), limits=c(0,max(a$mean))) +
	      scale_fill_viridis(discrete=TRUE) +
	      labs(fill='improvement fraction') +
	      xlab('tag') + ylab('average number of strains each gene introgressed in') +
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(angle = 45,vjust = 1,hjust=1,colour="black"), 
	      axis.text.y = element_text(colour="black"))
    
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/average_strains_per_gene.pdf', width = 12, height = 7)



##======
# plot: average fraction of each (introgressed) gene that's introgressed
##======

print('plotting average fraction of each introgressed gene that is introgressed')

#----
# for each tag, avg fraction vs. gene
#----

for (tag in tags) {
    print(tag)
    dir.create(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/',sep=''),showWarnings=FALSE)


    a = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plot_frac_introgressed_by_genes_',tag,'.txt',sep=''), sep='\t')
    names(a) = c('gene', 'avg_frac_introgressed', 'lower', 'upper', 'median', 'min', 'max')

    ggplot(a, aes(x=reorder(gene, -avg_frac_introgressed), y=avg_frac_introgressed, fill='a')) + 
    	      geom_bar(stat='identity',position='dodge') +
	      xlab('gene') + ylab('average fraction introgressed') +
	      scale_y_continuous(expand = c(0,0), limits=c(0,1)) +
	      scale_fill_viridis(discrete=TRUE) +
	      guides(fill=FALSE) +
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(size=.25, angle = 45,vjust = 1,hjust=1,colour="black"), 
	      axis.text.y = element_text(colour="black"))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/average_fraction_gene_introgressed_',tag,'.pdf',sep=''), width = 12, height = 7)

}
