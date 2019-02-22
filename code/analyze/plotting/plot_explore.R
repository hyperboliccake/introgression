library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)

genome_size = 12071326

predict_args = read.table('predict_args.txt', sep=' ', stringsAsFactors=F)
predict_args = predict_args[c(seq(1,19),seq(21,25),seq(27,36)),]
tags = predict_args[,1]

genes_by_region = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/plot_number_genes_by_region.txt',sep=''), sep='\t')
names(genes_by_region) = c('tag', 'improvement_frac', 'threshold', 'expected_length', 'expected_frac', 'mean', 'bs_lower', 'bs_upper', 'median', 'min', 'max')

region_lengths = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/plot_region_lengths.txt',sep=''), sep='\t')
names(region_lengths) = c('tag', 'improvement_frac', 'threshold', 'expected_length', 'expected_frac', 'mean', 'bs_lower', 'bs_upper', 'median', 'min', 'max', 'total_num_regions')

genes_by_strain_long = read.table('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/plot_number_introgressed_genes_by_strain.txt', sep='\t')
names(genes_by_strain_long) = c('tag', 'improvement_frac', 'threshold', 'expected_length', 'expected_frac', 'strain', 'number_genes')
genes_by_strain = dcast(genes_by_strain_long, tag + improvement_frac + threshold + expected_length + expected_frac ~ 'number_genes', value.var='number_genes', fun.aggregate=mean)

strains_by_gene = read.table('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/plot_number_strains_by_genes.txt', sep='\t')
names(strains_by_gene) = c('tag', 'improvement_frac', 'threshold', 'expected_length', 'expected_frac', 'mean', 'lower', 'upper', 'median', 'min', 'max', 'total_num_genes', 'total_num_genes_1', 'total_num_genes_g1')

bases_by_strain_long = read.table('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/plot_number_introgressed_bases_by_strain.txt', sep='\t')
names(bases_by_strain_long) = c('tag', 'improvement_frac', 'threshold', 'expected_length', 'expected_frac', 'strain', 'number_bases')
bases_by_strain = dcast(bases_by_strain_long, tag + improvement_frac + threshold + expected_length + expected_frac ~ 'number_bases', value.var='number_bases', fun.aggregate=mean)
bases_by_strain_frac_diffs = dcast(bases_by_strain_long, improvement_frac + threshold + expected_length + strain ~ expected_frac, value.var='number_bases')
bases_by_strain_frac_diffs$diff = bases_by_strain_frac_diffs$'0.1' - bases_by_strain_frac_diffs$'0.001'

##======
# does assumed length/fraction change results (hopefully not)
##======

print('')

#----
# for different lengths or fractions (connected by line, holding other parameters constant), plot:
# * avg bp introgressed
# * avg introgressed region size
# * total number of introgressed regions
# * avg number of genes introgressed/strain
# * total number of genes introgressed
#----

ggplot(bases_by_strain_frac_diffs, aes(diff/genome_size, fill=interaction(improvement_frac, threshold, expected_length))) + 
	      geom_histogram(position='identity', alpha=.2, binwidth = .001) +
	      xlab('difference in fraction of genome introgressed') + ylab('number of strains') +
	      scale_colour_viridis(discrete=TRUE) +
	      guides(colour=FALSE) +
	      #scale_y_continuous(limit=c(0,93)) + 
	      #scale_x_continuous(limit=c(0,1)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/bases_per_strain_vs_expected_length_hist.pdf', width=12, height=7)

sgadg

ggplot(bases_by_strain, aes(x=expected_length, y=number_bases/genome_size, group=interaction(improvement_frac, threshold, expected_frac), shape=as.factor(improvement_frac), colour=as.factor(threshold), linetype=as.factor(expected_frac))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('average fraction of genome introgressed') + xlab('expected introgressed region size') +
	      scale_colour_viridis(discrete=TRUE) +
	      guides(colour=FALSE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$expected_length)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/bases_per_strain_vs_expected_length.pdf', width=12, height=7)

ggplot(bases_by_strain, aes(x=expected_frac, y=number_bases/genome_size, group=interaction(improvement_frac, expected_length, threshold,strain), shape=as.factor(improvement_frac), linetype=as.factor(expected_length), colour=as.factor(threshold))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('average fraction of genome introgressed') + xlab('expected fraction of genome introgressed') +
	      scale_colour_viridis(discrete=TRUE) +
	      guides(colour=FALSE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$expected_fraction)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/bases_per_strain_vs_expected_fraction.pdf', width=12, height=7)

ggplot(region_lengths, aes(x=expected_length, y=total_num_regions, group=interaction(threshold, improvement_frac, expected_frac), colour=as.factor(threshold), shape=as.factor(improvement_frac), linetype=as.factor(expected_frac))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('number of introgressed regions') + xlab('expected introgressed region size') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$expected_length)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/num_regions_vs_expected_length.pdf', width=12, height=7)

ggplot(region_lengths, aes(x=expected_frac, y=total_num_regions, group=interaction(threshold, improvement_frac, expected_length), colour=as.factor(threshold), shape=as.factor(improvement_frac), linetype=as.factor(expected_length))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('number of introgressed regions') + xlab('expected fraction of genome introgressed') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$expected_frac)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/num_regions_vs_expected_frac.pdf', width=12, height=7)

ggplot(region_lengths, aes(x=expected_length, y=mean, group=interaction(threshold, improvement_frac, expected_frac), colour=as.factor(threshold), shape=as.factor(improvement_frac), linetype=as.factor(expected_frac))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('average introgressed region size') + xlab('expected introgressed region size') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$expected_length)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/average_region_size_vs_expected_length.pdf', width=12, height=7)

ggplot(region_lengths, aes(x=expected_frac, y=mean, group=interaction(threshold, improvement_frac, expected_length), colour=as.factor(threshold), shape=as.factor(improvement_frac), linetype=as.factor(expected_length))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('average introgressed region size') + xlab('expected fraction of genome introgressed') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$expected_frac)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/average_region_size_vs_expected_frac.pdf', width=12, height=7)

ggplot(genes_by_strain, aes(x=expected_length, y=number_genes, group=interaction(threshold, improvement_frac, expected_frac), colour=as.factor(threshold), shape=as.factor(improvement_frac), linetype=as.factor(expected_frac))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('number of genes per strain') + xlab('expected introgressed region size') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$expected_length)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/genes_per_strain_vs_expected_length.pdf', width=12, height=7)

ggplot(genes_by_strain, aes(x=expected_frac, y=number_genes, group=interaction(threshold, improvement_frac, expected_length), colour=as.factor(threshold), shape=as.factor(improvement_frac), linetype=as.factor(expected_length))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('number of genes per strain') + xlab('expected fraction of genome introgressed') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$expected_frac)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/genes_per_strain_vs_expected_frac.pdf', width=12, height=7)

ggplot(strains_by_gene, aes(x=expected_length, y=total_num_genes, group=interaction(threshold, improvement_frac, expected_frac), colour=as.factor(threshold), shape=as.factor(improvement_frac), linetype=as.factor(expected_frac))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('total number of introgressed genes') + xlab('expected introgressed region size') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$expected_length)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/num_genes_vs_expected_length.pdf', width=12, height=7)

ggplot(strains_by_gene, aes(x=expected_frac, y=total_num_genes, group=interaction(threshold, improvement_frac, expected_length), colour=as.factor(threshold), shape=as.factor(improvement_frac), linetype=as.factor(expected_length))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('total number of introgressed genes') + xlab('expected fraction of genome introgressed') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$expected_frac)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/num_genes_vs_expected_frac.pdf', width=12, height=7)

##======
# what happens when we train more
##======

print('')

#----
# for different likelihood increase thresholds
# * avg bp introgressed
# * total number of introgressed regions
# * avg introgressed region size 
# * median introgressed region size
# * avg number of genes/region
# * total number of genes introgressed
#----

ggplot(bases_by_strain, aes(x=improvement_frac, y=number_bases, group=interaction(threshold, expected_length, expected_frac), colour=as.factor(threshold), linetype=as.factor(expected_length), shape=as.factor(expected_frac))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('number of introgressed bases per strain') + xlab('training threshold') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$improvement_frac)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/bases_per_strain_vs_improvement_frac.pdf', width=12, height=7)

ggplot(region_lengths, aes(x=improvement_frac, y=total_num_regions, group=interaction(threshold, expected_length, expected_frac), colour=as.factor(threshold), linetype=as.factor(expected_length), shape=as.factor(expected_frac))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('total number of introgressed regions') + xlab('training threshold') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$improvement_frac)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/num_regions_vs_improvement_frac.pdf', width=12, height=7)

ggplot(region_lengths, aes(x=improvement_frac, y=mean, group=interaction(threshold, expected_length, expected_frac), colour=as.factor(threshold), linetype=as.factor(expected_length), shape=as.factor(expected_frac))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('average introgressed region size') + xlab('training threshold') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$improvement_frac)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/average_region_size_vs_improvement_frac.pdf', width=12, height=7)

ggplot(region_lengths, aes(x=improvement_frac, y=median, group=interaction(threshold, expected_length, expected_frac), colour=as.factor(threshold), linetype=as.factor(expected_length), shape=as.factor(expected_frac))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('median introgressed region size') + xlab('training threshold') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$improvement_frac)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/median_region_size_vs_improvement_frac.pdf', width=12, height=7)


ggplot(genes_by_region, aes(x=improvement_frac, y=mean, group=interaction(threshold, expected_length, expected_frac), colour=as.factor(threshold), linetype=as.factor(expected_length), shape=as.factor(expected_frac))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('average number of genes per introgressed region') + xlab('training threshold') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(genes_by_region$improvement_frac)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/genes_per_region_vs_improvement_frac.pdf', width=12, height=7)

ggplot(strains_by_gene, aes(x=improvement_frac, y=total_num_genes, group=interaction(threshold, expected_length, expected_frac), colour=as.factor(threshold), linetype=as.factor(expected_length), shape=as.factor(expected_frac))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('total number of introgressed genes') + xlab('training threshold') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(genes_by_region$improvement_frac)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/num_genes_vs_improvement_frac.pdf', width=12, height=7)



##======
# what happens when change posterior threshold
##======

print('')

#----
# for different likelihood increase thresholds
# * avg bp introgressed
# * total number of introgressed regions
# * avg introgressed region size
# * avg median region size
# * avg number of genes/region
# * total number of genes introgressed
# * number of genes introgressed in just one strain
#----

ggplot(bases_by_strain, aes(x=threshold, y=number_bases, group=interaction(improvement_frac, expected_length, expected_frac), colour=as.factor(improvement_frac), linetype=as.factor(expected_length), shape=as.factor(expected_frac))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('number of introgressed bases per strain') + xlab('posterior threshold') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$threshold)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/bases_per_strain_vs_threshold.pdf', width=12, height=7)

ggplot(region_lengths, aes(x=threshold, y=total_num_regions, group=interaction(improvement_frac, expected_length, expected_frac), colour=as.factor(improvement_frac), linetype=as.factor(expected_length), shape=as.factor(expected_frac))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('total number of introgressed regions') + xlab('posterior threshold') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$threshold)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/num_regions_vs_threshold.pdf', width=12, height=7)

ggplot(region_lengths, aes(x=threshold, y=mean, group=interaction(improvement_frac, expected_length, expected_frac), colour=as.factor(improvement_frac), linetype=as.factor(expected_length), shape=as.factor(expected_frac))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('average introgressed region size') + xlab('posterior threshold') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$threshold)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/average_region_size_vs_threshold.pdf', width=12, height=7)

ggplot(region_lengths, aes(x=threshold, y=median, group=interaction(improvement_frac, expected_length, expected_frac), colour=as.factor(improvement_frac), linetype=as.factor(expected_length), shape=as.factor(expected_frac))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('median introgressed region size') + xlab('posterior threshold') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$threshold)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/median_region_size_vs_threshold.pdf', width=12, height=7)

ggplot(genes_by_region, aes(x=threshold, y=mean, group=interaction(improvement_frac, expected_length, expected_frac), colour=as.factor(improvement_frac), linetype=as.factor(expected_length), shape=as.factor(expected_frac))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('average number introgressed genes per region') + xlab('posterior threshold') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$threshold)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/genes_per_region_vs_threshold.pdf', width=12, height=7)

ggplot(strains_by_gene, aes(x=threshold, y=total_num_genes, group=interaction(improvement_frac, expected_length, expected_frac), colour=as.factor(improvement_frac), linetype=as.factor(expected_length), shape=as.factor(expected_frac))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('total number of introgressed genes') + xlab('posterior threshold') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$threshold)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/num_genes_vs_threshold.pdf', width=12, height=7)

ggplot(strains_by_gene, aes(x=threshold, y=total_num_genes_1, group=interaction(improvement_frac, expected_length, expected_frac), colour=as.factor(improvement_frac), linetype=as.factor(expected_length), shape=as.factor(expected_frac))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('total number of genes introgressed in only one strain') + xlab('posterior threshold') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$threshold)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/num_genes_1_vs_threshold.pdf', width=12, height=7)

ggplot(strains_by_gene, aes(x=threshold, y=total_num_genes_g1, group=interaction(improvement_frac, expected_length, expected_frac), colour=as.factor(improvement_frac), linetype=as.factor(expected_length), shape=as.factor(expected_frac))) + 
    	      geom_point() +
	      geom_line() +
	      ylab('total number of genes introgressed in more than one strain') + xlab('posterior threshold') +
	      scale_colour_viridis(discrete=TRUE) +
	      scale_x_continuous(breaks=unique(bases_by_strain$threshold)) + 
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/num_genes_g1_vs_threshold.pdf', width=12, height=7)
















zbcxzcb








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

}

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
