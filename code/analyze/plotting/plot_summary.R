library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(dplyr)
source('../my_color_palette.R')

args = commandArgs(trailingOnly=TRUE)
tag = args[1]

strains3 = c('yjm1252', 'yjm1078', 'yjm248')

chrms = c('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI')

intd_refs = c("CBS432", "N_45", "DBVPG6304", "UWOPS91_917_1") 
intd_refs = c("par")

regions = data.frame()
for (ref in intd_refs) {
    #regions_ref = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/',
    #                               'introgression/results/analysis/', tag,
    #                               '/blocks_', ref, '_', tag, '_', suffix,
    #                               '.txt', sep=''), sep='\t', header=T,
    #                         stringsAsFactors=F)
    regions_ref = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/',
                                   'introgression/results/analysis/', tag,
                                   '/introgressed_blocks_filtered_', ref, '_',
                                   tag, '_summary_plus',
                                   '.txt', sep=''), sep='\t', header=T,
                             stringsAsFactors=F)
    print(nrow(regions_ref))
    regions_ref$variable = "variable" # make things a lil easier with dcast
    #regions[[ref]]$overlap_gene = regions$number_genes >= 1
    regions_ref$length = regions_ref$end - regions_ref$start + 1
    if (!"predicted_species" %in% names(regions_ref)) {
        regions_ref$predicted_species = ref
    }
    regions = rbind(regions, regions_ref)
}

# bar chart: genes introgressed per strain
ag = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/',
                      'introgression/results/analysis/', tag,
                      '/genes_for_each_strain_filtered_', tag,
                      '.txt', sep=''), sep='\t', header=T,
                stringsAsFactors=F)
ag2 = ag %>%
    group_by(strain) %>%
    summarise(total = sum(num_genes))
print(ag2)
ggplot(ag2, aes(x=reorder(toupper(strain), -total), y=total, fill='x')) + 
    geom_bar(stat='identity',position='dodge') + 
    xlab('Strain') + ylab('Number of genes introgressed') +
    scale_fill_manual(values =c(my_color_palette[['introgressed']]) ) +
    guides(fill=FALSE) +
    scale_y_continuous(expand = c(0,0), limits=c(0,max(ag2$total))) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.ticks=element_line(colour="black"),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(size=8,angle = 45,vjust = 1,hjust=1,colour="black"), 
          axis.text.y = element_text(colour="black"))
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/num_genes_bar','_',tag,'.pdf',sep=''), width = 12, height = 6)

## bar chart: fraction genome introgressed per strain
g = 12071326
a = dcast(regions, strain + predicted_species ~ variable, value.var="length", fun.aggregate=sum)
a$frac = a$variable/g
a$in3 = a$strain %in% strains3
print(a[1:10,])
ggplot(a, aes(x=reorder(toupper(strain), -frac), y=frac, fill=in3)) + 
    geom_bar(stat='identity',position='dodge') + 
    xlab('Strain') + ylab('Fraction of genome introgressed') +
    scale_fill_manual(values =c(my_color_palette[['introgressed']], "#ABBEAB") ) +
    guides(fill=FALSE) +
    scale_y_continuous(expand = c(0,0), limits=c(0,max(a$frac))) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.ticks=element_line(colour="black"),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(size=5, angle = 45,
                                     vjust = 1, hjust=1, colour="black"), 
          axis.text.y = element_text(colour="black"))
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/frac_bar','_',tag,'.pdf',sep=''), width = 10, height = 5)
print(mean(a$frac))
print(median(a$frac))
print(min(a$frac))
print(max(a$frac))

## histogram: region lengths, 3 strains

r3 = regions[which(regions$strain %in% strains3),]
ggplot(r3, aes(x=length, fill=strain)) + geom_histogram(binwidth=100) +
    xlab('Region length') + ylab('Number of regions') +
    #scale_fill_manual(values =c(my_color_palette[['introgressed']]) ) +
    guides(fill=FALSE) +
    #scale_y_continuous(expand = c(0,0), limits=c(0,1500)) + 
    scale_x_continuous(expand = c(0,0), limits=c(0,max(r3$length)*1.05)) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.line=element_line(),
          axis.ticks=element_line(colour="black"),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(colour="black",size=12), 
          axis.text.y = element_text(colour="black",size=12))
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/length_histr3','_',tag,'.pdf',sep=''), width = 9, height = 6)
print(mean(r3$length))
print(median(r3$length))
print(min(r3$length))
print(max(r3$length))
print(length(r3$length))

## histogram: region lengths, all but 3 strains

nr3 = regions[which(!regions$strain %in% strains3),]
ggplot(nr3, aes(x=length, fill=predicted_species)) + geom_histogram(binwidth=100) +
    xlab('Region length') +
    ylab('Number of regions') +
    scale_fill_manual(values =c(my_color_palette[['introgressed']]) ) +
    guides(fill=FALSE) +
    #scale_y_continuous(expand = c(0,0), limits=c(0,)) + 
    scale_x_continuous(expand = c(0,0), limits=c(0,max(nr3$length)*1.05)) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.line=element_line(),
          axis.ticks=element_line(colour="black"),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=18),
          axis.text.x = element_text(colour="black",size=12),
          axis.text.y = element_text(colour="black",size=12))
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/length_hist_notnr3','_',tag,'.pdf',sep=''), width = 9, height = 6)
print(mean(nr3$length))
print(median(nr3$length))
print(min(nr3$length))
print(max(nr3$length))
print(length(nr3$length))

# boxplot: region lengths for only 3 strains and for rest of strains
regions$in3 = regions$strain %in% strains3
regions[which(regions$in3 == TRUE),]$in3 = "Three highly-\nintrogressed strains"
regions[which(regions$in3 == FALSE),]$in3 = "Ninety other\nstrains"
regions$in3 = factor(regions$in3, levels = c("Three highly-\nintrogressed strains",
                                             "Ninety other\nstrains"))
ggplot(regions, aes(x=in3, y=length/1000, fill = "a")) +
    geom_boxplot() +
    ylab('Region length (kb)') +
    xlab('') +
    scale_fill_manual(values = c(my_color_palette[['introgressed']]) ) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray90"),
          panel.grid.major=element_line(colour="gray80"),
          legend.position = "none",
          axis.line=element_line(),
          axis.ticks=element_line(colour="black"),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=18),
          axis.text.x = element_text(colour="black",size=12),
          axis.text.y = element_text(colour="black",size=12))

ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/length_boxplot_r3','_',tag,'.pdf',sep=''), width = 4, height = 6)

# ecdf: region lengths for only 3 strains and for rest of strains
regions$in3 = regions$strain %in% strains3
regions[which(regions$in3 == TRUE),]$in3 = "Three highly-introgressed strains"
regions[which(regions$in3 == FALSE),]$in3 = "Ninety other strains"
regions$in3 = factor(regions$in3, levels = c("Three highly-introgressed strains",
                                             "Ninety other strains"))
ggplot(regions, aes(x=length/1000, colour=in3)) +
    stat_ecdf(size=2, alpha = .8) +
    xlab('Region length (kb)') +
    ylab('Fraction of regions shorter') +
    scale_x_continuous(expand=c(0,.4)) +
    scale_y_continuous(expand=c(0,.01)) +
    scale_colour_manual(values = c("#ABBEAB", my_color_palette[['introgressed']])) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray90"),
          panel.grid.major=element_line(colour="gray80"),
          legend.title=element_blank(),
          legend.position = c(.6,.1),
          legend.key = element_rect(fill = "transparent"),
          legend.text=element_text(size=18),
          axis.line=element_line(),
          axis.ticks=element_line(colour="black"),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          axis.text.x = element_text(colour="black",size=14),
          axis.text.y = element_text(colour="black",size=14))

ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/length_ecdf_r3','_',tag,'.png',sep=''), width = 7, height = 7)


####################
stop here

# histogram: region lengths
ggplot(regions, aes(x=length, fill=predicted_species)) + geom_histogram(binwidth=100) +
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
ggplot(regions, aes(x=length, fill=predicted_species)) + geom_histogram(binwidth=100) +
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
ggplot(regions, aes(x=length, fill=predicted_species)) + geom_histogram(binwidth=50) +
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

asdf

#####################################

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
