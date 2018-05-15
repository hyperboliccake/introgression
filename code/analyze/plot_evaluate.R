# this is going to make a bunch of plots to help evaluate which
# regions look real and which should be filtered out

# in particular, it should be useful to look at regions overlapping
# and not overlapping genes; the ones overlapping genes should on
# average be more accurate (or at least less often due to poor
# alignment)

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)


args = commandArgs(trailingOnly=TRUE)
tag = args[1]
suffix = ''
if (length(args) == 2) {
    suffix = args[2]
}


regions = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/introgressed_blocks', suffix,'_par_', tag, '_summary_plus.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)
regions$overlap_gene = regions$number_genes >= 1
regions$length = regions$end - regions$start + 1
regions$fraction_gap = regions$number_gaps / regions$aligned_length
regions$fraction_gap_masked = (regions$number_gaps + regions$number_masked_non_gap) / regions$aligned_length
regions$cer_id = regions$number_match_ref1 / (regions$aligned_length - regions$number_gaps)
regions$par_id = regions$number_match_ref2 / (regions$aligned_length - regions$number_gaps)


#regions_filtered = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/introgressed_blocks_filtered_par_', tag, '_summary_plus.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)

##=====
# par id vs cer id
##=====

ggplot(regions, (aes(x=cer_id, y=par_id, label=region_id))) + geom_point(size=.2, alpha=.5) +  coord_cartesian(xlim=c(.6, 1),ylim=c(.6,1)) + geom_abline(slope=1,intercept=0)
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/par_id_vs_cer_id_',tag,'.pdf',sep=''), width = 8, height = 8)

ggplot(regions, (aes(x=cer_id, y=par_id, label=region_id))) + geom_point(size=.2, alpha=.5) +  coord_cartesian(xlim=c(.6, 1),ylim=c(.6,1)) + geom_abline(slope=1,intercept=0) + geom_text(aes(label=as.character(region_id)),hjust=0,vjust=0, cex=.2)
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/par_id_vs_cer_id_labeled_',tag,'.pdf',sep=''), width = 8, height = 8)

asdg

##=====
# comparing patterns in regions different distances from telomeres
##=====

# fraction gaps vs distance from telomere
ggplot(regions, (aes(x=distance_from_telomere, y=fraction_gap, label=region_id))) + geom_point(size=.2, alpha=.5) + coord_cartesian(xlim=c(0,10000)) #+ geom_text(aes(label=ifelse(length>1000,as.character(region_id),'')),hjust=0,vjust=0, cex=.2)
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/frac_gaps_vs_dist_from_tel_',tag,'.pdf',sep=''), width = 12, height = 7)

# fraction gaps vs distance from centromere
ggplot(regions, (aes(x=distance_from_centromere, y=fraction_gap, label=region_id))) + geom_point(size=.2, alpha=.5) #+ geom_text(aes(label=ifelse(length>1000,as.character(region_id),'')),hjust=0,vjust=0, cex=.2)
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/frac_gaps_vs_dist_from_cen_',tag,'.pdf',sep=''), width = 12, height = 7)

# longest gap vs distance from telomere
ggplot(regions, (aes(x=distance_from_telomere, y=longest_gap, label=region_id))) + geom_point(size=.2, alpha=.5) #+ geom_text(aes(label=ifelse(length>1000,as.character(region_id),'')),hjust=0,vjust=0, cex=.2)
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/longest_gap_vs_dist_from_tel_',tag,'.pdf',sep=''), width = 12, height = 7)

##=====
# comparing patterns in regions that overlap and don't overlap genes
##=====

# scatter of frac gaps vs length
ggplot(regions, (aes(x=length, y=fraction_gap, colour=overlap_gene, label=region_id))) + geom_point(size=.2, alpha=.5) + geom_text(aes(label=ifelse(length>1000,as.character(region_id),'')),hjust=0,vjust=0, cex=.2)
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/gaps_vs_length_by_overlap_gene_',tag,'.pdf',sep=''), width = 12, height = 7)

## scatter of frac gaps+masked vs length
ggplot(regions, (aes(x=length, y=fraction_gap_masked, colour=overlap_gene, label=region_id))) + geom_point(size=.2, alpha=.5) + geom_text(aes(label=ifelse(length>1000,as.character(region_id),'')),hjust=0,vjust=0, cex=.2)
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/gaps_masked_vs_length_by_overlap_gene_',tag,'.pdf',sep=''), width = 12, height = 7)

# above with (no labels)
ggplot(regions, (aes(x=length, y=fraction_gap_masked, colour=overlap_gene))) + geom_point(size=.2, alpha=.5)
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/gaps_masked_vs_length_by_overlap_gene_nolab_',tag,'.pdf',sep=''), width = 12, height = 7)

## scatter of frac gaps+masked vs cer_id
ggplot(regions, (aes(x=cer_id, y=fraction_gap_masked, colour=overlap_gene, label=region_id))) + geom_point(size=.2, alpha=.5) + geom_text(aes(label=ifelse(length>1000,as.character(region_id),'')),hjust=0,vjust=0, cex=.2)
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/gaps_masked_vs_cer_id_by_overlap_gene_',tag,'.pdf',sep=''), width = 12, height = 7)

## scatter of frac gaps+masked vs cer_id
ggplot(regions, (aes(x=cer_id, y=fraction_gap_masked, colour=overlap_gene))) + geom_point(size=.2, alpha=.5)
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/gaps_masked_vs_cer_id_by_overlap_gene_nolab_',tag,'.pdf',sep=''), width = 12, height = 7)


## lengths
ggplot(regions, aes(x=as.factor(overlap_gene), y=(length))) + geom_violin()
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/violin_region_length_vs_number_genes_',tag,'.pdf',sep=''), width = 12, height = 7)

ggplot(regions, aes(x=as.factor(overlap_gene), y=number_non_gap)) + geom_violin()
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/violin_nongap_region_length_vs_number_genes_',tag,'.pdf',sep=''), width = 12, height = 7)

# number of gaps
ggplot(regions, aes(x=as.factor(overlap_gene), y=(number_gaps))) + geom_violin()
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/violin_num_gaps_vs_number_genes_',tag,'.pdf',sep=''), width = 12, height = 7)

# fraction of gaps
ggplot(regions, aes(x=as.factor(overlap_gene), y=(fraction_gap))) + geom_violin()
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/violin_frac_gaps_vs_number_genes_',tag,'.pdf',sep=''), width = 12, height = 7)

# longest gap stretch
ggplot(regions, aes(x=as.factor(overlap_gene), y=(longest_gap))) + geom_violin()
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/violin_longest_gap_vs_number_genes_',tag,'.pdf',sep=''), width = 12, height = 7)

# fraction gap+masked
ggplot(regions, aes(x=as.factor(overlap_gene), y=(fraction_gap_masked))) + geom_violin()
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/violin_gap_masked_vs_number_genes_',tag,'.pdf',sep=''), width = 12, height = 7)



# fraction gap+masked histogram
ggplot(regions, aes(x=(fraction_gap_masked))) + geom_histogram(binwidth=.01)
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/gap_masked_hist_',tag,'.pdf',sep=''), width = 12, height = 7)

# fraction gap+masked cdf
ggplot(regions, aes(x=(fraction_gap_masked))) + stat_ecdf()
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/gap_masked_cdf_',tag,'.pdf',sep=''), width = 12, height = 7)




# region length histogram

# non gap length histogram

# fraction gaps histogram

# fraction gaps vs total length

# fraction gaps vs overlap/not overlap gene

# num sites match only par vs fraction gaps

# num sites match only par vs length

# fraction sites match only par vs length

# num sites match only par vs non gap length

# fraction sites match only par vs non gap length

# num sites match only par histogram

# fraction sites match only par histogram

# num sites match only par vs overlap/not overlap gene

# ... vs no gene/gene with paralog/gene with no paralog

# ... vs distance from telomere
