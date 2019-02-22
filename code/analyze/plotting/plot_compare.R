# plotting things to help evaluate the regions I found vs the regions
# Strope et al found

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(VennDiagram)
source('../my_color_palette.R')


args = commandArgs(trailingOnly=TRUE)
tag = args[1]
suffix = ''
if (length(args) == 2) {
    suffix = args[2]
}

regions = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/introgressed_blocks', suffix,'_par_', tag,
                           '_summary_plus.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)
#regions$overlap_gene = regions$number_genes >= 1
#regions$length = regions$end - regions$start + 1
#regions$fraction_gap = regions$number_gaps / regions$aligned_length
#regions$fraction_gap_masked = (regions$number_gaps + regions$number_masked_non_gap) / regions$aligned_length
#regions$cer_id = regions$number_match_ref1 / (regions$aligned_length - regions$number_gaps)


genes_strope = read.table('compare_to_strope/genes_both.txt', header=F, stringsAsFactors=F)

gt = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/gene_summary', suffix,'_', tag, '.txt', sep=''),
                sep='\t', header=T, stringsAsFactors=F)

##=====
# venn diagram
##=====

x = nrow(read.table('compare_to_strope/genes_me_only.txt',
                    header=F, stringsAsFactors=F))
y = nrow(read.table('compare_to_strope/genes_strope_only.txt',
                    header=F, stringsAsFactors=F))
xy = nrow(read.table('compare_to_strope/genes_both.txt',
                     header=F, stringsAsFactors=F))
x = x + xy
y = y + xy

v = draw.pairwise.venn(x, y, xy,
                       #category = c('Our predictions', 'Strope et al. predictions'),
                       #cat.col = c(my_color_palette[['misc1']],
                       #            my_color_palette[['misc2']]),
                       fill=c(my_color_palette[['misc1']],
                              my_color_palette[['misc2']]),
                       alpha=.7, cex=4)
png(filename = paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/',
                     'results/analysis/', tag,
                     '/plots/compare_to_strope_venn.png', sep=''),
    width = 3, height = 3, units = 'in', res=300)
grid.draw(v)
dev.off()


##=====
# average fraction of gene introgressed vs avg containing region length,
# labeled by whether or not called by strope et al
##=====

gt$strope_found = gt$gene %in% genes_strope[,1]

ggplot(gt, (aes(x=avg_region_length, y=avg_frac_intd, colour=strope_found, label=gene))) + geom_point(size=1, alpha=.5)+
    xlab('average region length') +
    ylab('average fraction of gene introgressed') + 
        theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          legend.title = element_blank(),
          legend.text = element_text(size=18), 
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"))

ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/intd_vs_region_length_by_strope_found_',tag,'.pdf',sep=''), width = 12, height = 7)

ggplot(gt, (aes(x=avg_region_length, y=avg_frac_intd, colour=strope_found, label=gene))) + geom_point(size=1, alpha=.5) + geom_text(aes(label=gene),hjust=0,vjust=0, cex=.6)+
        theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"))

ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/intd_vs_region_length_by_strope_found_labeled_',tag,'.pdf',sep=''), width = 12, height = 7)


##=====
# gene identity with cer reference vs avg containing region length,
# labeled by whether or not called by strope et al
##=====

ggplot(gt, (aes(x=avg_region_length, y=average_cer_ref_id, colour=strope_found, label=gene))) + geom_point(size=1, alpha=.5) + geom_text(aes(label=gene),hjust=0,vjust=0, cex=.6) +
    xlab('average region length') +
    ylab('average introgressed gene identity to cerevisiae reference') + 
    scale_color_hue(labels = c("not found by Strope et al.", "found by Strope et al.")) +
        theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          legend.title = element_blank(),
          legend.text = element_text(size=18), 
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"))


ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/cer_ref_id_vs_region_length_by_strope_found_',tag,'.pdf',sep=''), width = 12, height = 7)

ggplot(gt, (aes(x=avg_region_length, y=average_cer_ref_id, colour=strope_found, label=gene))) + geom_point(size=1, alpha=.5) +
            theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"))

ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/cer_ref_id_vs_region_length_by_strope_found_labeled_',tag,'.pdf',sep=''), width = 12, height = 7)



##=====
# gene identity with cer reference vs avg frac introgressed,
# labeled by whether or not called by strope et al
##=====

ggplot(gt, (aes(x=avg_frac_intd, y=average_cer_ref_id,
                colour=strope_found, shape=strope_found))) +
    geom_point(size=3, alpha=.8) +
    ylab('Average identity with cerevisiae reference') +
    xlab('Average fraction of gene we call introgressed') +
    labs(colour="Also found by\nStrope et al.?",
         shape="Also found by\nStrope et al.?") +
    scale_color_manual(values = c(my_color_palette[['misc1']],
                                  my_color_palette[['misc2']]),
                       labels = c("no", "yes")) +
    scale_shape_manual(values=c(19,17),labels=c("no","yes")) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"),
          panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          legend.title = element_text(size=18),
          legend.text = element_text(size=16),
          legend.position = c(.85, .15),
          legend.key = element_rect(fill = "transparent"),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(colour="black",size=14), 
          axis.text.y = element_text(colour="black",size=14))

ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/cer_ref_id_vs_intd_by_strope_found_',tag,'.pdf',sep=''), width = 8.5, height = 7)

stop here

ggplot(gt, (aes(x=avg_frac_intd, y=average_cer_ref_id, colour=strope_found, label=gene))) + geom_point(size=1, alpha=.5) + geom_text(aes(label=gene),hjust=0,vjust=0, cex=.6) +
    ylab('average identity with cerevisiae reference') + xlab('average fraction of gene I call introgressed') + 
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"))


ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/cer_ref_id_vs_intd_by_strope_found_labeled_',tag,'.pdf',sep=''), width = 12, height = 7)

sagdasg


##=====
# not really for comparison but uses the same data: histrogram of average fraciton of gene introgressed
##=====

bw = .05
ggplot(gt, (aes(x=avg_frac_intd, fill='x', label=gene))) + geom_histogram(binwidth=.05, boundary=0) +
    ylab('number of genes') + xlab('average fraction of gene introgressed') + 
    scale_fill_viridis(discrete=TRUE) +
    guides(fill=FALSE) +
     scale_y_continuous(expand = c(0,0), limits=c(0,200)) + 
    scale_x_continuous(expand = c(0,0), limits=c(0,1)) +
           theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"))

print(min(gt$avg_frac_intd))
print(max(gt$avg_frac_intd))
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/frac_gene_intd_hist_',tag,'.pdf',sep=''), width = 12, height = 7)



##=====
# par id vs cer id,
# labeled by whether or not called by strope et al
##=====

#ggplot(gt, (aes(x=cer_id, y=par_id, label=region_id, colour=strope_found))) + geom_point(size=.2, alpha=.5) +  coord_cartesian(xlim=c(.6, 1),ylim=c(.6,1)) + geom_abline(slope=1,intercept=0) + geom_text(aes(label=as.character(region_id)),hjust=0,vjust=0, cex=.2)
#ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/par_id_vs_cer_id_by_strope_found',tag,'.pdf',sep=''), width = 8, height = 8)
