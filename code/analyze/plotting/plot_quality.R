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
library(hexbin)
library(data.table)
source('../my_color_palette.R')

args = commandArgs(trailingOnly=TRUE)
tag = args[1]
refs = c("S288c", "CBS432", "N_45", "DBVPG6304", "UWOPS91_917_1")

regions = list()
regions_f1 = list()
regions_f2 = list()
for (ref in refs[2:5]) {
    regions[[ref]] = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/blocks_', ref, '_', tag, '_quality.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)

    regions_f1[[ref]] = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/blocks_', ref, '_', tag, '_filtered1.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)

    regions_f2[[ref]] = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/blocks_', ref, '_', tag, '_filtered2.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)
}

# no quality file for unknown regions
regions_unk = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/blocks_unknown_', tag, '.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)

##=====
# predicted species id vs cer id, various kinds of filtering
##=====

for (ref in refs[2:5]) {

    ######
    ## plot identity in terms of HMM sites
    ######
    
    for (ref2 in refs[1:5]) {
        regions[[ref]][[paste('frac_hmm_', ref2, sep='')]] =
            regions[[ref]][[paste("match_hmm_", ref2, sep="")]] /
            regions[[ref]][["num_sites_hmm"]]
        regions[[ref]][[paste('frac_nongap_', ref2, sep='')]] =
            regions[[ref]][[paste("match_nongap_", ref2, sep="")]] /
            regions[[ref]][[paste("num_sites_nongap_", ref2, sep='')]]

    }

    # unlabeled
    
    ggplot(regions[[ref]], aes_string(x=paste('frac_hmm_', refs[1], sep=''), y=paste('frac_hmm_', ref, sep=''))) +
        geom_point(alpha=.15) +
        geom_abline(slope=1,intercept=0, linetype='dashed') +
        xlab(paste('Identity with ', refs[1], ' reference', sep='')) + 
        ylab(paste('Identity with ', ref, ' reference', sep='')) + 
        scale_colour_manual(values = c(my_color_palette[['introgressed']], my_color_palette[['nonintrogressed']])) +
        theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
              axis.line=element_line(),
              #legend.position = "none",
              axis.title.x = element_text(size=18), 
              axis.title.y = element_text(size=18), 
              axis.text.x = element_text(colour="black"), 
              axis.text.y = element_text(colour="black"))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', ref, '_hmm_id_vs_', refs[1], '_hmm_id_',tag,'.png', sep=''),
           width = 8, height = 8, dpi=300)

    # labeled

    ggplot(regions[[ref]], aes_string(x=paste('frac_hmm_', refs[1], sep=''), y=paste('frac_hmm_', ref, sep=''))) +
        geom_text(aes(label=as.character(region_id)),hjust=0,vjust=0, cex=1) +
        geom_point(alpha=.15) +
        geom_abline(slope=1,intercept=0, linetype='dashed') +
        xlab(paste('Identity with ', refs[1], ' reference', sep='')) + 
        ylab(paste('Identity with ', ref, ' reference', sep='')) + 
        scale_colour_manual(values = c(my_color_palette[['introgressed']], my_color_palette[['nonintrogressed']])) +
        theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
              axis.line=element_line(),
              #legend.position = "none",
              axis.title.x = element_text(size=18), 
              axis.title.y = element_text(size=18), 
              axis.text.x = element_text(colour="black"), 
              axis.text.y = element_text(colour="black"))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', ref, '_hmm_id_vs_', refs[1], '_hmm_id_labeled_',tag,'.png', sep=''),
           width = 8, height = 8, dpi=300)

    ######
    ## plot identity in terms of all nongap sites
    ######
    
    # unlabeled
    
    ggplot(regions[[ref]], aes_string(x=paste('frac_nongap_', refs[1], sep=''), y=paste('frac_nongap_', ref, sep=''))) +
        geom_point(alpha=.15) +
        geom_abline(slope=1,intercept=0, linetype='dashed') +
        xlab(paste('Identity with ', refs[1], ' reference', sep='')) + 
        ylab(paste('Identity with ', ref, ' reference', sep='')) + 
        scale_colour_manual(values = c(my_color_palette[['introgressed']], my_color_palette[['nonintrogressed']])) +
        theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
              axis.line=element_line(),
              #legend.position = "none",
              axis.title.x = element_text(size=18), 
              axis.title.y = element_text(size=18), 
              axis.text.x = element_text(colour="black"), 
              axis.text.y = element_text(colour="black"))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', ref, '_nongap_id_vs_', refs[1], '_nongap_id_',tag,'.png', sep=''),
           width = 8, height = 8, dpi=300)

    # labeled
    
    ggplot(regions[[ref]], aes_string(x=paste('frac_nongap_', refs[1], sep=''), y=paste('frac_nongap_', ref, sep=''))) +
        geom_text(aes(label=as.character(region_id)),hjust=0,vjust=0, cex=1) +
        geom_point(alpha=.15) +
        geom_abline(slope=1,intercept=0, linetype='dashed') +
        xlab(paste('Identity with ', refs[1], ' reference', sep='')) + 
        ylab(paste('Identity with ', ref, ' reference', sep='')) + 
        scale_colour_manual(values = c(my_color_palette[['introgressed']], my_color_palette[['nonintrogressed']])) +
        theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
              axis.line=element_line(),
              #legend.position = "none",
              axis.title.x = element_text(size=18), 
              axis.title.y = element_text(size=18), 
              axis.text.x = element_text(colour="black"), 
              axis.text.y = element_text(colour="black"))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', ref, '_nongap_id_vs_', refs[1], '_nongap_id_labeled_',tag,'.png', sep=''),
           width = 8, height = 8, dpi=300)

    ######
    ## first and second sets of filtering
    ######

    regions[[ref]]$passes_filter1 = regions[[ref]]$region_id %in%
        regions_f1[[ref]]$region_id
    regions[[ref]]$passes_filter2 = regions[[ref]]$region_id %in%
        regions_f2[[ref]]$region_id
    regions[[ref]]$passes_filters = paste(regions[[ref]]$passes_filter1,  regions[[ref]]$passes_filter2)

    # unlabeled
    
    ggplot(regions[[ref]], aes_string(x=paste('frac_nongap_', refs[1], sep=''), y=paste('frac_nongap_', ref, sep=''), colour='passes_filters')) +
        geom_point(alpha=.15) +
        geom_abline(slope=1,intercept=0, linetype='dashed') +
        xlab(paste('Identity with ', refs[1], ' reference', sep='')) + 
        ylab(paste('Identity with ', ref, ' reference', sep='')) + 
        #scale_colour_manual(values = c(my_color_palette[['introgressed']], my_color_palette[['nonintrogressed']])) +
        theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
              axis.line=element_line(),
              #legend.position = "none",
              axis.title.x = element_text(size=18), 
              axis.title.y = element_text(size=18), 
              axis.text.x = element_text(colour="black"), 
              axis.text.y = element_text(colour="black"))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', ref, '_nongap_id_vs_', refs[1], '_nongap_id_color_filtered_',tag,'.png', sep=''),
           width = 8, height = 8, dpi=300)

    # labeled
    
    ggplot(regions[[ref]], aes_string(x=paste('frac_nongap_', refs[1], sep=''), y=paste('frac_nongap_', ref, sep=''), colour='passes_filters')) +
        geom_text(aes(label=as.character(region_id)),hjust=0,vjust=0, cex=1) +
        geom_point(alpha=.15) +
        geom_abline(slope=1,intercept=0, linetype='dashed') +
        xlab(paste('Identity with ', refs[1], ' reference', sep='')) + 
        ylab(paste('Identity with ', ref, ' reference', sep='')) + 
        #scale_colour_manual(values = c(my_color_palette[['introgressed']], my_color_palette[['nonintrogressed']])) +
        theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
              axis.line=element_line(),
              #legend.position = "none",
              axis.title.x = element_text(size=18), 
              axis.title.y = element_text(size=18), 
              axis.text.x = element_text(colour="black"),
              axis.text.y = element_text(colour="black"))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', ref, '_nongap_id_vs_', refs[1], '_nongap_id_color_filtered_labeled_',tag,'.png', sep=''),
           width = 8, height = 8, dpi=300)

    ######
    ## predicted par id vs max other par id
    ######

    regions[[ref]]$max_frac_nongap_other = as.numeric(apply(regions[[ref]], 1, function(x) max(x[names(x) %like% 'frac_nongap_' & ! names(x) %like% refs[1] & ! names(x) %like% ref])))

    regions[[ref]]$max_frac_nongap_other_name = NULL
    for (row in 1:nrow(regions[[ref]])) {
        best_matches = c()
        for (ref2 in refs[2:5]) {
            if ((ref2 != ref) & isTRUE(all.equal(regions[[ref]]$max_frac_nongap_other[row], regions[[ref]][[paste('frac_nongap_', ref2, sep='')]][row], .0000001))) {
                best_matches = c(best_matches, ref2)
            }

        }

        if (length(best_matches) == 1) {
            regions[[ref]]$max_frac_nongap_other_name[row] = best_matches[1]
        }
        else {            
            regions[[ref]]$max_frac_nongap_other_name[row] = "multiple"
        }
        if (length(best_matches) == 0) {
            regions[[ref]]$max_frac_nongap_other_name[row] = "none"
        }

    }
    
    ## unlabeled and colored by which other paradoxus is best match

    ggplot(regions[[ref]], aes_string(x=paste('max_frac_nongap_', 'other', sep=''), y=paste('frac_nongap_', ref, sep=''), colour='max_frac_nongap_other_name')) +
        geom_point(alpha=.15) +
        geom_abline(slope=1,intercept=0, linetype='dashed') +
        xlab(paste('Max dentity with other paradoxus reference', sep='')) + 
        ylab(paste('Identity with ', ref, ' reference', sep='')) + 
        #scale_colour_manual(values = c(my_color_palette[['introgressed']], my_color_palette[['nonintrogressed']])) +
        theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
              axis.line=element_line(),
              #legend.position = "none",
              axis.title.x = element_text(size=18), 
              axis.title.y = element_text(size=18), 
              axis.text.x = element_text(colour="black"),
              axis.text.y = element_text(colour="black"))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', ref, '_nongap_id_vs_max_other_nongap_id_color_',tag,'.png', sep=''),
           width = 8, height = 8, dpi=300)
    
    ## labeled and colored by filtering
    
    
    ggplot(regions[[ref]], aes_string(x=paste('max_frac_nongap_', 'other', sep=''), y=paste('frac_nongap_', ref, sep=''), colour='passes_filters')) +
        geom_text(aes(label=as.character(region_id)),hjust=0,vjust=0, cex=1) +
        geom_point(alpha=.15) +
        geom_abline(slope=1,intercept=0, linetype='dashed') +
        xlab(paste('Max dentity with other paradoxus reference', sep='')) + 
        ylab(paste('Identity with ', ref, ' reference', sep='')) + 
        #scale_colour_manual(values = c(my_color_palette[['introgressed']], my_color_palette[['nonintrogressed']])) +
        theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
              axis.line=element_line(),
              #legend.position = "none",
              axis.title.x = element_text(size=18), 
              axis.title.y = element_text(size=18), 
              axis.text.x = element_text(colour="black"),
              axis.text.y = element_text(colour="black"))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', ref, '_nongap_id_vs_max_other_nongap_id_color_filtered_labeled_',tag,'.png', sep=''),
           width = 8, height = 8, dpi=300)
    
    
}



#####################
asdf

#print(regions[["CBS432"]][which(regions[["CBS432"]]$match_hmm_CBS432/regions[["CBS432"]]

for (ref in refs) {

    #####
    # color by whether region passes first round of filtering: at
    # least one HMM site that matches predicted reference but not
    # S288c; shape by whether passes second round of filtering: not on
    # a strain x chromosome that contains any of above regions (idea
    # being that if a prediction is that far off, probably the
    # training just didn't work well at all)
    #####

    regions[[ref]]$strainchrm = paste(regions[[ref]]$strain, regions[[ref]]$chromosome, sep="")
    regions[[ref]]$passes_filter1 = regions[[ref]]$count_P > 0
    bad_strainchrms = regions[[ref]][which(regions[[ref]]$passes_filter1 == FALSE),]$strainchrm
    regions[[ref]]$passes_filter2 = !(regions[[ref]]$strainchrm %in% bad_strainchrms)

    ## plot identity in terms of HMM sites
    regions[[ref]]$frac_hmm_cer = regions[[ref]][[paste("match_hmm_", refs[1], sep="")]] / regions[[ref]][["num_sites_hmm"]]
    regions[[ref]]$frac_hmm = regions[[ref]][[paste("match_hmm_", ref, sep="")]] / regions[[ref]][["num_sites_hmm"]]

    ggplot(regions[[ref]], aes(x=frac_hmm_cer, y=frac_hmm, colour=passes_filter1, shape=passes_filter2)) +
        geom_text(aes(label=as.character(region_id)),hjust=0,vjust=0, cex=1) +
        geom_point(alpha=.15) +
        geom_abline(slope=1,intercept=0, linetype='dashed') +
        xlab(paste('Identity with ', refs[1], ' reference', sep='')) + 
        ylab(paste('Identity with ', ref, ' reference', sep='')) + 
        scale_colour_manual(values = c(my_color_palette[['introgressed']], my_color_palette[['nonintrogressed']])) +
        theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
              axis.line=element_line(),
              #legend.position = "none",
              axis.title.x = element_text(size=18), 
              axis.title.y = element_text(size=18), 
              axis.text.x = element_text(colour="black"), 
              axis.text.y = element_text(colour="black"))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', ref, '_hmm_id_vs_', refs[1], '_hmm_id_',tag,'.png', sep=''),
           width = 8, height = 8, dpi=300)

    ## plot identity in terms of HMM sites, only those that pass filters
    ggplot(regions[[ref]][which(regions[[ref]]$passes_filter2),], aes(x=frac_hmm_cer, y=frac_hmm, size=num_sites_hmm)) +
        geom_text(aes(label=as.character(region_id)),hjust=0,vjust=0, cex=1) +
        geom_point(alpha=.15) +
        scale_size_continuous(range = c(1, 5)) +
        geom_abline(slope=1,intercept=0, linetype='dashed') +
        xlab(paste('Identity with ', refs[1], ' reference', sep='')) + 
        ylab(paste('Identity with ', ref, ' reference', sep='')) + 
        scale_colour_manual(values = c(my_color_palette[['introgressed']], my_color_palette[['nonintrogressed']])) +
        theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
              axis.line=element_line(),
              #legend.position = "none",
              axis.title.x = element_text(size=18), 
              axis.title.y = element_text(size=18), 
              axis.text.x = element_text(colour="black"), 
              axis.text.y = element_text(colour="black"))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', ref, '_hmm_id_vs_', refs[1], '_hmm_id_',tag,'_filtered.png', sep=''),
           width = 8, height = 8, dpi=300)
    

    ## plot identity in terms of all nongap sites
    regions[[ref]]$frac_nongap_cer = regions[[ref]][[paste("match_nongap_", refs[1], sep="")]] / regions[[ref]][[paste("num_sites_nongap_", refs[1], sep='')]]
    regions[[ref]]$frac_nongap = regions[[ref]][[paste("match_nongap_", ref, sep="")]] / regions[[ref]][[paste("num_sites_nongap_", ref, sep='')]]

    ggplot(regions[[ref]], aes(x=frac_nongap_cer, y=frac_nongap, colour=passes_filter1, shape=passes_filter2)) +
        geom_text(aes(label=as.character(region_id)),hjust=0,vjust=0, cex=1) +
        geom_point(alpha=.15) +
        geom_abline(slope=1,intercept=0, linetype='dashed') +
        xlab(paste('Identity with ', refs[1], ' reference', sep='')) + 
        ylab(paste('Identity with ', ref, ' reference', sep='')) + 
        scale_colour_manual(values = c(my_color_palette[['introgressed']], my_color_palette[['nonintrogressed']])) +
        theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
              axis.line=element_line(),
              #legend.position = "none",
              axis.title.x = element_text(size=18), 
              axis.title.y = element_text(size=18), 
              axis.text.x = element_text(colour="black"), 
              axis.text.y = element_text(colour="black"))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', ref, '_nongap_id_vs_', refs[1], '_nongap_id_',tag,'.png', sep=''),
           width = 8, height = 8, dpi=300)

    ## plot identity in terms of only nongap sites, only those that pass filters
    ggplot(regions[[ref]][which(regions[[ref]]$passes_filter2),], aes(x=frac_nongap_cer, y=frac_nongap, size=end-start+1)) +
        geom_text(aes(label=as.character(region_id)),hjust=0,vjust=0, cex=1) +
        geom_point(alpha=.15) +
        scale_size_continuous(range = c(1, 5)) +
        geom_abline(slope=1,intercept=0, linetype='dashed') +
        xlab(paste('Identity with ', refs[1], ' reference', sep='')) + 
        ylab(paste('Identity with ', ref, ' reference', sep='')) + 
        scale_colour_manual(values = c(my_color_palette[['introgressed']], my_color_palette[['nonintrogressed']])) +
        theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
              axis.line=element_line(),
              #legend.position = "none",
              axis.title.x = element_text(size=18), 
              axis.title.y = element_text(size=18), 
              axis.text.x = element_text(colour="black"), 
              axis.text.y = element_text(colour="black"))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', ref, '_nongap_id_vs_', refs[1], '_nongap_id_',tag,'_filtered.png', sep=''),
           width = 8, height = 8, dpi=300)
}


##=====
# predicted species id vs cer id
##=====

for (ref in refs) {
    regions[[ref]]$frac_hmm_cer = regions[[ref]][[paste("match_hmm_", refs[1], sep="")]] / regions[[ref]][["num_sites_hmm"]]
    regions[[ref]]$frac_hmm = regions[[ref]][[paste("match_hmm_", ref, sep="")]] / regions[[ref]][["num_sites_hmm"]]
    ggplot(regions[[ref]], aes(x=frac_hmm_cer, y=frac_hmm, colour='x')) +
        geom_point(alpha=.15) +
        geom_abline(slope=1,intercept=0, linetype='dashed') +
        xlab(paste('Identity with ', refs[1], ' reference', sep='')) + 
        ylab(paste('Identity with ', ref, ' reference', sep='')) + 
        scale_colour_manual(values = c(my_color_palette[['introgressed']])) +
        theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
              axis.line=element_line(),
              legend.position = "none",
              axis.title.x = element_text(size=18), 
              axis.title.y = element_text(size=18), 
              axis.text.x = element_text(colour="black"), 
              axis.text.y = element_text(colour="black"))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', ref, '_hmm_id_vs_', refs[1], '_hmm_id_',tag,'.png', sep=''),
           width = 8, height = 8)
}










#################################






asfd

ggplot(regions_S288c, (aes(x=match_hmm_S288c/num_sites_hmm, y=match_hmm_N_45/num_sites_hmm, label=region_id, colour='x', size=num_sites_hmm))) + geom_point(alpha=.15) +
 geom_text(aes(label=as.character(region_id)),hjust=0,vjust=0, cex=.2) +
    scale_size_continuous(range = c(1, 5)) +
    #coord_cartesian(xlim=c(.6, 1),ylim=c(.6,1))
    geom_abline(slope=1,intercept=0, linetype='dashed') +
    xlab('Identity with S288c reference') + 
    ylab('Identity with CBS432 reference') + 
    scale_colour_manual(values = c(my_color_palette[['introgressed']])) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          legend.position = "none",
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"))
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/S-N_45_hmm_id_vs_S288c_hmm_id_',tag,'.pdf',sep=''), width = 8, height = 8)

ggplot(regions_N, (aes(x=match_hmm_S288c/num_sites_hmm, y=match_hmm_N_45/num_sites_hmm, label=region_id, colour='x', size=num_sites_hmm))) + geom_point(alpha=.15) +
 geom_text(aes(label=as.character(region_id)),hjust=0,vjust=0, cex=.2) +
    scale_size_continuous(range = c(1, 5)) +
    #coord_cartesian(xlim=c(.6, 1),ylim=c(.6,1))
    geom_abline(slope=1,intercept=0, linetype='dashed') +
    xlab('Identity with S288c reference') + 
    ylab('Identity with CBS432 reference') + 
    scale_colour_manual(values = c(my_color_palette[['introgressed']])) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          legend.position = "none",
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"))
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/N_45_hmm_id_vs_S288c_hmm_id_',tag,'.pdf',sep=''), width = 8, height = 8)


ggplot(regions_C, (aes(x=match_hmm_S288c/num_sites_hmm, y=match_hmm_CBS432/num_sites_hmm, label=region_id, colour='x', size=num_sites_hmm))) + geom_point(alpha=.15) +
 geom_text(aes(label=as.character(region_id)),hjust=0,vjust=0, cex=.2) +
    scale_size_continuous(range = c(1, 5)) +
    #coord_cartesian(xlim=c(.6, 1),ylim=c(.6,1))
    geom_abline(slope=1,intercept=0, linetype='dashed') +
    xlab('Identity with S288c reference') + 
    ylab('Identity with CBS432 reference') + 
    scale_colour_manual(values = c(my_color_palette[['introgressed']])) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          legend.position = "none",
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"))
ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/CBS432_hmm_id_vs_S288c_hmm_id_',tag,'.pdf',sep=''), width = 8, height = 8)
