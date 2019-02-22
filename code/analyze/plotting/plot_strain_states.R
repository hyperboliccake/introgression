library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(hexbin)
library(data.table)
library(dplyr)
library(VennDiagram)
source('../my_color_palette.R')

args = commandArgs(trailingOnly=TRUE)
tag = args[1]
refs = c("CBS432", "N_45", "DBVPG6304", "UWOPS91_917_1")

## generate dataframe a which has a column for all alternative state
## names (individual paradoxus reference as well as all possible
## combinations), a column for the strain, and a column (total) for
## the total number of bases in each category
a = data.frame()
for (ri in 1:length(refs)) {
    ax = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/',
                          'results/analysis/', tag, '/blocks_', refs[ri], '_', tag,
                          '_filtered2intermediate.txt',
                          sep=''), sep='\t', header=T, stringsAsFactors=F)

    acur = ax %>%
        group_by(.dots=c('alternative_states', 'strain')) %>%
        summarize(total = sum(end-start+1))

    a = bind_rows(a, acur) %>%
        group_by(.dots=c('alternative_states', 'strain')) %>%
        summarize(total = sum(total))
}

a = a %>%
    group_by(strain) %>%
    mutate(frac = total/sum(total))

grand_total = a %>%
    group_by(strain) %>%
    summarize(grand_total = sum(total))

a = merge(a, grand_total)

## now we can make something like a structure plot, but showing the
## amount of introgression from each paradoxus reference in each
## strain

## east = refs 2&3 red-orange-gold
## west = refs 4&5 blue-teal-green

## these are all the individuals reference plus the two combinations
## we care about
sl = c(refs[1], paste(refs[1], refs[2], sep=','), refs[2],
       refs[3], paste(refs[3], refs[4], sep=','), refs[4])
ax = unique(a$alternative_states)
## these are all the remaining combinations which we don't care about
ay = setdiff(ax, sl)
## put them together, with the ones we care about ordered up front
sl = c(sl, ay)
a$alternative_states = factor(a$alternative_states, levels = rev(sl))
## color the categories we care about appropriately, and make all the
## other ones the same gray
cols = c(c("#E13939", "#DA7921", "#E1A939", "#007CEB", "#00A89C", "#009E2A"),
         rep('gray75', length(sl) - 6))
cols = rev(cols)

a$index1 = 1
a$index2 = 1
for (row in 1:nrow(a)) {
    a[row,]$index1 = sum(a[which(a$strain == a[row,]$strain &
                                (a$alternative_states == refs[1] |
                                 a$alternative_states == refs[2] |
                                 a$alternative_states ==
                                 paste(refs[1], refs[2], sep=','))),]$frac)
    a[row,]$index2 = a[row,]$index1 /
        sum(a[which(a$strain == a[row,]$strain &
                                (a$alternative_states == refs[3] |
                                 a$alternative_states == refs[4] |
                                 a$alternative_states ==
                                 paste(refs[3], refs[4], sep=','))),]$frac)
}

## plot totals
ggplot(a, aes(x = reorder(strain, -grand_total), y = total / 12071326,
              fill = alternative_states)) +
    geom_bar(stat='identity') +
    xlab('strain') +
    ylab('fraction of genome introgressed') + 
    scale_fill_manual(values = cols) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"),
          panel.grid.major=element_line(colour="gray92"),
          legend.position="none",
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(size=7, angle=90, vjust=1, hjust=1,
                                     colour="black"), 
          axis.text.y = element_text(colour="black"))

ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',
             tag, '/plots/strain_states.png',
             sep=''),
       width = 16, height = 8, dpi=300)

## plots fraction
ggplot(a, aes(x = reorder(strain, index1), y = frac,
              fill = alternative_states)) +
    geom_bar(stat='identity') +
    xlab('strain') +
    ylab('fraction of introgression') + 
    #geom_label(aes(population, frac_total + shift,
    #               label = base_total, fill=NULL), data = totals) +
    scale_fill_manual(values = cols) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
          legend.position="none",
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(size=7, angle=90, vjust=1, hjust=1,
                                     colour="black"), 
          axis.text.y = element_text(colour="black"))

ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',
             tag, '/plots/strain_states_frac_order1.png',
             sep=''),
       width = 16, height = 8, dpi=300)

## plots fraction
ggplot(a, aes(x = reorder(strain, index2), y = frac,
              fill = alternative_states)) +
    geom_bar(stat='identity') +
    xlab('strain') +
    ylab('fraction of introgression') + 
    #geom_label(aes(population, frac_total + shift,
    #               label = base_total, fill=NULL), data = totals) +
    scale_fill_manual(values = cols) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
          legend.position="none",
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(size=7, angle=90, vjust=1, hjust=1,
                                     colour="black"), 
          axis.text.y = element_text(colour="black"))

ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',
             tag, '/plots/strain_states_frac_order2.png',
             sep=''),
       width = 16, height = 8, dpi=300)


asdf



a = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/', tag, '/state_counts_by_strain.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)

##=====
# something sort of like structure plot - bar chart showing fraction
# of introgression assigned to each reference
##=====

a$id = paste(a$population, a$strain)

## for each paradoxus reference, get the fraction of all introgressed
## bases that we can assign specifically to that reference, at
## different stages of filtering
for (ri in 1:length(refs)) {
    a[[paste("frac_bases_", refs[ri], sep='')]] = a[[paste("num_bases_", refs[ri], sep='')]] / a[["num_bases_total"]]
    a[[paste("frac_bases_", refs[ri], '_filtered1', sep='')]] = a[[paste("num_bases_", refs[ri], '_filtered1', sep='')]] / a[["num_bases_total_filtered1"]]
    a[[paste("frac_bases_", refs[ri], '_filtered2i', sep='')]] = a[[paste("num_bases_", refs[ri], '_filtered2', sep='')]] / a[["num_bases_total_filtered1"]]
    a[[paste("frac_bases_", refs[ri], '_filtered2', sep='')]] = a[[paste("num_bases_", refs[ri], '_filtered2', sep='')]] / a[["num_bases_total_filtered2"]]

    a[[paste("frac_bases_3_filtered2i", sep='')]] = a[[paste("num_bases_3_filtered2i", sep='')]] / a[["num_bases_total_filtered1"]]
    a[[paste("frac_bases_4_filtered2i", sep='')]] = a[[paste("num_bases_4_filtered2i", sep='')]] / a[["num_bases_total_filtered1"]]
    for (rj in (ri+1):length(refs)) {
        ## fuck this
        if (ri < rj && rj <= length(refs)) {
            a[[paste("frac_bases_", refs[ri], '_or_', refs[rj], '_filtered2i', sep='')]] = a[[paste("num_bases_", refs[ri], '_or_', refs[rj], '_filtered2i', sep='')]] / a[["num_bases_total_filtered1"]]
    }
    }
}

#a[["frac_bases_west"]] = a[["frac_bases_DBVPG6304"]] + a[["frac_bases_UWOPS91_917_1"]]
#a[["frac_bases_west_filtered1"]] = a[["frac_bases_DBVPG6304_filtered1"]] + a[["frac_bases_UWOPS91_917_1_filtered1"]]
#a[["frac_bases_west_filtered2"]] = a[["frac_bases_DBVPG6304_filtered2"]] + a[["frac_bases_UWOPS91_917_1_filtered2"]]


##=====
# scatter plot of euro+east introgressed bases vs american+hawaiian introgressed bases
##=====

a[["num_bases_west"]] =
    a[["num_bases_DBVPG6304_filtered2"]] +
    a[["num_bases_DBVPG6304_or_UWOPS91_917_1_filtered2i"]] +
    a[["num_bases_UWOPS91_917_1_filtered2"]]
a[["num_bases_east"]] =
    a[["num_bases_CBS432_filtered2"]] +
    a[["num_bases_CBS432_or_N_45_filtered2i"]] +
    a[["num_bases_N_45_filtered2"]]
#b = melt(a, id.var=c("strain", "population", "geographic_origin", "environmental_origin"))

ggplot(a, aes(x = log(num_bases_west), y=log(num_bases_east), colour=population)) + geom_point(size=3, alpha=.4) +
        geom_abline(slope=1,intercept=0, linetype='dashed') +


ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', 'east_vs_west_',tag,'.png', sep=''),
       width = 12, height = 8, dpi=300)

asdf
######
## fraction assigned to each reference after filtering stage 1, only non-mosaic populations
######

b = melt(a, id.var=c("strain", "population", "geographic_origin", "environmental_origin", "id"))
b = b[which(b$variable %in% c('frac_bases_CBS432', 'frac_bases_DBVPG6304', 'frac_bases_UWOPS91_917_1', 'frac_bases_N_45')),]
b = b[which(!(b$population == "mosaic")),] 

ggplot(b, aes(x=id, y=value, fill=variable)) +
    geom_bar(stat="identity") +
    xlab('') +
    ylab('fraction of introgressed bases (before filtering)') + 
    #scale_colour_manual(values = c(my_color_palette[['introgressed']], my_color_palette[['nonintrogressed']])) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(size=8,angle = 90,vjust = 1,hjust=1,colour="black"), 
          axis.text.y = element_text(colour="black"))

ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', 'strain_states_nonmosaic_',tag,'.png', sep=''),
       width = 12, height = 8, dpi=300)

######
## fraction assigned to each reference after filtering stage 1, only non-mosaic populations
######

b = melt(a, id.var=c("strain", "population", "geographic_origin", "environmental_origin", "id"))
b = b[which(b$variable %in% c('frac_bases_CBS432_filtered1', 'frac_bases_DBVPG6304_filtered1', 'frac_bases_UWOPS91_917_1_filtered1', 'frac_bases_N_45_filtered1')),]
b = b[which(!(b$population == "mosaic")),] 

ggplot(b, aes(x=id, y=value, fill=variable)) +
    geom_bar(stat="identity") +
    xlab('') +
    ylab('fraction of introgressed bases (after filtering stage 1)') + 
    #scale_colour_manual(values = c(my_color_palette[['introgressed']], my_color_palette[['nonintrogressed']])) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(size=8,angle = 90,vjust = 1,hjust=1,colour="black"), 
          axis.text.y = element_text(colour="black"))

ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', 'strain_states_1_nonmosaic_',tag,'.png', sep=''),
       width = 12, height = 8, dpi=300)

######
## fraction assigned to each reference after filtering stage 2, only non-mosaic populations
######

b = melt(a, id.var=c("strain", "population", "geographic_origin", "environmental_origin", "id"))
b = b[which(b$variable %in% c('frac_bases_CBS432_filtered2', 'frac_bases_DBVPG6304_filtered2', 'frac_bases_UWOPS91_917_1_filtered2', 'frac_bases_N_45_filtered2')),]
b = b[which(!(b$population == "mosaic")),] 

ggplot(b, aes(x=id, y=value, fill=variable)) +
    geom_bar(stat="identity") +
    xlab('') +
    ylab('fraction of introgressed bases (after filtering stage 2)') + 
    #scale_colour_manual(values = c(my_color_palette[['introgressed']], my_color_palette[['nonintrogressed']])) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(size=8,angle = 90,vjust = 1,hjust=1,colour="black"), 
          axis.text.y = element_text(colour="black"))

ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', 'strain_states_2_nonmosaic_',tag,'.png', sep=''),
       width = 12, height = 8, dpi=300)

######
## fraction assigned to each reference after filtering stage 1, only non-mosaic populations, but further specified by filtering stage 2
######

b = melt(a, id.var=c("strain", "population", "geographic_origin", "environmental_origin", "id"))
l = rev(c('frac_bases_CBS432_filtered2i',
      'frac_bases_N_45_filtered2i',
      'frac_bases_DBVPG6304_filtered2i',
      'frac_bases_UWOPS91_917_1_filtered2i',
      'frac_bases_CBS432_or_N_45_filtered2i',
      'frac_bases_DBVPG6304_or_UWOPS91_917_1_filtered2i',
      'frac_bases_CBS432_or_DBVPG6304_filtered2i',
      'frac_bases_CBS432_or_UWOPS91_917_1_filtered2i',
      'frac_bases_DBVPG6304_or_N_45_filtered2i',
      'frac_bases_N_45_or_UWOPS91_917_1_filtered2i',
      'frac_bases_3_filtered2i',
      'frac_bases_4_filtered2i'
      ))

b = b[which(b$variable %in% l),]
b = b[which(!(b$population == "mosaic")),] 

b$variable = factor(b$variable, levels = l)

cols = rev(c('#3E3EEE', '#499749', '#B34242', '#D6A73A', '#187994', '#E58211', 'gray50', 'gray60', 'gray70', 'gray80', 'antiquewhite2', 'white'))


a = a[order(a$id),]
xcols = rep("black", nrow(a))
xcols[which(a$environmental_origin %like% "Clinical")] = "red"
xcols = xcols[which(a$population != "mosaic")]
print(xcols)

ggplot(b, aes(x=id, y=value, fill=variable)) +
    geom_bar(stat="identity") +
    xlab('') +
    ylab('fraction of introgressed bases (after filtering stage 1)') + 
    scale_fill_manual(values = cols) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(size=8,angle = 90,vjust = 1,hjust=1,colour=xcols), 
          axis.text.y = element_text(colour="black"))

ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', 'strain_states_2i_nonmosaic_',tag,'.png', sep=''),
       width = 12, height = 8, dpi=300)


#############################
asdf

a[["frac_bases_west"]] = a[["frac_bases_DBVPG6304"]] + a[["frac_bases_UWOPS91_917_1"]]
b = melt(a, id.var=c("strain", "population", "geographic_origin", "environmental_origin"))
b = b[which(!(b$population == "mosaic")),] 
b$id = paste(b$population, b$strain)

ggplot(b, aes(x=id, y=value, fill=variable)) +
    geom_bar(stat="identity") +
    #xlab(paste('Identity with ', refs[1], ' reference', sep='')) + 
    #ylab(paste('Identity with ', ref, ' reference', sep='')) + 
    #scale_colour_manual(values = c(my_color_palette[['introgressed']], my_color_palette[['nonintrogressed']])) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(size=8,angle = 45,vjust = 1,hjust=1,colour="black"), 
          axis.text.y = element_text(colour="black"))

ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', 'strain_states_west_',tag,'.png', sep=''),
       width = 12, height = 8, dpi=300)

asdf

#b = melt(a, id.var=c("strain", "population", "geographic_origin", "environmental_origin"))
a$clinical = FALSE
#print(which("Clinical" %in% a$environmental_origin))
a[a$environmental_origin %like% "Clinical",]$clinical=TRUE
                                        #b = a[which(a$population != "mosaic"),]
ggplot(a, aes(x=reorder(id, -frac_bases_west), fill = clinical, y=frac_bases_west)) +
    geom_bar(stat="identity") +
    #xlab(paste('Identity with ', refs[1], ' reference', sep='')) + 
    #ylab(paste('Identity with ', ref, ' reference', sep='')) + 
    #scale_colour_manual(values = c(my_color_palette[['introgressed']], my_color_palette[['nonintrogressed']])) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray92"), panel.grid.major=element_line(colour="gray92"),
          axis.line=element_line(),
          axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          axis.text.x = element_text(size=8,angle = 45,vjust = 1,hjust=1,colour="black"), 
          axis.text.y = element_text(colour="black"))

ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/', 'strain_states_west_',tag,'.png', sep=''),
       width = 12, height = 8, dpi=300)

######
## fraction assigned to each reference after filtering 2 (noninclusive)
######




venn.plot = draw.quad.venn(category = refs,
                           direct.area=T,
                           area.vector=c(
                               sum(a$num_bases_CBS432_filtered2),
                               sum(a$num_bases_DBVPG6304_filtered2),
                               sum(a$num_bases_N_45_filtered2),
                               sum(a$num_bases_UWOPS91_917_1_filtered2),
                               sum(a$num_bases_CBS432_or_DBVPG6304_filtered2i),
                               sum(a$num_bases_CBS432_or_N_45_filtered2i),
                               sum(a$num_bases_CBS432_or_UWOPS91_917_1_filtered2i),
                               sum(a$num_bases_DBVPG6304_or_N_45_filtered2i),
                               sum(a$num_bases_DBVPG6304_or_UWOPS91_917_1_filtered2i),
                               sum(a$num_bases_N_45_or_UWOPS91_917_1_filtered2i),
                               sum(a$num_bases_CBS432_or_DBVPG6304_or_N_45_filtered2i),
                               sum(a$num_bases_CBS432_or_DBVPG6304_or_UWOPS91_917_1_filtered2i),
                               sum(a$num_bases_CBS432_or_N_45_or_UWOPS91_917_1_filtered2i),
                               sum(a$num_bases_DBVPG6304_or_N_45_or_UWOPS91_917_1_filtered2i),
                               sum(a$num_bases_CBS432_or_DBVPG6304_or_N_45_or_UWOPS91_917_1_filtered2i)))

