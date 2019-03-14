library(ggplot2)
library(viridis)
require(grDevices)
library(grid)
source('../my_color_palette.R')

#####
# plotting multiple sets of results for same method
#####

args = commandArgs(trailingOnly=TRUE)

id = args[1]

pred_ids = args[2:length(args)]

pred_args = read.table('predict_compare_args.txt', header=F, sep=' ')
names(pred_args) = c('sim_id', 'pred_id', 'training_threshold', 'posterior_threshold', 'expected_length', 'expected_tracts', 'ref')

sim_args = read.table('sim_compare_args.txt', header=F, sep=' ')
names(sim_args) = c('sim_id', 'tree', 'species_to', 'n_species_to', 'n0_species_to', 'species_from', 'n_species_from', 'n0_species_from', 'migration_rate', 'nsites', 'recomb_rate', 'outcross_rate', 'nreps', 'name_species_to', 'name_species_from')

a = data.frame()

for (pred_id in pred_ids) {
    
    sim_id = pred_args[which(pred_args$pred_id == pred_id),]$sim_id

    print(sim_id)
    fn = paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/sim/sim_out_', sim_id, '_roc_predicted_', 
               pred_id, '.txt', sep='')
    print(fn)
    b = read.table(fn, sep='\t', header=T)
    b$migration_rate = sim_args[which(sim_args$sim_id == sim_id),]$migration_rate
    ## b$method = 'posterior'

    a = rbind(a, b)

    #fn = paste('../../results/sim/sim_out_', sim_id, '_introgressed_predicted_viterbi_', 
    #           pred_id, '.txt', sep='')
    #print(fn)
    #b = read.table(fn, sep='\t', header=T)
    #b$migration_rate = sim_args[which(sim_args$sim_id == sim_id),]$migration_rate
    #b$method = 'posterior'

    #a = rbind(a, b)
}


ggplot(a, aes(x=fpr,y=tpr,colour=as.factor(migration_rate))) + geom_line() + geom_point(alpha=.5) +
          xlab('FPR') + ylab('TPR') +
          scale_colour_viridis(discrete=TRUE) +
	  geom_abline(intercept = 0, slope = 1, linetype='dashed') + 
	  labs(colour="migration rate") +
	  coord_cartesian(xlim=c(0,1),ylim=c(0,1)) + 
          theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.line=element_line(),
          axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"))

ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/sim/roc_plots/roc_', id, '.pdf', sep=''), height=9, width=9)

ggplot(a, aes(x=fpr,y=tpr,colour=as.factor(migration_rate))) + geom_line() + geom_point() +
          xlab('FPR') + ylab('TPR') +
          scale_colour_viridis(discrete=TRUE) +
	  geom_abline(intercept = 0, slope = 1, linetype='dashed') + 
	  labs(colour="migration rate") +
	  coord_cartesian(xlim=c(0,.05),ylim=c(0,.05)) + 
          theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.line=element_line(),
          axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))

ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/sim/roc_plots/roc_', id, '_zoomed.pdf', sep=''), height=9, width=9)



#####
# plotting results from different methods
#####

args = commandArgs(trailingOnly=TRUE)

sim_id = args[1]
pred_id = args[2]

fn1 = paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/sim/sim_out_', sim_id, '_roc_predicted_', 
            pred_id, '.txt', sep='')
a1 = read.table(fn1, sep='\t', header=T)
a1$method = 'HMM'
a1 = a1[order(a1$tpr),]

fn2 = paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/sim/sim_out_', sim_id, '_roc_predicted_phylohmm_', 
            pred_id, '.txt', sep='')
a2 = read.table(fn2, sep='\t', header=T)
a2$method = 'PhyloNet-HMM'
a2 = a2[order(a2$tpr),]




a = rbind(a1, a2)

ggplot(a, aes(x=fpr,y=tpr,shape=method,colour=method)) +
    geom_line(size=2) +
    geom_point(size=4, alpha=.7) +
    xlab('False positive rate') + ylab('True positive rate') +
    scale_colour_manual(values = c(my_color_palette[['misc1']],
                                   my_color_palette[['misc2']])) +
    scale_shape_manual(values=c(19,17)) +
    geom_abline(intercept = 0, slope = 1, linetype='dashed') + 
    labs(colour="Method", shape="Method") +
          scale_x_continuous(expand=c(0,.01)) + 
	  scale_y_continuous(expand=c(0,.01)) + 
          theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          axis.line=element_line(),
          axis.text.x = element_text(colour="black", size=16),
          axis.text.y = element_text(colour="black", size=16),
          axis.title.x = element_text(colour="black", size=20),
          axis.title.y = element_text(colour="black", size=20),
          plot.margin=unit(c(.5,.5,.5,.5),"in"),
          legend.position = c(.8,.15),
          legend.title = element_text(size=18),
          legend.text = element_text(size=16), 
          legend.key = element_rect(fill = "transparent"))

ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/sim/roc_plots/roc_', sim_id, '_', pred_id, '.png', sep=''), height=6, width=6)

