library(ggplot2)
library(viridis)
require(grDevices)
library(grid)
source('../my_color_palette.R')

#####
# plotting multiple sets of results, from same prediction method for
# different simulation parameters (specifically migration rate)
#####

pred_args = read.table('predict_compare_args.txt', header=F, sep=' ')
names(pred_args) = c('sim_id', 'pred_id', 'training_threshold', 'posterior_threshold', 'expected_length', 'expected_tracts', 'ref')

sim_args = read.table('sim_compare_args.txt', header=F, sep=' ')
names(sim_args) = c('sim_id', 'tree', 'species_to', 'n_species_to', 'n0_species_to', 'species_from', 'n_species_from', 'n0_species_from', 'migration_rate', 'nsites', 'recomb_rate', 'outcross_rate', 'nreps', 'name_species_to', 'name_species_from')

a = data.frame()

for (s in sim_args$sim_id) {

    # assume only one
    pred_id = pred_args[which(pred_args$sim_id == s),]$pred_id
    
    fn = paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/sim/sim_out_', s, '_roc_predicted_', 
               pred_id, '.txt', sep='')
    print(fn)
    b = read.table(fn, sep='\t', header=T, stringsAsFactors=FALSE)
    b$migration_rate = sim_args[which(sim_args$sim_id == s),]$migration_rate[1]

    if (b$migration_rate < .00000001 & b$migration_rate > 0) {
        a = rbind(a, b)
    }
}

ggplot(a, aes(x=fpr,y=tpr,colour=as.factor(migration_rate))) +
    geom_line(size=2) + geom_point(size=4, alpha=.7) +
          xlab('False positive rate') + ylab('True positive rate') +
          scale_colour_viridis(discrete=TRUE) +
	  geom_abline(intercept = 0, slope = 1, linetype='dashed') + 
	  labs(colour="Migration rate") +
	  coord_cartesian(xlim=c(0,1),ylim=c(0,1)) + 
          theme(panel.background=element_rect(fill="white"),
          legend.key = element_rect(fill = "white"),
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.line=element_line(),
          legend.position = c(.8,.2),
          legend.title=element_text(size=18),
          legend.text=element_text(size=16),
          plot.margin=unit(c(.5,.5,.5,.5),"in"),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          axis.text.x = element_text(colour="black", size=16),
          axis.text.y = element_text(colour="black", size=16))

ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/sim/roc_plots/roc_compare_sims.png', sep=''), height=6, width=6)
