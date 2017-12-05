library(ggplot2)
library(viridis)
require(grDevices)

args = commandArgs(trailingOnly=TRUE)

sim_id = args[1]
pred_id = args[2]

fn1 = paste('../../results/sim/sim_out_', sim_id, '_roc_predicted_', 
            pred_id, '.01.txt', sep='')
a1 = read.table(fn1, sep='\t', header=T)
a1$method = 'HMM'

fn2 = paste('../../results/sim/sim_out_', sim_id, '_roc_predicted_phylohmm_', 
            pred_id, '.txt', sep='')
a2 = read.table(fn2, sep='\t', header=T)
a2$method = 'PhyloNet-HMM'

a = rbind(a1, a2)

print(dim(a))
print(a)

ggplot(a, aes(x=fpr,y=tpr,colour=method)) + geom_line()

ggsave(paste('../../results/sim/roc_plots/roc_', sim_id, '_', pred_id, '.pdf', sep=''), height=9, width=9)

