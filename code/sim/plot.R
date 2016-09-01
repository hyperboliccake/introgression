library(ggplot2)
library(reshape2)
library(RColorBrewer)


make_plots = function(d)
{
    # 

}

ids = 1:18
results_dir = '../../results/sim/run3/'
prefix = 'sim_out_'
suffix = '_summary.txt'
args = read.table('sim_multi_model_args.txt', sep = ' ', header = F)
for (id in ids)
{
    a = read.table(results_dir + prefix + id + suffix, sep = '\t', header = T)
    
}
     
