library(ggplot2)
library(reshape2)
library(RColorBrewer)

barplot = function(d, xvar, yvar, fillvar, lower, upper, fn)
{
    ggplot(d, aes_string(x=xvar, y=yvar, fill=fillvar)) + geom_bar(stat="identity", position="dodge") +
        geom_errorbar(aes_string(ymin = lower, ymax = upper), position='dodge') +
        theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))
    
    ggsave(fn, width = 12, height = 7)
    
}

results_dir = '../../results/sim/run_3/'
prefix = 'sim_out_'
suffix = '_summary.txt'
a = read.table(paste(results_dir, prefix, 'all', suffix, sep=''), sep = '\t', header = T)
n = names(a)
b = melt(a, id = 1:12)
d = reshape(b, timevar='row_type', idvar = c(names(b)[1:11], names(b)[13]), direction = 'wide')
#d$value.lower_se = d$value.mean - d$value.std_err
#d$value.upper_se = d$value.mean + d$value.std_err

for (v in levels(d$variable))
{
    barplot(d[d$variable == v,], 'tag', 'value.mean', 'factor(par_cer_migration)', 'value.lower', 'value.upper', paste('../../results/sim/run_3/plots/',v,'.pdf',sep=''))
}
