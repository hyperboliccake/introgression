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

results_dir = '../../results/sim/analyze/'
prefix = 'sim_out_'
suffix = '_summary.txt'
a = read.table(paste(results_dir, prefix, 'all', suffix, sep=''), sep = '\t', header = T)
n = names(a)
i = which(n == "row_type") 
b = melt(a, id = 1:i)
d = reshape(b, timevar='row_type', idvar = c(names(b)[1:i-1], names(b)[i+1]), direction = 'wide')
#d$value.lower_se = d$value.mean - d$value.std_err
#d$value.upper_se = d$value.mean + d$value.std_err

for (v in levels(d$variable))
{
    print(v)
    barplot(d[d$variable == v,], 'tag', 'value.mean', 'factor(mig_par)', 'value.lower', 'value.upper', paste('../../results/sim/plots/',v,'.pdf',sep=''))
}
