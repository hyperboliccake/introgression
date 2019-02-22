library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(grid)
source('../my_color_palette.R')

a = read.table('ils_rosenberg_out.txt', sep=' ', header=F, stringsAsFactors=FALSE)
names(a) = c('T', 'prob')

a = a[which(a$prob > 0),]

ggplot(a, aes(x=log10(T), y=log10(prob), colour='a')) + 
    geom_point(size=3, alpha=.7) +
    geom_line(size=1.5) +
    xlab('log10(T)') +
    ylab('log10(Probability of ILS)') +
    scale_colour_manual(values = c("#9E1042")) +
    #scale_x_continuous(expand=c(0,0)) +
    geom_vline(xintercept = log10(46.875), linetype = 'dashed', colour='black',
               size = 1) +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(colour="gray90"),
          panel.grid.major=element_line(colour="gray80"),
          axis.ticks=element_line(colour="black"),
          axis.line=element_line(),
          legend.position = "none",
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=18),
          axis.text.x = element_text(size=18, colour="black"),
          axis.text.y = element_text(size=18, colour="black"),
          plot.margin=unit(c(.5,.5,.5,.5),"in"))
ggsave('ils_rosenberg.png', width = 8, height = 5)
