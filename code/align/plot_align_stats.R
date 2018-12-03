library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)

a = read.table('/tigress/anneec/projects/introgression/alignments/par4/mafft_stats_summary.txt', header=T, stringsAsFactors=F)

## plot number of alignment columns with fewer than all sequences
a$not_all_aligned = a$num_align_columns_1 + a$num_align_columns_2 +
    a$num_align_columns_3 + a$num_align_columns_4 +
    a$num_align_columns_5
ggplot(a, aes(x=chromosome, y=not_all_aligned,  label=strain)) + geom_text(hjust=0,vjust=0,cex=2)
ggsave('/tigress/anneec/projects/introgression/alignments/par4/mafft_stats_summary.pdf', width = 12, height = 7)

ggplot(a, aes(x=chromosome, y=frac_S288c_DBVPG6304, label=strain)) + geom_text(hjust=0,vjust=0,cex=2)
ggsave('/tigress/anneec/projects/introgression/alignments/par4/mafft_stats_summary.pdf', width = 12, height = 7)

