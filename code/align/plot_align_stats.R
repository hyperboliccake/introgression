library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)

a = read.table('/tigress/anneec/projects/introgression/alignments/genbank/mafft_stats_summary.txt', header=T, stringsAsFactors=F)
a$not_all_aligned = a$num_align_columns_1 + a$num_align_columns_2

# plot number of alignment columns with fewer than 3 sequences
ggplot(a, aes(x=chromosome, y=not_all_aligned,  label=strain)) + geom_text(hjust=0,vjust=0,cex=2)
ggsave('/tigress/anneec/projects/introgression/alignments/genbank/mafft_stats_summary.pdf', width = 12, height = 7)
