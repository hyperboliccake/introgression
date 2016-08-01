# get alignments between each strain and the cerevisiae and paradoxus
# references; align each chromosome separately
#python align/run_mugsy.sh

# predicted introgressed regions for each chromosome of each strain
sh analyze/run_analyze.sh

# extract alignments of introgressed regions and annotate genes in
# those regions
python analyze/process.py

# find predicted introgressed genes that are the same/different between
# 100-genomes paper and my sets; also notate which sites match cer/par
# references or both/neither for each region alignment
python analyze/compare_all.py

# do the above but interactively for individual genes
python analyze/compare.py

