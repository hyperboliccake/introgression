# TODO extract all directories into some file so they can be set separately

git status --porcelain > uncommitted.txt
if [[ -s uncommitted.txt ]]
then
    echo 'uncommitted changes...exiting!'
    exit 1
fi
git rev-parse HEAD > sha.txt


# simulations and analysis
cd sim
fn = "run_sim_multi_model"
qsub -N $fn run_sim_multi_model.sh
# while grep returns something, i.e. job still running
while [qstat | grep -q $fn]
do
sleep 300
done
python sim_analyze_hmm_bw.py
python aggregate.py
Rscript plot.R
cp ../sha.txt ../../results/sim/
cd ..

# get alignments between each strain and the cerevisiae and paradoxus
# references; align each chromosome separately
cd align
fn = "run_mugsy"
# while grep returns something, i.e. job still running
qsub -N $fn run_mugsy.sh
while [qstat | grep -q $fn]
do
sleep 300
done
cp ../sha.txt ../../alignments/genbank/
cd ..

# predicted introgressed regions for each chromosome of each strain
# note: this requires ~12G memory
cd analyze
fn = "run_analyze"
# while grep returns something, i.e. job still running
qsub -N $fn run_analyze.sh
while [qstat | grep -q $fn]
do
sleep 300
done
cp ../sha.txt ../../results/
cd ..

# extract alignments of introgressed regions and annotate genes in
# those regions
cd analyze
python process.py
cp ../sha.txt ../../results/regions/
cd ..

# find predicted introgressed genes that are the same/different between
# 100-genomes paper and my sets; also notate which sites match cer/par
# references or both/neither for each region alignment
cd analyze
python compare_all.py
cp ../sha.txt ../../results/gene_alignments/
cd ..

# do the above but interactively for individual genes
cd analyze
python compare.py
cd ..
