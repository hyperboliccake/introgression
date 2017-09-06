# 
# Use the bash shell to interpret this job script 
#$ -S /bin/bash -t 1-31 -tc 31 -l m_mem_free=8G
# 
# Send an e-mail to the address 
# specified in .sge_request when this job ends. 
##$ -m e 
# 
# Only submit this job to nodes 
# that have at least 8GB of RAM free.
##$ -l mem_requested=8G, h_vmem=8G

# Put the hostname, current directory, and start date 
# into variables, then write them to the SGE standard output file.
# Having this information in your output file will help track down
# any errors which might occur during the job's run.
WRKHOST=`/bin/hostname` 
WRKDIR=`/bin/pwd` 
STRTDATE=`/bin/date` 
echo "**** JOB STARTED ON $WRKHOST AT $STRTDATE" 
echo "**** JOB RUNNING IN $WRKDIR" 

# SGE jobs do not run in your login environment, so you'll need to 
# load your environment, or atleast the modules your script needs to run
#source /net/gs/vol1/home/aclark4/.bashrc
# also load any modules as needed
source /etc/profile.d/modules.sh
module load modules modules-init modules-gs
module load python/2.7.3
#module load numpy/latest
#module load scipy/latest

# Scrip or command(s) to run via SGE

cd /net/gs/vol1/home/aclark4/projects/introgression/code/sim
ARGS=$(head -n $SGE_TASK_ID sim_compare_args.txt | tail -n 1)
#TAG=$(echo $ARGS | awk -F" " '{print $1}') # lol i'm so bad at bash

########
# simulate sequences
########

#python sim_multi_model.py $ARGS

########
# convert to nucleotides
########

cd /net/gs/vol1/home/aclark4/projects/introgression/code/phylo-hmm

#python gen_sim_seqs.py $ARGS

########
# predict with phylo-hmm
########

#python gen_phylo_hmm_input_file.py $ARGS
#PHYLO_INPUT="autoinput_" + $TAG + ".txt"
#java -jar ~/software/phylo_hmm/phmm-0.1/dist/lib/phmm.jar < $PHYLO_INPUT
python sim_analyze_phylo_main.py $ARGS

########
# predict with my method
########

#cd /net/gs/vol1/home/aclark4/projects/introgression/code/sim
#python sim_analyze_hmm_bw_main.py $ARGS
#Rscript plot.R

########
# compare methods
########

# NOTE THAT THE TWO METHODS DON'T QUITE OUTPUT THE SAME STUFF IN THE SUMMARY FILES-PHYLO METHOD ONLY GIVES RESULTS FOR ONE INDIVIDUAL, MY METHOD FOR ALL

#cd /net/gs/vol1/home/aclark4/projects/introgression/code/phylo-hmm
#python compare_methods.py $ARGS
#Rscript plot_comparison.R


