######commands sparcc
#Correlation Calculation:
#First, we'll quantify the correlation between all OTUs, using SparCC, Pearson, and Spearman correlations:

#!/bin/bash

#PBS -l nodes=1:ppn=12
#PBS -l walltime=05:00:00
#PBS -N sparcc12h
cd $PBS_O_WORKDIR

module load devel/perl/5.26
module load lib/anaconda3/4.2.0
#conda create --name sparcc python=2.6.9 numpy=1.9.2 pandas=0.16.2 scipy
source  activate sparcc

# Open the Python environment for sparcc:
#source activate Networking

# Make the permutation tables:
# Need more permutations
mkdir file_Permutations_1000
mkdir Sparcc_Default_1000

python ../commands/MakeBootstraps.py file.txt -n 1000 -t permutation_#.txt -p file_Permutations_1000/

perl MakeRunScript_1000.pl
chmod 755 RunParallel_1000.s
nohup sh RunParallel_1000.s

# Summarize p-values
python ../commands/PseudoPvals.py ./Sparcc_Default_1000/file_corsparcc.txt ./Sparcc_Default_1000/perm_cor_#.txt 1000 -o ./Sparcc_Default_1000/pvals.two_sided.txt -t two_sided
#python commands/PseudoPvals.py ./Sparcc_Pearson_1000/Alltables_Trimmed_Combined_corsparcc.txt ./Sparcc_Pearson_1000/perm_cor_#.txt 1000 -o ./Sparcc_Pearson_1000/pvals.two_sided.txt -t two_sided
#python commands/PseudoPvals.py ./Sparcc_Spearman_1000/Alltables_Trimmed_Combined_corsparcc.txt ./Sparcc_Spearman_1000/perm_cor_#.txt 1000 -o ./Sparcc_Spearman_1000/pvals.two_sided.txt -t two_sided

# cleanup
#perl Clean_Reformat_SparCCNetwork_1000.pl
