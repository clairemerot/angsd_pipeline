#!/bin/bash
#SBATCH -J "10_local_pca"
#SBATCH -o log_%j
#SBATCH -c 1 
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=7-00:00
#SBATCH --mem=30G

###this script will work on all bamfiles and calculate saf, maf & genotype likelihood
#maybe edit
NB_CPU=1 #change accordingly in SLURM header
window_size=10000 #nb of Snps per window

#the path to a file where there are informations to correlate with  such as PCs extracted from pca by LG, phenotypês, etc...
path_info=02_info/LG_pc_info.txt #if we want to look at correlation with a given variable
#path_info="FALSE" 

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR
ulimit -S -n 2048

#prepare variables - avoid to modify
source 01_scripts/01_config.sh
echo "run local pca on covariance matrix made by window with window size = $window_size"

mkdir "10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"/$window_size/analyse_by_window"

ls 10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"/"$window_size"/cov_by_window/ > 10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"/"$window_size"/analyse_by_window/list_window.txt

#the last is the path to a file where there are informations to correlate with  such as PCs extracted from pca by LG, phenotypês, etc...

Rscript 01_scripts/Rscripts/local_pca.R "$MIN_MAF" "$PERCENT_IND" "$window_size" "$path_info"

