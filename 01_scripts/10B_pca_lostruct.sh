#!/bin/bash
#SBATCH -J "10B_lostruct_pca_all_LG"
#SBATCH -o log_%j
#SBATCH -c 1 
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=XX
#SBATCH --time=7-00:00
#SBATCH --mem=200G

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#A script to analyse local pca done on many windows with all LG/scaffodl together
#will use the lostruct method to make a mds
#will output figures of mds and mds scores per window

#### to edit ###
window_size=100 #nb of Snps per window
N_IND=1446 #nb of individuals included in the analysis
N_PC=2 #nb of PC to consider when comparing windows
N_MDS=50 #nb of MDS dimension in the output

#prepare variables - avoid to modify
source 01_scripts/01_config.sh

echo "run local pca on covariance matrix made by window with window size = $window_size"

mkdir "10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/$window_size/analyse_by_window"
ls 10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$window_size"/cov_by_window/ > 10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$window_size"/analyse_by_window/list_window.txt

#if cov files have been zipped
ls -1 10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$window_size"/cov_by_window/*.gz | 
    while read i
    do
	gunzip $i
done

Rscript 01_scripts/Rscripts/local_pca_mds_all_chr_maxdepth_speedy_version.R "$MIN_MAF" "$PERCENT_IND" "$window_size" "$MAX_DEPTH_FACTOR" "$N_IND" "$N_PC" "$N_MDS"
