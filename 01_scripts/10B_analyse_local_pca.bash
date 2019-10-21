#!/bin/bash
#SBATCH -J "10B_plot_pca_by_window"
#SBATCH -o log_%j
#SBATCH -c 1 
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=7-00:00
#SBATCH --mem=30G

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#A script to analyse local pca done on scaffold or window
#will output a figure by scaffold/by window of local pca with colours coming form an info file
#will output correlation with the main PCs
# to edit
#which window size was used for pca by scaffold?
WINDOW_SIZE=760000
#make sure to have factors on which we want to colour the points
PATH_INFO_FILE=02_info/1446_info.txt
#which column on the file do we want to use as variable?
var1=9
var2=10
var3=3
var4=4


#prepare other variables - avoid to modify
source 01_scripts/01_config.sh




mkdir 10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$WINDOW_SIZE"/pdf_by_window

#make the list of cov files on which we want to do pca
ls -1 10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$WINDOW_SIZE"/cov_by_window > 10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$WINDOW_SIZE"/pdf_by_window/window_cov_files.txt
#edit that file to remove small scaffolds?

#run Rscript
Rscript 01_scripts/Rscripts/visu_local_pca.r "$MIN_MAF" "$PERCENT_IND" "$WINDOW_SIZE" "$PATH_INFO_FILE" "$MAX_DEPTH_FACTOR" "$var1" "$var2" "$var3" "$var4"
