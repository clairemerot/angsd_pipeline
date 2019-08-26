#!/bin/bash
#SBATCH -J "10_pca_by_window"
#SBATCH -o log_%j
#SBATCH -c 4 
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=7-00:00
#SBATCH --mem=30G

###this script will work on all bamfiles and calculate saf, maf & genotype likelihood
#maybe edit
NB_CPU=4 #change accordingly in SLURM header
window_size=10000 #nb of Snps per window

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

ulimit -S -n 2048

#prepare variables - avoid to modify
source 01_scripts/01_config.sh
echo "split beagle into windows within each scaffold"

####split beagle into non overlapping windows
mkdir "10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND
mkdir "10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"/$window_size"
mkdir "10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"/$window_size/beagle_by_window"
mkdir "10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"/$window_size/cov_by_window"

python 01_scripts/utility_scripts/beagle_sliding_window.py 03_saf_maf_gl_all/all_maf"$MIN_MAF"_pctind"$PERCENT_IND".beagle.gz $window_size 10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"/"$window_size"/beagle_by_window/

#run pcangsd on each window
#this is the input file for the pca
ls -1 10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"/"$window_size"/beagle_by_window/ | 
    sort -u | 
    while read i
    do
        echo "GL for $i"
INPUT=10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"/"$window_size"/beagle_by_window/$i

echo "analyse covariance matrix on all individuals"
python2 $PCA_ANGSD_PATH/pcangsd.py -threads $NB_CPU -beagle $INPUT -o 10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"/"$window_size"/cov_by_window/$i

done
