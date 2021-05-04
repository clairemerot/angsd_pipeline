#!/bin/bash
#SBATCH -J "10A_cov_by_window"
#SBATCH -o log_%j
#SBATCH -c 4 
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=XXX
#SBATCH --time=7-00:00
#SBATCH --mem=1G

###this script will work on all bamfiles and calculate saf, maf & genotype likelihood
#maybe edit
NB_CPU=4 #change accordingly in SLURM header
#window_size=760000 #nb of Snps per window - this is for a window of the size of my biggest chromosome
window_size=1000 # for small window along a chromosome for instance
# adjust memory accordingly - small window = few memory - big window more memory (1G for 1000 snp - 50G for a chromosome

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

ulimit -S -n 2048

#prepare variables - avoid to modify
source 01_scripts/01_config.sh
echo "split beagle into windows within each scaffold"

####split beagle into non overlapping windows
mkdir "10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR
mkdir "10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$window_size
mkdir "10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$window_size"/beagle_by_window"
mkdir "10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$window_size"/cov_by_window"

#if your reference contigs names are "Chr1" use this line
python 01_scripts/utility_scripts/beagle_sliding_window.py 03_saf_maf_gl_all/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".beagle.gz $window_size 10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$window_size"/beagle_by_window/

#if your reference contigs names are "Chr_1" use thi line (same but different python script)
#python 01_scripts/utility_scripts/beagle_sliding_window_bis.py 03_saf_maf_gl_all/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".beagle.gz $window_size 10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$window_size"/beagle_by_window/

#run pcangsd on each window
#this is the input file for the pca
ls -1 10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$window_size"/beagle_by_window/ | 
    sort -u | 
    while read i
    do
        echo "GL for $i"
INPUT=10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$window_size"/beagle_by_window/$i

echo "analyse covariance matrix on all individuals"
python2 $PCA_ANGSD_PATH/pcangsd.py -threads $NB_CPU -beagle $INPUT -o 10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$window_size"/cov_by_window/$i

done
