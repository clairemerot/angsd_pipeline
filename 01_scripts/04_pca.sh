#!/bin/bash
#SBATCH -J "04_pca"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=X
#SBATCH --time=7-00:00
#SBATCH --mem=200G

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR
###this script will work on all individuals using the beagle genotype likelihood and calculate a covariance matrix with angsd & a pca with R
#this requires pcangsd to be cloned and a version of Python v2 with alias python2

#maybe edit
NB_CPU=1 #change accordingly in SLURM header

##load the adequate python environment that you have created (1st intitialize conda with YOUR PATH to miniconda, except if you already have it in your .bashrc)

##activate a conda environnement in which you have install the necessary library and python version
#conda create --name pcangsd_test python=2.7 
#conda activate pcangsd_test 
#conda install ipython numpy scipy six pandas python-dateutil cython numba

#you may need to edit the name of the environnemnt depending on what you chose
conda activate pcangsd_test 

#prepare variables - avoid to modify
source 01_scripts/01_config.sh

#this is the list of file we are working on
BAM_LIST=02_info/bam.filelist

#this is the input file for the pca
INPUT=03_saf_maf_gl_all/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".beagle.gz

echo "analyse covariance matrix on all individuals"
python2 $PCA_ANGSD_PATH/pcangsd.py -threads $NB_CPU \
	-beagle $INPUT -o 04_pca/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"

echo "transform covariance matrix into PCA"
COV_MAT=04_pca/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".cov
Rscript 01_scripts/Rscripts/make_pca_simple.r "$COV_MAT" "$BAM_LIST"

