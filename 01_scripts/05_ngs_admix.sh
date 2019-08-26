#!/bin/bash
#SBATCH -J "05_ngs_admix"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p large
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=21-00:00
#SBATCH --mem=200G

###this script will work on all individuals using the beagle genotype likelihood and perform an admixture analysis 
#this requires NGSadmix to be installed and its path export in the bashrc
#NGS admiw will explore all number of population between K_MIN and K_MAX as provided in the 01_config.sh

#maybe edit
NB_CPU=1 #change accordingly in SLURM header


# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#prepare variables - avoid to modify
source 01_scripts/01_config.sh


##boucle de MIN_K à MAX_K

for i in $(seq $K_MIN $K_MAX)
	do 
	echo $i
	NGSadmix -P $NB_CPU -likes 03_saf_maf_gl_all/all_maf"$MIN_MAF"_pctind"$PERCENT_IND".beagle.gz -minMaf $MIN_MAF -K $i -o 05_ngs_admix/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_K"$i"
	done
