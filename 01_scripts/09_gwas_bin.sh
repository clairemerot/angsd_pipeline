#!/bin/bash
#SBATCH -J "09_gwas_bin"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=3-00:00
#SBATCH --mem=40960

NB_CPU=1 #change accordingly in SLURM header
PHENO=02_info/bin_pheno.txt #this file must be one single column with phenotype coded as 1 or 2, each line is one individual in the same order as bamfile

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#prepare variables - avoid to modify
source 01_scripts/01_config.sh
N_IND=$(wc -l 02_info/bam.filelist | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 

# perform gwas on all samples
angsd -P $NB_CPU -nQueueSize 50 -ybin $PHENO -doAsso 2 -GL 2 -doPost 1 -doMajorMinor 1 -doMaf 1 -remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -minMaf $MIN_MAF -b 02_info/bam.filelist -out 09_gwas/all_maf"$MIN_MAF"_pctind"$PERCENT_IND".binary.gwas



