#!/bin/bash
#SBATCH -J "09_gwas_bin"
#SBATCH -o log_%j
#SBATCH -c 3
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=7-00:00
#SBATCH --mem=150G

NB_CPU=3 #change accordingly in SLURM header
PHENO=02_info/bin_pheno.txt #this file must be one single column with phenotype coded as 1 or 2, each line is one individual in the same order as bamfile
REGIONS="-rf 02_info/regions_25kb_100snp.txt" #optional edit with your region selected file

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR
module load angsd
ulimit -S -n 2048

#prepare variables - avoid to modify
source 01_scripts/01_config.sh
N_IND=$(wc -l 02_info/bam.filelist | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)


# perform gwas on all samples
angsd -P $NB_CPU -nQueueSize 50 -ybin $PHENO -doAsso 2 -GL 2 \
-doPost 1 -doMajorMinor 1 -doMaf 1 -doCounts 1 \
-remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -minMaf $MIN_MAF -setMaxDepth $MAX_DEPTH -setMinDepthInd $MIN_DEPTH \
$REGIONS -b 02_info/bam.filelist -out 09_gwas/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".binary.gwas



