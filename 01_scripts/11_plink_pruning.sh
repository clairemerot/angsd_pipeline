#!/bin/bash
#SBATCH -J "plink_all"
#SBATCH -o log_%j
#SBATCH -c 4
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=7-00:00
#SBATCH --mem=200G

###this script will work on all bamfiles and output a plink file for each chromosome
#maybe edit
NB_CPU=4 #change accordingly in SLURM header
REGIONS="-rf 02_info/regions_25kb_100snp.txt" 
#REGIONS="-rf 02_info/chromosomes.txt" 

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

module load angsd
module load plink
ulimit -S -n 2048

#prepare variables - avoid to modify
source 01_scripts/01_config.sh

N_IND=$(wc -l 02_info/bam.filelist | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

echo " output genotype in plink format for all individuals listed in 02_info/bam.filelist"
echo "keep loci with at leat one read for n individuals = $MIN_IND, which is $PERCENT_IND % of total $N_IND individuals"
echo "filter on allele frequency = $MIN_MAF"


####Calculate the SAF, MAF and GL

angsd -P $NB_CPU -nQueueSize 50 -doPlink 2 -doMaf 1 -doCounts 1 \
-GL 2 -doMajorMinor 1 -doGeno -4 -doPost 1 -postCutoff 0.8 \
-anc 02_info/genome.fasta -remove_bads 1 -minMapQ 30 -minQ 20 $REGIONS \
-minInd $MIN_IND -minMaf $MIN_MAF -setMaxDepth $MAX_DEPTH -b 02_info/bam.filelist \
-out 11_plink/all_plink_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"

echo "provide a list of LD-pruned snp"
WINDOW=100
SNP=5
VIF=2

plink --tped 11_plink/all_plink_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".tped \
--tfam 11_plink/all_plink_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".tfam --indep $WINDOW $SNP $VIF --allow-extra-chr \
--out 11_plink/all_plink_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".pruned --threads $NB_CPU 


