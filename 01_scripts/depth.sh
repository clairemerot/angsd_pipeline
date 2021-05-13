#!/bin/bash
#SBATCH -J "depth"
#SBATCH -o log_%j
#SBATCH -c 5 
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=7-00:00
#SBATCH --mem=50G

###this script will work on all bamfiles and calculate saf, maf & genotype likelihood
#maybe edit
NB_CPU=5 #change accordingly in SLURM header
#REGIONS="-rf 02_info/regions_25kb_100snp.txt" #optional edit with your region selected file
#REGIONS="-r LG6" # to remove the options to focus on a limited number of regions

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

module load angsd
ulimit -S -n 2048

#prepare variables - avoid to modify
source 01_scripts/01_config.sh


#if you want to have depth on position after some filtering:
N_IND=$(wc -l 02_info/bam.filelist | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

echo " Calculate the SAF, MAF and GL for all individuals listed in 02_info/bam.filelist"
echo "keep loci with at leat one read for n individuals = $MIN_IND, which is $PERCENT_IND % of total $N_IND individuals"
echo "filter on allele frequency = $MIN_MAF"

####Calculate the SAF, MAF and GL and export depth at the same time
angsd -P $NB_CPU -nQueueSize 50 \
-doCounts 1 \
-doDepth 1 -maxDepth 1000 -dumpCounts 2 \
-anc 02_info/genome.fasta -remove_bads 1 -minMapQ 30 -minQ 20 \
-b 02_info/bam.filelist \
$REGIONS -out 03_saf_maf_gl_all/depth

#if you want to have depth on position after some filtering include some of the following filters
#-minInd $MIN_IND -minMaf $MIN_MAF -setMaxDepth $MAX_DEPTH -setMinDepthInd $MIN_DEPTH

##the output is hard to understand. please refer to the manual
#http://www.popgen.dk/angsd/index.php/Depth