#!/bin/bash
#SBATCH -J "12A_prepngsLD_2LG"
#SBATCH -o log_%j
#SBATCH -c 6
#SBATCH -p large
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=21-00:00
#SBATCH --mem=30G

###this script calculate LD within and between 2 chromosomes
#maybe edit
NB_CPU=6 #change accordingly in SLURM header

LG_FILE1=02_info/chrs_fusion.txt #work on a chosen subset of LG 
#this is a text file with the name of the chromosomes. 1 chr by line. 
#Try to do a maximum of 2 chromosomes as this scrip will calculate intra-chr and inter-chr LD

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#module
module load angsd
ulimit -S -n 2048

#we need more stringent filter to assure good data for LD but not too many SNPs.
MIN_MAF=0.10 #filter : will keep SNP above this allele frequency (over all individuals)
PERCENT_IND=0.90 #filter : will keep SNP with at least one read for this percentage of individuals 
MAX_DEPTH_FACTOR=3 #filter : will keep SNP with less than X time the nb of individuals total coverage
MIN_DEPTH=2


N_IND=$(wc -l 02_info/bam.filelist | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)



#make a beagle and filter for maf an min ind
angsd -P $NB_CPU \
-GL 2 -doGlf 2 -doMajorMinor 3 -doMaf 1 -doCounts 1 \
-anc 02_info/genome.fasta  -minMaf $MIN_MAF \
-remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -setMinDepthInd $MIN_DEPTH -setMaxDepth $MAX_DEPTH \
-rf 02_info/chrs_fusion.txt \
-b 02_info/bam.filelist -out 12_ngsLD/all

##input position file
gunzip -c 12_ngsLD/all.mafs.gz | cut -f -2 | awk '{if(NR>1)print}' > 12_ngsLD/sites
POS_FILE=12_ngsLD/sites
#calculate the number of sites
N_SITES=$(wc -l $POS_FILE | cut -d " " -f 1)
N_IND=$(wc -l 02_info/bam.filelist | cut -d " " -f 1)

echo "beagle successfully built, it includes $N_SITES sites"
echo "if this number is too big (> 100,000) re-run with more stringent filters"
