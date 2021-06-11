#!/bin/bash
#SBATCH -J "12A_ngsLDprep_byLG"
#SBATCH -o log_%j
#SBATCH -c 4
#SBATCH -p large
#SBATCH --mail-type=ALL
#SBATCH --mail-user=X
#SBATCH --time=21-00:00
#SBATCH --mem=100G

###this script will work on all bamfiles and output a reduced beagle for each chromosome in preparation for LD calculation
#maybe edit
NB_CPU=4 #change accordingly in SLURM header

LG_FILE1=02_info/region_LG.txt #work on a chosen subset of LG

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#module
module load angsd
ulimit -S -n 2048

#we use more stringent filter as usually to be more confident in the SNPs and have a smaller number of sites
MIN_MAF=0.10 #filter : will keep SNP above this allele frequency (over all individuals)
PERCENT_IND=0.75 #filter : will keep SNP with at least one read for this percentage of individuals 
MAX_DEPTH_FACTOR=3 #filter : will keep SNP with less than X time the nb of individuals total coverage



N_IND=$(wc -l 02_info/bam.filelist | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

#loop over LG to output a beagle on which we will run ngsLD


cat $LG_FILE1 | while read i
do
echo $i
mkdir 12_ngsLD/$i

#make a beagle and filter for maf an min ind
angsd -P $NB_CPU -nQueueSize 100 -underFlowProtect 1 \
-dosaf 1 -GL 2 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 \
-anc 02_info/genome.fasta -r $i \
-remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -setMaxDepth $MAX_DEPTH -minMaf $MIN_MAF -setMinDepthInd $MIN_DEPTH \
-b 02_info/bam.filelist -out 12_ngsLD/$i/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$i"

#input position file
gunzip -c 12_ngsLD/$i/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$i".mafs.gz | cut -f -2 | awk '{if(NR>1)print}' > 12_ngsLD/$i/sites_maf"$MIN_MAF"_pctind_"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$i" 
POS_FILE=12_ngsLD/$i/sites_maf"$MIN_MAF"_pctind_"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$i"
#calculate the number of sites
N_SITES=$(wc -l $POS_FILE | cut -d " " -f 1)

echo "the beagle include " $N_SITES "sites for" $i
echo "if this number is to big (>100000-200000à, re-run with more stringent filters

done
