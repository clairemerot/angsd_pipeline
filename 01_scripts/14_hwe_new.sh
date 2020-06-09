#!/bin/bash
#SBATCH -J "14_hwe_by_group"
#SBATCH -o log_%j
#SBATCH -c 6
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=7-00:00
#SBATCH --mem=80G

###this script will work on bamfiles by population and calculate saf  & maf 
#maybe edit
NB_CPU=6 #change accordingly in SLURM header
POP_FILE1=02_info/sex.txt #choose on which list of pop run the analyses

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR
module load angsd
ulimit -S -n 2048

#prepare variables - avoid to modify
source 01_scripts/01_config.sh

# Do hwe for all population listed
cat $POP_FILE1 | while read i
do
echo $i

angsd -P $NB_CPU -nQueueSize 50 \
-doHWE 1 -GL 2 -remove_bads 1 -minMapQ 30 -minQ 20 -doMajorMinor 3 \
-sites 02_info/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR" \
-b 02_info/"$i"bam.filelist -out 14_HWE/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"



#Rscript to infer Hobs and do sliding window output

gunzip 14_HWE/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".hwe.gz

FILE=14_HWE/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".hwe

#window size and window step are chosen in the 01_config or can be edited here
#WINDOW=25000
#WINDOW_STEP=5000
Rscript 01_scripts/Rscripts/Hobs_sliding.r "$FILE" "$WINDOW" "$WINDOW_STEP"

gzip $FILE
gzip "$FILE".Hobs
done
