#!/bin/bash
#SBATCH -J "06_saf_maf_by_pop"
#SBATCH -o log_%j
#SBATCH -c 4
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=7-00:00
#SBATCH --mem=30G

###this script will work on bamfiles by population and calculate saf  & maf 
#maybe edit
NB_CPU=4 #change accordingly in SLURM header
POP_FILE1=02_info/pop.txt #choose on which list of pop run the analyses

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR
module load angsd
ulimit -S -n 2048

#prepare variables - avoid to modify
source 01_scripts/01_config.sh

# Do saf/maf for all population listed
cat $POP_FILE1 | while read i
do
echo $i

mkdir 06_saf_maf_by_pop/$i

N_IND=$(wc -l 02_info/"$i"bam.filelist | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 

echo "working on pop $i, $N_IND individuals, will use the sites file provided"
echo "will filter for sites with at least one read in $MIN_IND individuals, which is $PERCENT_IND of the total"

angsd -P $NB_CPU -nQueueSize 50 \
-doMaf 1 -dosaf 1 -GL 2 -doMajorMinor 3 \
-anc 02_info/genome.fasta \
-remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -setMinDepthInd $MIN_DEPTH \
-sites 02_info/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR" \
-rf 02_info/regions_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR" \
-b 02_info/"$i"bam.filelist -out 06_saf_maf_by_pop/$i/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"

gunzip 06_saf_maf_by_pop/$i/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".mafs.gz
done

