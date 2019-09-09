#!/bin/bash
#SBATCH -J "07_FST_by_pop_pair"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=3-00:00
#SBATCH --mem=50G

###this script will use the saf by population calculated at step 07 and calculate SFS and FST
#maybe edit
NB_CPU=1 #change accordingly in SLURM header
POP_FILE1=02_info/pop.txt #choose on which list of pop run the analyses


# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR
module load angsd
ulimit -S -n 2048

#prepare variables - avoid to modify
source 01_scripts/01_config.sh
num_pops=$(wc -l "$POP_FILE1" | cut -d " " -f 1)

# Estimate pairwise FST for all populations listed

for i in $(seq $num_pops)
do
	pop1=$(cat "$POP_FILE1" | head -"$i" | tail -1)
	for j in $(seq $[ $i + 1 ] $num_pops)
	do
		pop2=$(cat "$POP_FILE1" | head -"$j" | tail -1)
		echo "FST between $pop1 and $pop2"
		echo "$pop1"
		echo "$pop2"
		
		echo "calcualte the 2dsfs priors"
		realSFS  06_saf_maf_by_pop/$pop1/"$pop1"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".saf.idx \
06_saf_maf_by_pop/$pop2/"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".saf.idx \
-P $NB_CPU > 07_fst_by_pop_pair/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".ml
		
		echo " prepare the fst for easy window analysis etc"
		realSFS fst index 06_saf_maf_by_pop/$pop1/"$pop1"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".saf.idx \
06_saf_maf_by_pop/$pop2/"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".saf.idx \
-sfs 07_fst_by_pop_pair/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".ml \
-fstout 07_fst_by_pop_pair/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"

		echo "print SFS priori for each position"
		realSFS fst print 07_fst_by_pop_pair/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".fst.idx \
-P $NB_CPU > 07_fst_by_pop_pair/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".sfs
		
		echo "get the global estimate of FST throughout the genome"
		realSFS fst stats 07_fst_by_pop_pair/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".fst.idx \
-P $NB_CPU > 07_fst_by_pop_pair/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".fst
		
		echo "calculate FST by slidingwindow, window size=$WINDOW and step=$WINDOW_STEP, as given in 01_config.sh"
		realSFS  fst stats2 07_fst_by_pop_pair/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".fst.idx \
-win $WINDOW -step $WINDOW_STEP -P $NB_CPU > 07_fst_by_pop_pair/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".slidingwindow
	done
done
