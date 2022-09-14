#!/bin/bash
#SBATCH -J "07_FST_by_group"
#SBATCH -o log_%j
#SBATCH -c 6
#SBATCH -p large
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=21-00:00
#SBATCH --mem=50G

###this script will 
#1 make a subset of bamlist to ahev equal nb of samples
#2calculate the saf by pop
#3 calculate 2dSFS and FST
#it uses the last version of angsd (unfold for saf and fold for realSFS)

#maybe edit
NB_CPU=6 #change accordingly in SLURM header
NSITES=500000 #to make realSFS goes faster -reduce the number of sites considered
GROUP=pop #the subgroup on whcih we are making the fst comparison -> it should be a file like GROUP.txt in the folder 02_info
POP_FILE1=02_info/"$GROUP".txt #choose on which list of pop run the analyses

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR
module load angsd/0.931 #only with this version the SFS/FST script runs well (edit in sept 2022)
#make sure you index the sites file with the same version 
ulimit -S -n 2048

#prepare variables - avoid to modify
source 01_scripts/01_config.sh

#make a folder in which write new results
mkdir 07_fst_by_pop_pair/$GROUP

#1 subset bamfilelist
Rscript 01_scripts/Rscripts/subset_random_Nind.r "$GROUP"

#2 do saf for all population listed
cat $POP_FILE1 | while read i
do
echo $i

N_IND=$(wc -l 07_fst_by_pop_pair/$GROUP/"$i"subsetbam.filelist | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 

echo "working on pop $i, $N_IND individuals, will use the sites file provided"
echo "will filter for sites with at least one read in $MIN_IND individuals, which is $PERCENT_IND of the total"

angsd -P $NB_CPU \
-dosaf 1 -GL 2 -doMajorMinor 3 \
-anc 02_info/genome.fasta \
-remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -setMinDepthInd $MIN_DEPTH \
-sites 02_info/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR" \
-rf 02_info/regions_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR" \
-b 07_fst_by_pop_pair/$GROUP/"$i"subsetbam.filelist -out 07_fst_by_pop_pair/$GROUP/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"
done



# 3 calculate FST
#prepare variables - avoid to modify
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
		realSFS  07_fst_by_pop_pair/$GROUP/"$pop1"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".saf.idx \
07_fst_by_pop_pair/$GROUP/"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".saf.idx \
-P $NB_CPU -maxIter 30 -nSites $NSITES > 07_fst_by_pop_pair/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"."$NSITES"

file=07_fst_by_pop_pair/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"."$NSITES"

Rscript 01_scripts/Rscripts/sum_sites_2dsfs.r "$file"
		
		echo " prepare the fst for easy window analysis etc"
		realSFS fst index 07_fst_by_pop_pair/$GROUP/"$pop1"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".saf.idx \
07_fst_by_pop_pair/$GROUP/"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".saf.idx \
-sfs 07_fst_by_pop_pair/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"."$NSITES".2dsfs \
-P $NB_CPU -fstout 07_fst_by_pop_pair/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"

		echo "print SFS priori for each position"
		realSFS fst print 07_fst_by_pop_pair/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".fst.idx \
-P $NB_CPU > 07_fst_by_pop_pair/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".bypos.sfs
		
		echo "get the global estimate of FST throughout the genome"
		realSFS fst stats 07_fst_by_pop_pair/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".fst.idx \
-P $NB_CPU > 07_fst_by_pop_pair/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".fst
		
		echo "calculate FST by slidingwindow, window size=$WINDOW and step=$WINDOW_STEP, as given in 01_config.sh"
		realSFS  fst stats2 07_fst_by_pop_pair/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".fst.idx \
-win $WINDOW -step $WINDOW_STEP -P $NB_CPU > 07_fst_by_pop_pair/$GROUP/"$pop1"_"$pop2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".slidingwindow
	done
done
