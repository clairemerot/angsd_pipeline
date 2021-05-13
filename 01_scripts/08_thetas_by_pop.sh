#!/bin/bash
#SBATCH -J "08_thetas_bypop"
#SBATCH -o log_%j
#SBATCH -c 10
#SBATCH -p large
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=21-00:00
#SBATCH --mem=80G

###this script will work on bamfiles by population and calculate saf  & maf 
#maybe edit
NB_CPU=10 #change accordingly in SLURM header
NSITES=5000000 #to make realSFS goes faster -reduce the number of sites considered
GROUP=pop #the subgroup on whcih we are making the fst comparison -> it should be a file like GROUP.txt in the folder 02_info
POP_FILE1=02_info/"$GROUP".txt #choose on which list of pop run the analyses


# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR
module load angsd
ulimit -S -n 2048

#prepare variables - avoid to modify
source 01_scripts/01_config.sh


# Do thetas for all population listed
#note that as bamlist I used teh ones with similar sample size drawn at the step 07 when doing the Fst.
#if samplesize differences is not a pb, one can replace 07_fst_by_pop_pair/$GROUP/"$i"subsetbam.filelist by 02_info/"$i"bam.filelist
cat $POP_FILE1 | while read i
do
	echo $i
	N_IND=$(wc -l 07_fst_by_pop_pair/$GROUP/"$i"subsetbam.filelist | cut -d " " -f 1)
	MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
	MIN_IND=${MIN_IND_FLOAT%.*} 
	MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

	#unsure whether we should specify the sites & filter for the thetas too??
	echo "working on pop $i, $N_IND individuals"
	echo "will filter for sites with at least one read in $MIN_IND individuals, which is $PERCENT_IND of the total"
	
	#we need to re-do doSaf because we don't want to filter on maf for thetas calculation
	#I don't use the fold option anymore - but be aware that only T watterson and Taj D are interpretable if anc is the ref genome
	angsd -P $NB_CPU -underFlowProtect 1 \
	-dosaf 1 -GL 2 -doMajorMinor 1 -doCounts 1 \
	-anc 02_info/genome.fasta \
	-remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -setMaxDepth $MAX_DEPTH -setMinDepthInd $MIN_DEPTH \
	-b 07_fst_by_pop_pair/$GROUP/"$i"subsetbam.filelist -out 08_thetas/"$i"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"
	
	echo "estimate real sfs for pop $i"
	realSFS 08_thetas/"$i"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".saf.idx -P $NB_CPU -nSites $NSITES  -maxIter 50 > 08_thetas/"$i"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"."$NSITES"
	file=08_thetas/"$i"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"."$NSITES"
	
	Rscript 01_scripts/Rscripts/sum_sites_sfs.r "$file"
	
	echo "estimate thetas for pop $i"
	angsd -P $NB_CPU -dosaf 1 -doThetas 1 -GL 2 -doMajorMinor 1 -underFlowProtect 1 \
	-anc 02_info/genome.fasta -remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -setMinDepthInd $MIN_DEPTH \
	-b 07_fst_by_pop_pair/$GROUP/"$i"subsetbam.filelist \
	-pest 08_thetas/"$i"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"."$NSITES".dsfs \
	-out 08_thetas/"$i"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"
	
	#Estimate for every Chromosome/scaffold
	thetaStat do_stat 08_thetas/"$i"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".thetas.idx
	#Do a sliding window analysis based on the output from the make_bed command.
	thetaStat do_stat 08_thetas/"$i"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".thetas.idx -win $WINDOW -step $WINDOW_STEP \
	-outnames 08_thetas/"$i"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".thetaswindow
done
