#!/bin/bash
#SBATCH -J "08_thetas_bypop"
#SBATCH -o log_%j
#SBATCH -c 10
#SBATCH -p large
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=21-00:00
#SBATCH --mem=80G

###this script will work on bamfiles by population and calculate saf  and then thetas
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
#it will loop on populations
	#to work on a signle population, comment the loop and use
	#i="popname"
cat $POP_FILE1 | while read i
do
	echo $i
	#note that as bamlist I used teh ones with similar sample size drawn at the step 07 when doing the Fst.
	#if samplesize differences is not a pb, one can replace 07_fst_by_pop_pair/$GROUP/"$i"subsetbam.filelist by 02_info/"$i"bam.filelist
	BAM_LIST=07_fst_by_pop_pair/$GROUP/"$i"subsetbam.filelist
	#BAM_LIST=02_info/"$i"bam.filelist
	
	#this will calculate filter for the initial saf
	N_IND=$(wc -l $BAM_LIST | cut -d " " -f 1)
	MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
	MIN_IND=${MIN_IND_FLOAT%.*} 
	MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

	#we do not filter on maf nor provide a site list
	#we need to re-do doSaf because we don't want to filter on maf for thetas calculation
	#I don't use the fold option anymore - but be aware that only T watterson and Taj D are interpretable if anc is the ref genome
	echo "working on pop $i, $N_IND individuals"
	echo "will filter for sites with at least one read in $MIN_IND individuals, which is $PERCENT_IND of the total and output saf"
		
	angsd -P $NB_CPU -underFlowProtect 1 \
	-dosaf 1 -GL 2 -doMajorMinor 1 -doCounts 1 \
	-anc 02_info/genome.fasta \
	-remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -setMaxDepth $MAX_DEPTH -setMinDepthInd $MIN_DEPTH \
	-b $BAM_LIST -out 08_thetas/"$i"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"
	
	#the output if a saf (for pop i)
	
	#now we use realSFS to calculate a 1dimension sfs
	#because it is long to calculate we ask to do it on sub section of the genome whose size is define by $NSITES"
	echo "estimate real sfs for pop $i by chunks of size $NSITES"
	realSFS 08_thetas/"$i"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".saf.idx -P $NB_CPU -nSites $NSITES  -maxIter 50 > 08_thetas/"$i"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"."$NSITES"
	
	#this gives a file with several lines while the expected sfs has just one line. So we use the R script to sum accross all columns
	#the R script will output a file whose name is the same but without the suffix "NSITES" (500000 in our defaut settings)
	file=08_thetas/"$i"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"."$NSITES"
	Rscript 01_scripts/Rscripts/sum_sites_sfs.r "$file"
	
	#now we use the sfs and the function dotheta to get the theta output
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
