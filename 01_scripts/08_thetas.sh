#!/bin/bash
#SBATCH -J "08_thetas"
#SBATCH -o log_%j
#SBATCH -c 4
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claier.merot@gmail.com
#SBATCH --time=7-00:00
#SBATCH --mem=80G

###this script will work on bamfiles by population and calculate saf  & maf 
#maybe edit
NB_CPU=4 #change accordingly in SLURM header
POP_FILE1=02_info/pop2.txt #choose on which list of pop run the analyses

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR
module load angsd
ulimit -S -n 2048

#prepare variables - avoid to modify
source 01_scripts/01_config.sh


#Do thetas on all samples
N_IND=$(wc -l 02_info/bam.filelist | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 

echo "estimate thetas for all samples"
#we need to re-do doSaf because we don't want to filter on maf for thetas calculation
#I don't use the fold option anymore - but be aware that only T watterson and Taj D are interpretable if anc is the ref genome
#angsd -P $NB_CPU -nQueueSize 50 \
#-dosaf 1 -GL 2 -doMajorMinor 1 \
#-anc 02_info/genome.fasta \
#-remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND \
#-b 02_info/bam.filelist -out 08_thetas/all_pctind"$PERCENT_IND"

#realSFS 08_thetas/all_pctind"$PERCENT_IND".saf.idx -P $NB_CPU > 08_thetas/all_pctind"$PERCENT_IND".sfs

#angsd -P $NB_CPU -nQueueSize 50 -dosaf 1 -doThetas 1 -GL 2 -doMajorMinor 3 -anc 02_info/genome.fasta \
#-remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -b 02_info/bam.filelist \
#-pest 08_thetas/all_pctind"$PERCENT_IND".sfs -out 08_thetas/all_pctind"$PERCENT_IND"

#Estimate for every Chromosome/scaffold
#thetaStat do_stat 08_thetas/all_pctind"$PERCENT_IND".thetas.idx

#Do a sliding window analysis based on the output from the make_bed command.
#thetaStat do_stat 08_thetas/all_pctind"$PERCENT_IND".thetas.idx -win $WINDOW -step $WINDOW_STEP \
-outnames 08_thetas/all_pctind"$PERCENT_IND".thetaswindow



# Do thetas for all population listed
cat $POP_FILE1 | while read i
do
	echo $i
	N_IND=$(wc -l 02_info/"$i"bam.filelist | cut -d " " -f 1)
	MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
	MIN_IND=${MIN_IND_FLOAT%.*} 

	#unsure whether we should specify the sites & filter for the thetas too??
	echo "working on pop $i, $N_IND individuals"
	echo "will filter for sites with at least one read in $MIN_IND individuals, which is $PERCENT_IND of the total"
	
	#we need to re-do doSaf because we don't want to filter on maf for thetas calculation
	#I don't use the fold option anymore - but be aware that only T watterson and Taj D are interpretable if anc is the ref genome
	angsd -P $NB_CPU -nQueueSize 50 \
	-dosaf 1 -GL 2 -doMajorMinor 1 \
	-anc 02_info/genome.fasta \
	-remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND \
	-b 02_info/"$i"bam.filelist -out 08_thetas/"$i"_pctind"$PERCENT_IND"
	
	echo "estimate real sfs for pop $i"
	realSFS 08_thetas/"$i"_pctind"$PERCENT_IND".saf.idx -P $NB_CPU > 08_thetas/"$i"_pctind"$PERCENT_IND".sfs
	
	echo "estimate thetas for pop $i"
	angsd -P $NB_CPU -nQueueSize 50 -dosaf 1 -doThetas 1 -GL 2 -doMajorMinor 3 \
	-anc 02_info/genome.fasta -remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND \
	-b 02_info/"$i"bam.filelist \
	-pest 08_thetas/"$i"_pctind"$PERCENT_IND".sfs -out 08_thetas/"$i"_pctind"$PERCENT_IND"
	
	#Estimate for every Chromosome/scaffold
	thetaStat do_stat 08_thetas/"$i"_pctind"$PERCENT_IND".thetas.idx
	#Do a sliding window analysis based on the output from the make_bed command.
	thetaStat do_stat 08_thetas/"$i"_pctind"$PERCENT_IND".thetas.idx -win $WINDOW -step $WINDOW_STEP \
	-outnames 08_thetas/"$i"_pctind"$PERCENT_IND".thetaswindow
done
