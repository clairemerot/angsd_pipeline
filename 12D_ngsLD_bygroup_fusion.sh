#!/bin/bash
#SBATCH -J "12D_ngsLD_bygroup"
#SBATCH -o log_%j
#SBATCH -c 10
#SBATCH -p large
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=21-00:00
#SBATCH --mem=30G

###this script will take a beagle by chromosome and cut it by subgroups of equal size to calculate LD on the same genotype original file

#maybe edit
NB_CPU=10 #change accordingly in SLURM header
GROUP=fusion_homo #the subgroup on whcih we are making the fst comparison -> it should be a file like GROUP.txt in the folder 02_info
POP_FILE1=02_info/"$GROUP".txt #choose on which list of pop run the analyses

#do not edit if 12A has been run for each LG under the same filters. if not, please run before 
#because we are here re-using the beagle produced with all samples, on the given LG

BAM_ALL=02_info/bam.filelist # where and what is the name of the big bam list on which the beagle has been built


# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR
#module
module load angsd
ulimit -S -n 2048


#make a folder in which write new results
mkdir 12_ngsLD/$GROUP

#1 subset bamfilelist
#and make a file with the position of the individuals for each subgroup in the bigger beagle
Rscript 01_scripts/Rscripts/subset_random_Nind_12ngsLD.r "$GROUP" "$BAM_ALL"


#2 decompresss the beagle produced by the script 12A on the whole dataset for the given LG
gunzip -k 12_ngsLD/all.beagle.gz
#get the position file
POS_FILE=12_ngsLD/sites

#calculate the number of sites
#check this number. Do not run if > 200 000. it will be long and make big files if around 100 000
N_SITES=$(wc -l $POS_FILE | cut -d " " -f 1)
	


#3 loop on populations to make sub-beagle, calculate LD & plot it for each subgroup
cat $POP_FILE1 | while read i
do
echo $i

N_IND=$(wc -l 12_ngsLD/"$GROUP"/"$i"subsetbam.filelist | cut -d " " -f 1)

#we are not running anymore angsd again because this was calling a different set of snp and diff MAj/min alleles
#but this would still be possible using the subsetbam.filelist
#just trying to subset the global beagle
	
	cat 12_ngsLD/"$GROUP"/"$i"_column.pos | while read k 
	do 
	echo $k
	cat 12_ngsLD/all.beagle | cut -f $k > 12_ngsLD/$GROUP/"$i".beagle
	gzip 12_ngsLD/$GROUP/"$i".beagle
	done
	
	
	
echo "calculating LD on " $i
#run ngsLD
ngsLD --geno 12_ngsLD/$GROUP/"$i".beagle.gz --probs --n_ind $N_IND --n_sites $N_SITES --pos $POS_FILE \
--n_threads $NB_CPU --max_kb_dist 0 --max_snp_dist 0 --out 12_ngsLD/$GROUP/"$i".ld --extend_out

gzip 12_ngsLD/$GROUP/"$i".ld

done

# LD files can be processed as done with the 12C file to simplify by windows of 500kb