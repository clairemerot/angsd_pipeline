#!/bin/bash
#SBATCH -J "ngsLD_bypop3"
#SBATCH -o log_%j
#SBATCH -c 3
#SBATCH -p large
#SBATCH --mail-type=ALL
#SBATCH --mail-user=X
#SBATCH --time=3-00:00
#SBATCH --mem=100G

###this script will take a beagle by chromosome and cut it by subgroups of equal size to calculate LD on the same genotype original file

#maybe edit
NB_CPU=3 #change accordingly in SLURM header
LG=LG3 #work on a subset of LG
GROUP=geno_LG3_homozyg #the subgroup on whcih we are making the fst comparison -> it should be a file like GROUP.txt in the folder 02_info
POP_FILE1=02_info/"$GROUP".txt #choose on which list of pop run the analyses

#do not edit if 12A has been run for each LG under the same filters. if not, please run before 
#because we are here re-using the beagle produced with all samples, on the given LG

MIN_MAF=0.10 #filter : will keep SNP above this allele frequency (over all individuals)
PERCENT_IND=0.75 #filter : will keep SNP with at least one read for this percentage of individuals 
MAX_DEPTH_FACTOR=3 #filter : will keep SNP with less than X time the nb of individuals total coverage

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
gunzip 12_ngsLD/$LG/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$LG".beagle.gz
#get the position file
gunzip -c 12_ngsLD/$LG/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$LG".mafs.gz | cut -f -2 | awk '{if(NR>1)print}' > 12_ngsLD/$LG/all_sites_maf"$MIN_MAF"_pctind_"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$LG" 
POS_FILE=12_ngsLD/$LG/all_sites_maf"$MIN_MAF"_pctind_"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$LG" 
#calculate the number of sites
N_SITES=$(wc -l $POS_FILE | cut -d " " -f 1)
	


#3 loop on populations to make sub-beagle, calculate LD & plot it for each subgroup
cat $POP_FILE1 | while read i
do
echo $i

N_IND=$(wc -l 12_ngsLD/"$GROUP"/"$i"subsetbam.filelist | cut -d " " -f 1)

# we are not running anymore angsd again because this was calling a different set of snp and diff MAj/min alleles
#but this would still be possible using the subsetbam.filelist
#just trying to subset the global beagle
	
	cat 12_ngsLD/"$GROUP"/"$i"_column.pos | while read k 
	do 
	echo $k
	cat 12_ngsLD/$LG/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$LG".beagle | cut -f $k > 12_ngsLD/$GROUP/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$LG".beagle
	gzip 12_ngsLD/$GROUP/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$LG".beagle
	done
	
	
	
	echo "calculating LD on " $i $LG
	#run ngsLD
	ngsLD --geno 12_ngsLD/$GROUP/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$LG".beagle.gz --probs --n_ind $N_IND --n_sites $N_SITES --pos $POS_FILE \
	--n_threads $NB_CPU --max_kb_dist 0 --max_snp_dist 0 --out 12_ngsLD/$GROUP/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$LG".ld --extend_out
	
	#smooth by window of X kb
	echo "smooth by window "
	cat 12_ngsLD/header 12_ngsLD/$GROUP/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$LG".ld > 12_ngsLD/$GROUP/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$LG"_header.ld
	
	python3 01_scripts/utility_scripts/ld_by_blocks.py \
	12_ngsLD/$GROUP/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$LG"_header.ld \
	1000 \
	12_ngsLD/$GROUP/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$LG"_header_by_1000.ld
	
	python3 01_scripts/utility_scripts/ld_by_blocks.py \
        12_ngsLD/$GROUP/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$LG"_header.ld \
        500 \
        12_ngsLD/$GROUP/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$LG"_header_by_500.ld

	python3 01_scripts/utility_scripts/ld_by_blocks.py \
        12_ngsLD/$GROUP/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$LG"_header.ld \
        250 \
        12_ngsLD/$GROUP/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$LG"_header_by_250.ld

	#make space
	rm 12_ngsLD/$GROUP/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$LG"_header.ld
	gzip 12_ngsLD/$GROUP/"$i"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$LG".ld
	
	
done
