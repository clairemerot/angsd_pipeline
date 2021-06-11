#!/bin/bash
#SBATCH -J "12B_ngsLD_byLG"
#SBATCH -o log_%j
#SBATCH -c 4
#SBATCH -p large
#SBATCH --mail-type=ALL
#SBATCH --mail-user=X
#SBATCH --time=21-00:00
#SBATCH --mem=100G

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

###this script will use a beagle by chromosome, calculate LD by chromosome and output a reduce matrix
#IMPORTANT: you should check how many sites you have per chromosome.
#if this is >100000 sites DO NOT RUN: this will be very very big matrices that may crash the server or fill it quicker than you think (LD matrix is big as TB!!) 
#this can be easily checked for each chromosome by doing
# wc -l 12_ngsLD/$i/sites_maf"$MIN_MAF"_pctind_"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$i" where $i is the chromosome name

#new option that can be used if you have too many sites [warning: not run yet]
RANDOM_FRACTION=0.10 #randomly subsampling this fraction of sites to calculate LD


NB_CPU=4 #change accordingly in SLURM header
LG_FILE1=02_info/region_LG.txt #work on a chosen subset of LG
#same filters as used for preparation
MIN_MAF=0.10 #filter : will keep SNP above this allele frequency (over all individuals)
PERCENT_IND=0.75 #filter : will keep SNP with at least one read for this percentage of individuals 
MAX_DEPTH_FACTOR=3 #filter : will keep SNP with less than X time the nb of individuals total coverage




N_IND=$(wc -l 02_info/bam.filelist | cut -d " " -f 1)


#loop over LG to run ngsLD

cat $LG_FILE1 | while read i
do
echo $i

BEAGLE_FILE=12_ngsLD/$i/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$i".beagle.gz
POS_FILE=12_ngsLD/$i/sites_maf"$MIN_MAF"_pctind_"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$i"
LDOUT_FILE=12_ngsLD/$i/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$i".ld
#calculate the number of sites
N_SITES=$(wc -l $POS_FILE | cut -d " " -f 1)

echo "calculating LD on " $i

#run ngsLD
ngsLD --geno $BEAGLE_FILE --probs --n_ind $N_IND --n_sites $N_SITES --pos $POS_FILE \
--n_threads $NB_CPU --max_kb_dist 0 --max_snp_dist 0 --rnd_sample $RANDOM_FRACTION --out $LDOUT_FILE --extend_out


#run Eric's optimized script on compressed file
gzip $LDOUT_FILE

cat 12_ngsLD/header.gz 12_ngsLD/$i/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$i".ld.gz > 12_ngsLD/$i/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$i"_header.ld.gz


##run Eric's script to keep quartile values per windows of 1000kb, 500 kb, 250 kb
python3 01_scripts/utility_scripts/ld_by_blocks_optimized_gzinput.py 12_ngsLD/$i/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$i"_header.ld.gz 1000 12_ngsLD/$i/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$i"_header_by_1000.ld
python3 01_scripts/utility_scripts/ld_by_blocks_optimized_gzinput.py 12_ngsLD/$i/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$i"_header.ld.gz 500 12_ngsLD/$i/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$i"_header_by_500.ld
python3 01_scripts/utility_scripts/ld_by_blocks_optimized_gzinput.py 12_ngsLD/$i/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$i"_header.ld.gz 250 12_ngsLD/$i/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$i"_header_by_250.ld

#make space : remove the duplicate with the header
rm 12_ngsLD/$i/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_"$i"_header.ld.gz

done
