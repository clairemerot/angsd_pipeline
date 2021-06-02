#!/bin/bash
#SBATCH -J "12B_run_ngsLD_2LG"
#SBATCH -o log_%j
#SBATCH -c 6
#SBATCH -p large
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=21-00:00
#SBATCH --mem=30G

###this script calculate LD within and between 2 chromosomes
#It follows the prep step script 12A
###WARNING: before running ngsLD check the nb of sites (wc -l 12_ngsLD/sites)
#avoid running if this is > 200 000. 100 000 sites will already make super big files.
# ideally below 100 000

#maybe edit
NB_CPU=6 #change accordingly in SLURM header

LG_FILE1=02_info/chrs_fusion.txt #work on a chosen subset of LG 
#this is a text file with the name of the chromosomes. 1 chr by line. 
#Try to do a maximum of 2 chromosomes as this scrip will calculate intra-chr and inter-chr LD

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#module
module load angsd
ulimit -S -n 2048


#make a beagle and filter for maf an min ind
#this step was done by the script 12A

##input position file
POS_FILE=12_ngsLD/sites
#calculate the number of sites
N_SITES=$(wc -l $POS_FILE | cut -d " " -f 1)
N_IND=$(wc -l 02_info/bam.filelist | cut -d " " -f 1)


echo "calculating LD with $N_SITES SNPs "
#run ngsLD
ngsLD --geno 12_ngsLD/all.beagle.gz --probs --n_ind $N_IND --n_sites $N_SITES --pos $POS_FILE \
--n_threads $NB_CPU --max_kb_dist 0 --max_snp_dist 0 --out 12_ngsLD/all.ld --extend_out

gzip 12_ngsLD/all.ld