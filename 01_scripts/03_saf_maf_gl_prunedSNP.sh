#!/bin/bash
#SBATCH -J "03_saf_maf_gl_prunedSNP"
#SBATCH -o log_%j
#SBATCH -c 4 
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=XXX
#SBATCH --time=7-00:00
#SBATCH --mem=150G

###this script will work on all bamfiles and calculate saf, maf & genotype likelihood
#maybe edit
NB_CPU=4 #change accordingly in SLURM header
#REGIONS="-rf 02_info/regions_25kb_100snp.txt" #optional edit with your region selected file
REGIONS=" " # to remove the options to focus on a limited number of regions

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

module load angsd
ulimit -S -n 2048

#prepare variables - avoid to modify
source 01_scripts/01_config.sh
N_IND=$(wc -l 02_info/bam.filelist | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

echo " Calculate the SAF, MAF and GL for all individuals listed in 02_info/bam.filelist"
echo "keep loci with at leat one read for n individuals = $MIN_IND, which is $PERCENT_IND % of total $N_IND individuals"
echo "filter on allele frequency = $MIN_MAF"

####Calculate the SAF, MAF and GL
#we now use a list of snp corresponding to the LD-pruned snp by plink and removing the areas of inversions
angsd -P $NB_CPU -nQueueSize 50 \
-dosaf 1 -GL 2 -doGlf 2 -doMajorMinor 3 \
-anc 02_info/genome.fasta -remove_bads 1 -minMapQ 30 -minQ 20 \
-sites 02_info/sites_"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_LDpruned \
-rf 02_info/regions_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR" \
-b 02_info/bam.filelist \
-out 03_saf_maf_gl_all/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_LDpruned

#main features
#-P nb of threads -nQueueSize maximum waiting in memory (necesary to optimize CPU usage
# -doMaf 1 (allele frequencies)  -dosaf (prior for SFS) -GL (Genotype likelihood 2 GATK method - export GL in beagle format  -doGLF2) 
# -doMajorMinor 1 use the most frequent allele as major
# -anc provide a ancestral sequence = reference in our case -fold 1 (car on utilise la ref comme ancestral
# -rf (file with the region written) work on a defined region : OPTIONAL
# -b (bamlist) input file
# -out  output file

