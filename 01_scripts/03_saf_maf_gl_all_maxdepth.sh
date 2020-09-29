#!/bin/bash
#SBATCH -J "03_saf_maf_gl_all_maxdepth"
#SBATCH -o log_%j
#SBATCH -c 4 
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=7-00:00
#SBATCH --mem=100G

###this script will work on all bamfiles and calculate saf, maf & genotype likelihood
#maybe edit
NB_CPU=4 #change accordingly in SLURM header
#REGIONS="-rf 02_info/regions_25kb_100snp.txt" #optional edit with your region selected file
REGIONS="" # to remove the options to focus on a limited number of regions

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
angsd -P $NB_CPU -nQueueSize 50 \
-doMaf 1 -dosaf 1 -GL 2 -doGlf 2 -doMajorMinor 1 -doCounts 1 \
-anc 02_info/genome.fasta -remove_bads 1 -minMapQ 30 -minQ 20 \
-minInd $MIN_IND -minMaf $MIN_MAF -setMaxDepth $MAX_DEPTH \
-b 02_info/bam.filelist \
$REGIONS -out 03_saf_maf_gl_all/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"

#main features
#-P nb of threads -nQueueSize maximum waiting in memory (necesary to optimize CPU usage
# -doMaf 1 (allele frequencies)  -dosaf (prior for SFS) -GL (Genotype likelihood 2 GATK method - export GL in beagle format  -doGLF2) 
# -doMajorMinor 1 use the most frequent allele as major
# -anc provide a ancestral sequence = reference in our case -fold 1 (car on utilise la ref comme ancestral
# -rf (file with the region written) work on a defined region : OPTIONAL
# -b (bamlist) input file
# -out  output file

#main filters
#filter on bam files -remove_bads (remove files with flag above 255) -minMapQ minimum mapquality -minQ (minimum quality of reads?)
#filter on frequency -minInd (minimum number of individuals with at least one read at this locus) we set it to 50%
#filter on allele frequency -minMaf, set to 0.05 

#extract SNP which passed the MIN_MAF and PERCENT_IND filters & their Major-minor alleles
#order the sites file by chromosome names 
#makes a region file matching the sites files and with same order
#index sites file
echo "from the maf file, extract a list of SNP chr, positoin, major all, minor all"
gunzip 03_saf_maf_gl_all/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".mafs.gz 
Rscript 01_scripts/Rscripts/make_sites_list_maxdepth.R "$MIN_MAF" "$PERCENT_IND" "$MAX_DEPTH_FACTOR"

angsd sites index 02_info/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"

#Eric propose a much shorter version using bash to cut the 4 columns of the mafs.gz. but then angsd is unable to index it
#I can't find the problem, so I came back to older solution with R which works
#gunzip -c 03_saf_maf_gl_all/all_maf"$MIN_MAF"_pctind_"$PERCENT_IND".mafs.gz | cut -f -4 > 02_info/sites_all_maf"$MIN_MAF"_pctind_"$PERCENT_IND"



