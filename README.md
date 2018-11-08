# angsd_pipeline

#run all commands from the angsd_pipeline folder

## 00_DEPENDANCIES
install angsd & associated programs
http://www.popgen.dk/angsd/index.php/ANGSD

add angsd to the path in .bashrc

add the misc folder (containing RealSFS, theta stat etc to the path in .bashrc

export PATH="/home/camer78/Softwares/angsd2/angsd:$PATH"
export PATH="/home/camer78/Softwares/angsd2/angsd/misc:$PATH"

install NGSAdmix (maybe in the misc folder, else export its path)
http://www.popgen.dk/software/index.php/NgsAdmix

install pcangsd (maybe in the misc folder) & check if you have python2
http://www.popgen.dk/software/index.php/PCAngsd
copy the path into 01_config.sh PCA_ANGSD_PATH=~/Softwares/pcangsd

for all script file, you may edit the header to put your email adress and adjust cpu/memory/time/allocation and slurm partition 

## 01_PREPARE_DATA

input: 
- bam-files
They must be aligned to the reference, indexed and sorted, named like "id_sex_pop_group_blablabla.sorted.bam.
They can  be kept in original folder, just know the path to their location from the angsd_pipeline folder

insert this path in BAM_PATH= in the 01_config.sh

- info.txt in 02_info folder: a file listing the bamfile names ordered  with a column for any relevant information on the individuals
for follow-up analyses with R ideally: col1=bam_filename, col2=id, col3=sex, col4=pop, col5=group, col6=group ...
- pop.txt in 02_info folder : a file listing population names with one item by line (there can be several files if we aimed at analysing different grouping, group.txt)
- genome.fasta in 02_info folder: the reference genome on which bam have been aligned

if it is not indexed run:

module load samtools

samtools faidx 02-info/genome.fasta

- region.txt: a file listing the regions of the genome (chromosome or scaffolds) to be included in the analysis.
NOTE THIS POSSIBILITY HAS BEEN TEMPORALLY REMOVED

For initial tests, put manually just a few, one per line

To list all scaffolds of the genome

grep -e ">" 02_info/genome.fasta | awk 'sub(/^>/, "")' | sort -k1 > 02_info/region.txt

- edit the 01_config.sh file

choose MIN_MAF, PERCENT_IND filters, WINDOW and WINDOW_STEP for sliding-windows analyses, and K_MIN, K_MAX for admixture analysis

## 02_LIST_BAMFILES_AND_LIST_BY_POP
this script will make a list of all bamfiles and several list by population, based on information in bamfile names and pop.txt

edit script if you want to add another way of grouping (group.txt)

./01_scripts/02_list_bamfiles.sh

## 03_RUN_INITIAL_ANALYSIS_ON_WHOLE_DATASET
this script will work on all bamfiles and calculate saf, maf & genotype likelihood on the whole dataset. It will output in 
02_info folder the list of SNP which passed the MIN_MAF and PERCENT_IND filters & their Major-minor alleles (sites_*)

maybe edit the number of cpu NB_CPU and allocated memory/time

sbatch 01_scripts/03_saf_maf_gl_all.sh


## 04_PCA_VISUALISE_CHECK_WHETHER_YOU_WANT_TO_EXCLUDE_OUTLIERS
this script will work on all individuals using the beagle genotype likelihood and calculate a covariance matrix with angsd.
it also output the pca with R, and visualisation in pdf

check for outliers (or duplicates), one may to re-run step 03 and 04 with an edited bamlist

this requires pcangsd to be cloned and a version of Python v2 with alias python2

maybe edit NB_CPU and memory (sometimes require a lot of memory >100 G)

sbatch 01_scripts/04_pca.sh

#for further visualisation using information from info.txt

source 01_scripts/01_config.sh

Rscript 01_scripts/Rscripts/visualise_pca.r "$MIN_MAF" "$PERCENT_IND"

## 05_ADMIXTURE_ANALYSIS
this script will work on all individuals using the beagle genotype likelihood and perform an admixture analysis. 
this requires NGSadmix to be installed and its path export in the bashrc. 
NGS admiw will explore all number of population between K_MIN and K_MAX as provided in the 01_config.sh.

maybe edit NB_CPU=1 & edit K_MIN and K_MAX in the 01_config.sh

sbatch 01_scripts/05_ngs_admix.sh

#for further visualisation using information from info.txt

source 01_scripts/01_config.sh

BAM_LIST=02_info/bam.filelist

Rscript 01_scripts/Rscripts/visualise_admix.r "$MIN_MAF" "$PERCENT_IND" "$K_MIN" "$K_MAX" "$BAM_LIST"

## 06_CALCULATE_ALLELES_FREQUENCIES_BY_POP
this script will work on bamfiles by population and calculate saf  & maf. 
It will run on the list of sites determined at step 3 (filter on global population). 
In addition it will filter for sites with at least one read in a minimum proportion of individuals within each pop

maybe edit cpu & choose on which list of pop run the analyses

NB_CPU=1 & POP_FILE1=02_info/pop.txt 

sbatch 01_scripts/06_saf_maf_by_pop.sh

## 07_CALCULATE_PAIRWISE_FST
This script will use the saf by population calculated at step 07 and calculate SFS and FST

maybe edit 

NB_CPU=1 #change accordingly in SLURM header

POP_FILE1=02_info/pop.txt #choose on which list of pop run the analyses

sbatch 01_scripts/07_fst_by_pop_pair.sh

#for further visualisation (requires the corrplot package)

install.packages("corrplot")

POP1_FILE=02_info/pop.txt

Rscript 01_scripts/Rscripts/visualise_fst.r "$MIN_MAF" "$PERCENT_IND" "$POP1_FILE"

## 08_CALCULATE_THETAS
this script will use the saf on all individuals and the saf by population calculated at step 07
to calculate thetas statistics.Beware if ancestral sequenc is the reference (folded spectrum), not all stats are meaningful.

sbatch 01_scripts/08_thetas.sh

## 09_MAKE_GWAS
those two scripts can do a gwas either with binary phenoty "_bin" or quantitative phenotype "_quant".

phenotype files should be put in 02_info
- bin_pheno.txt #this file must be one single column with phenotype coded as 1 or 2, each line is one individual in the same order as bamfile
- quant_pheno.txt #this file must be one single column with quantitative phenotype, each line is one individual in the same order as bamfiles. 
Beware, the quantitative gwas is made to include a covariable (for instance sex) coded as binary # change the $COV if needed

sbatch 01_scripts/09_gwas_bin.sh

sbatch 01_scripts/09_gwas_quant.sh

## 10_ANALYSING_MAF_SELECTION_TESTS_ETC
See selection_pipeline

## 11_ANALYSING_LD
possibility to output plink format from angsd and then further analyses in plink?
