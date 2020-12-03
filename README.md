# angsd_pipeline

#run all commands from the angsd_pipeline folder
#for an overview of the pipeline please see
https://drive.google.com/file/d/14bmwOkdbdfSsfNDrYNxR2V8kHZyNuIlm/view?usp=sharing

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

Useful for most analysis:

- info.txt in 02_info folder: a file listing the bamfile names ordered  with a column for any relevant information on the individuals
for follow-up analyses with R ideally: col1=bam_filename, col2=id, col3=sex, col4=pop, col5=group, col6=group ...
- pop.txt in 02_info folder : a file listing population names with one item by line (there can be several files if we aimed at analysing different grouping, pop1.txt, pop_geo.txt, etc)

Necessary:

- genome.fasta in 02_info folder: the reference genome on which bam have been aligned

if it is not indexed run:
```
module load samtools

samtools faidx 02-info/genome.fasta
```
- region.txt: a file listing the regions of the genome (chromosome or scaffolds) to be included in the analysis.

For initial tests, put manually just a few, one per line

To list all scaffolds of the genome
You may want to reduce this list to improve computation time and exclude short unanchored scaffold. 
```
grep -e ">" 02_info/genome.fasta | awk 'sub(/^>/, "")' | sort -k1 > 02_info/region.txt
```
- edit the 01_config.sh file

choose MIN_MAF, PERCENT_IND filters, MAX_DEPTH filter for intials steps
For later steps, choose WINDOW and WINDOW_STEP for sliding-windows analyses, and K_MIN, K_MAX for admixture analysis

## 02_LIST_BAMFILES_AND_LIST_BY_POP
this script will make a list of all bamfiles and several list by population, based on information in bamfile names and pop.txt

edit script if you want to add another way of grouping (group.txt)
```
./01_scripts/02_list_bamfiles.sh
```
## 03_RUN_INITIAL_ANALYSIS_ON_WHOLE_DATASET
this script will work on all bamfiles and calculate saf, maf & genotype likelihood on the whole dataset. It will output in 
02_info folder the list of SNP which passed the MIN_MAF and PERCENT_IND filters & their Major-minor alleles (sites_*)

maybe edit the number of cpu NB_CPU and allocated memory/time

```
sbatch 01_scripts/03_saf_maf_gl_all.sh
```
## 04_PCA_VISUALISE_CHECK_WHETHER_YOU_WANT_TO_EXCLUDE_OUTLIERS
this script will work on all individuals using the beagle genotype likelihood and calculate a covariance matrix with angsd.
it also output the pca with R, and visualisation in pdf

check for outliers (or duplicates), one may to re-run step 03 and 04 with an edited bamlist

this requires pcangsd to be cloned and a version of Python v2 with alias "python2"

maybe edit NB_CPU and memory (sometimes require a lot of memory >100 G)
```
sbatch 01_scripts/04_pca.sh
```
for further visualisation using information from info.txt, the script 01_scripts/Rscripts/visualise_pca.r can be useful


## 05_ADMIXTURE_ANALYSIS
this script will work on all individuals using the beagle genotype likelihood and perform an admixture analysis. 
this requires NGSadmix to be installed and its path export in the bashrc. 
NGS admiw will explore all number of population between K_MIN and K_MAX as provided in the 01_config.sh.

maybe edit NB_CPU=1 & edit K_MIN and K_MAX in the 01_config.sh

sbatch 01_scripts/05_ngs_admix.sh

for further visualisation using information from info.txt, the script 01_scripts/Rscripts/visualise_admix.r can be useful

## 06_CALCULATE_ALLELES_FREQUENCIES_BY_POP
this script will work on bamfiles by population and calculate saf  & maf. 
It will run on the list of sites determined at step 3 (filter on global population). Major and minor alleles are polarized by the list of SNPs from step 3 which means 
that an allele can be minor at the scale of all populations but at frequency >50% in a given population. Keeping this polarisation is important because we want to have the frequency of the same allele accross popualtions.
In addition it will filter for sites with at least one read in a minimum proportion of individuals within each pop

maybe edit cpu & choose on which list of pop run the analyses
NB_CPU=1 & POP_FILE1=02_info/pop.txt 
```
sbatch 01_scripts/06_saf_maf_by_pop.sh
```
The resulting MAF by population are the data used by the selection_analysis pipeline which does environmental associations.

## 07_CALCULATE_PAIRWISE_FST
This script will calculate the unfold saf by population; then the 2dSFS and FST for each pair of populations
It starts with a R script that subset the population bamlist ot have the same number of individuals as this factor can strongly influence Fst values.

maybe edit
```
NB_CPU=1 #change accordingly in SLURM header

POP_FILE1=02_info/pop.txt #choose on which list of pop run the analyses

sbatch 01_scripts/07_fst_by_group.sh
```
for further visualisation (requires the corrplot package), you may use 01_scripts/Rscripts/visualise_fst.r 

## 08_CALCULATE_THETAS
this script will NOT filter on Maf as we want to keep all positions to calculate thetas statistics. It will simply filter on coverage with the same parameters fixed in 01_config.sh. It calculates the saf, 1DSFS and thetas statistics by population. I have tried on all populations together but it does not really make sense and it is impossible to run on thousands of individuals.

Beware if ancestral sequenc is the reference (folded spectrum), not all stats are meaningful.
```
sbatch 01_scripts/08_thetas_by_pop.sh
```
## 09_MAKE_GWAS
those two scripts can do a gwas either with binary phenoty "_bin" or quantitative phenotype "_quant".

phenotype files should be put in 02_info (see some examples herein)
- bin_pheno.txt #this file must be one single column with phenotype coded as 1 or 2, each line is one individual in the same order as bamfile
- quant_pheno.txt #this file must be one single column with quantitative phenotype, each line is one individual in the same order as bamfiles. 
Beware, the quantitative gwas is made to include a covariable (for instance sex) coded as binary # change the $COV if needed
Missing data must be coded -999

In the quant script, There are three ways of implementing the GWAS. See Angsd help about it http://www.popgen.dk/angsd/index.php/Association
```
sbatch 01_scripts/09_gwas_bin.sh

sbatch 01_scripts/09_gwas_quant.sh
```
## 10 MAKING PCA BY WINDOW ALONG GENOME
This step include several modules that should be run successively
important: edit window size - this is a number of SNPs
windows od 100 to 10 000 SNPs allows analysis along the genome, while a large window can be chosen if one wants to make a pca by chromosome (see also 10C to directly do a local PCA on a given window)

A- This script will call a python script written by Eric Normandeau to split the begale files into windows of a given size within 
each chromosome/scaffold. They are stored into beagle_by_window folder
Then, it will run pcangsd on each window of X SNPs. Covariances matrices are stored into cov_by_window folder
```
sbatch 01_scripts/10A_cov_by_window.sh
```

B- This script will call a R script using the lostruct package https://github.com/petrelharp/local_pca
It applies the method proposed in Li, H., & Ralph, P. (2019). Local PCA shows how the effect of population structure differs 
along the genome. Genetics, 211(1), 289-304.

In addition, the R script test correlation between the PC1 scores of each window and the PCs of the global pca done at step 04

IMPORTANT: edit window size to fit what was given to the A script, edit the nb of individuals and the nombre of PC to consider in lostruct, and the nb of axis of the MDS to look at and save

```
window_size=100 #nb of Snps per window
N_IND=1446 #nb of individuals included in the analysis
N_PC=2 #nb of PC to consider when comparing windows
N_MDS=50 #nb of MDS dimension in the output
```

```
sbatch 01_scripts/10B_pca_lostruct.sh
```

C- After the analysis of the MDS which may highlight regions of interest, or if one is interested in doing PCA on a chosen window, this script will call select the region in the beagle using a python script written by Eric Normandeau, call Pcangsd to get the covariance matrix, and a R script to get the PCa and plot it.

important: this requires a file with two column with the start and stop of the wanted region separated by a tab
for instance
```
LG1_125 LG1_456
```

It will output pca, eigen-values and images.
```
sbatch 01_scripts/10C_pca_chosen_regions.sh
```

## 11 LD pruning with PLINK
This step allow outputing genotype likelihoods in plink format and run plink to extract a list of LD-pruned SNPs.

by default values for Plink pruning are as follow and can be adjusted directly in the script 
```
WINDOW=100
SNP=5
VIF=2

sbatch 01_scripts/11_plink_pruning.sh
```
## 12 LD calculation by NGSLD
This step allow outputing genotype likelihoods in plink format and run plink to extract a list of LD-pruned SNPs.
```
sbatch 01_scripts/11_plink_pruning.sh
```
## 13 Hardy-Weinberg statistics and Hobs
This step use the -hwe module of angsd to output HW statistics per population/group.

I have add a Rscript it to extract the number of heterozygotes observed at each SNP.

We later uses the windowscannr library in R to calculate Hobs along the genome for each group
```
sbatch 01_scripts/13_hwe.sh
```
## ANALYSING_MAF_SELECTION_TESTS_ETC
See selection_pipeline

