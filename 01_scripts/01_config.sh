#path to bam folder
BAM_PATH=../wgs_sample_preparation/09_no_overlap

#path to pcaangsd
#this is the former version of ANGSD - it uses python2 and may require some libraries.
PCA_ANGSD_PATH=/project/lbernatchez/programs/pcangsd 

#filter : will keep SNP above this allele frequency (over all individuals)
MIN_MAF=0.05 

#filter : will keep positions with at least MIN_DEPTH reads for each individual, This is not necessaily for all individuals, we consider a PERCENT_IND (percentage of individuals over all individuals in step 03, and within each pop at step 07)
#advice: as min depth use a value that is a bit below what you expected. we use 1 for 1X of coverage but if you have 5X it may make sense to put the bar a bit higher to 2 or 3
#advice: as percentage, avoid going below 50% and also consider the whole number of individuals. (it may makes sense to use 50% with 100 ind/pop, but you may want 90% with 9 ind/pop
PERCENT_IND=0.5
MIN_DEPTH=1

#filter: will keep SNP with at least a coverage of this factor multiplied by the number of ind - across all ind. 
# advice: we usually set it at 2-4 times the expected coverage to remove repeated regions
MAX_DEPTH_FACTOR=3

#window size for sliding window FST & Thetas
WINDOW=25000 

#window step
WINDOW_STEP=5000 

#min nb of pop to consider for NGS admix
K_MIN=2 

#maximum nb of pop to consider for NGS admix
K_MAX=5 
