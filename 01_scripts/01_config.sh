#path to bam folder
BAM_PATH=../wgs_sample_preparation/09_no_overlap

#path to pcaangsd
PCA_ANGSD_PATH=/project/lbernatchez/programs/pcangsd 

#filter : will keep SNP above this allele frequency (over all individuals)
MIN_MAF=0.05 

#filter : will keep SNP with at least one read for this percentage of individuals (over all individuals in step 03, and within each pop at step 07)
PERCENT_IND=0.5 

#window size for sliding window FST & Thetas
WINDOW=25000 

#window step
WINDOW_STEP=5000 

#min nb of pop to consider for NGS admix
K_MIN=2 

#maximum nb of pop to consider for NGS admix
K_MAX=5 
