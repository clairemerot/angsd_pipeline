#this R script uses the mafs file provided by the analyses on all individuals, with filter for quality and maf 0.05
#it simply extract the first columns with chr, position of each SNP and Major/minor alleles as determined in step 03
#output is a "sites" files that will allwo restraining subsequent analyses 
#(for instance maf by pop, FSt by pop pairs, etc to a limited number of loci

argv <- commandArgs(T)
MIN_MAF <- argv[1]
PERCENT_IND <- argv[2]
maf<-read.table( paste0("all_maf",MIN_MAF,"_pctind",PERCENT_IND,".mafs"), header=T)
head (maf)
write.table(maf[,1:4], paste0("sites_all_maf", MIN_MAF, "_pctind", PERCENT_IND), row.names=F, col.names=F, quote=F)