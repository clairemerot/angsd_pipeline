#this R script pextract the pruned SNPs to make a list of sites fur subsequent analysis in angsd

argv <- commandArgs(T)
INPUT_plink <- argv[1]
INPUT_angsd <- argv[2]

library(dplyr)

pruned<-read.table(INPUT_plink)
head(pruned)
colnames(pruned)<-"LG_pos"
sites<-read.table(INPUT_angsd)
head(sites)
sites$LG_pos<-paste0(sites[,1],"_", sites[,2])


sites_pruned<-inner_join(sites, pruned)
print(paste("there is a total of ", dim(sites)[1], "sites"))
print(paste("we keep ", dim(sites_pruned)[1], "pruned sites"))

write.table(sites_pruned[,1:4], paste0(INPUT_angsd,"_pruned"), col.names=F, row.names=F, quote=F, sep="\t")