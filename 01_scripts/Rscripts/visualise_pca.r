#this R script performs a first visualisation in a pdf file with PC1-2 and PC3-4, and colours corresponding to groups
library(dplyr)
library(ade4)
argv <- commandArgs(T)
MIN_MAF <- argv[1]
PERCENT_IND <- argv[2]


#read info file
x_info<-read.table("02_info/info.txt", header=T)
pca.mat<-read.table(paste0("04_pca/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"cov.pca"), header=T)

if (dim(x_info)[1]!=dim(pca.mat)[1]){print ("warning : not the same number of indivuals in the pca and the info file")}

#just in case order is different, use rownames of the pca (coming from bam.list) to order the info matrix
bam_name<-as.data.frame(rownames(pca.mat))
colnames(bam_name)<-colnames(x_info)[1]
bam_x_info<-left_join(bam_name,x_info)
head(bam_x_info)

#for each var, make a pdf with PC1-2 & PC3-4 coloured 
nvar<-dim(x_info)[2]-2 #we assume that col1 is bamfile name & col 2 is ind id

for (i in 1:nvar)
{
#plot pca
group<-as.factor(bam_x_info[,i+2])
group_name<-colnames(bam_x_info)[i+2]
pdf(file=paste0("04_pca/all_maf",MIN_MAF,"_pctind",PERCENT_IND,".pca.",group_name,".pdf"))
par(mfrow=c(1,2))
s.class(pca.mat, group,  xax=1, yax=2, cellipse=T,grid=F, col=as.numeric(group), cstar=F,sub="PC1-2", pch=20)
s.class(pca.mat, group,  xax=3, yax=4, cellipse=T,grid=F, col=as.numeric(group), cstar=F,sub="PC3-4", pch=20)
dev.off()
}

