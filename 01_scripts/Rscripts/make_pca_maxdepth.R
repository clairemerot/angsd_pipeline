#this R script performs a PCa on the covariance matrix and export it labelling rows with bamfile name and columns as PC
#it also makes a first visualisation in a pdf file with PC1-2 and PC3-4

argv <- commandArgs(T)
MIN_MAF <- argv[1]
PERCENT_IND <- argv[2]
BAM <- argv[3]
MAX_DEPTH_FACTOR <-argv[4]

#perform a cpa on covariance matrix
cov_mat<-as.matrix(read.table(paste0("04_pca/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,".cov"), header=F))
pca<-prcomp(cov_mat)

#add column names
nPC<-dim(pca$x)[2]
col_PC<-vector(length=nPC)
for (i in 1 : nPC) {col_PC[i]<-paste0("PC",i)}
colnames(pca$x)<-c(col_PC)

#add rownames
bam_names<-read.table(BAM,header=F)
rownames(pca$x)<-bam_names$V1

#plot pca
pdf(file=paste0("04_pca/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,".pca.pdf"))
par(mfrow=c(1,2))
plot(pca$x[,1], pca$x[,2], pch=20, ylab="PC2", xlab="PC1")
plot(pca$x[,3], pca$x[,4], pch=20, ylab="PC4", xlab="PC3")
dev.off()

write.table(pca$x, paste0("04_pca/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,".pca"), quote=F)
