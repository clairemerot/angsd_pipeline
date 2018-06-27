#this R script performs a first visualisation in a pdf file with PC1-2 and PC3-4, and colours corresponding to groups
library(dplyr)


argv <- commandArgs(T)
MIN_MAF <- argv[1]
PERCENT_IND <- argv[2]
K_MIN<-argv[3]
K_MAX<-argv[4]
BAM <- argv[5]




#read info file
x_info<-read.table("02_info/info.txt", header=T)
bam_name<-read.table(BAM,header=F)
colnames(bam_name)<-colnames(x_info)[1]
if (dim(x_info)[1]!=dim(bam_name)[1]){print ("warning : not the same number of indivuals in the pca and the info file")}
#just in case order or nb is different, use bam.list to order the info matrix
bam_x_info<-left_join(bam_name,x_info)
head(bam_x_info)

nvar<-dim(x_info)[2]-2 #we assume that col1 is bamfile name & col 2 is ind id

admix_info<-bam_x_info

print(K_MIN)
print(K_MAX)

for (k in K_MIN : K_MAX)
	{
	print(paste0("k=",k))
	#read admixture file for k
	admix_value<-read.table(paste0("05_ngs_admix/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_K",k,".qopt"), header=F)
	admix<-t(as.matrix(admix_value))
	colnames(admix_value)<-rep(paste0("k=",k),k)
	
	#keep admixture Q value with info
	admix_info<-cbind(admix_info, admix_value)

	#open a graphic
	pdf(file=paste0("05_ngs_admix/all_maf",MIN_MAF,"_pctind",PERCENT_IND,".admix.k",k,".pdf"))
	par(mfrow=c((nvar+1),1))

	#structure plot unordered
	barplot(admix,col=1:k,space=0,border=NA,xlab="Individuals",ylab="admixture", main="unordered")
	
	#for each var, make a pdf with barplot ordered follwing this group
	for (i in 1:nvar)
		{
		group<-as.factor(bam_x_info[,i+2])
		group_name<-colnames(bam_x_info)[i+2]
		admixi<-admix[,order(group)]
		barplot(admixi,col=1:k,space=0,border=NA,xlab="Individuals",ylab="admixture", main=paste0("ordered by ",group_name))
		}
	dev.off()
}
write.table(admix_info, paste0("05_ngs_admix/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"all_K.admix"), quote=F, row.names=F)
