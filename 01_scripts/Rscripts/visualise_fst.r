#this R script performs a first visualisation in a pdf file of pairwise
library(corrplot)

argv <- commandArgs(T)
MIN_MAF <- argv[1]
PERCENT_IND <- argv[2]
POP<-argv[3]


#read pop file
print(POP)
pop<-read.table(POP, header=F)
npop<-dim(pop)[1]
pop_group<-unlist(strsplit(unlist(strsplit(POP,"/"))[2],".txt"))

FST_pop<-matrix(nrow=npop, ncol=npop)
rownames(FST_pop)=colnames(FST_pop)=pop[,1]

for (j in 1:(npop-1))
{
  for (i in (j+1):npop)
  {
    pi<-pop[i,1]
	pj<-pop[j,1]
	FSTi<-read.delim(paste0("07_fst_by_pop_pair/",pj, "_", pi,"_maf",MIN_MAF,"_pctind",PERCENT_IND,".fst",sep=""), header=F, row.names = NULL, sep=" ")
    FST_pop[j,j]<-0
    FST_pop[i,i]<-0
    FST_pop[i,j]<-FSTi[1,2]
    FST_pop[j,i]<-FSTi[1,2]
	}
}
print(FST_pop)
write.table (FST_pop, paste0("07_fst_by_pop_pair/all",pop_group, "_maf",MIN_MAF,"_pctind",PERCENT_IND,".fstmatrix"), quote=F)
	
#visualise the fst matrix
pdf(file=paste0("07_fst_by_pop_pair/all",pop_group, "_maf",MIN_MAF,"_pctind",PERCENT_IND,".fstmatrix.pdf"))
par(mfrow=c(1,1))
corrplot(FST_pop, is.corr = F, method="number", diag=T, type="lower", number.digits = 3)
dev.off()







