#this R script output bamfilelist by downsampling the different populations to compare population with the same sample size


argv <- commandArgs(T)
GROUP <- argv[1]

pop<-read.table(paste0("02_info/",GROUP,".txt"))
n_pop<-dim(pop)[1]

dim_pop<-vector(length=n_pop)
list_bamlist<-list()

for (i in 1 : n_pop)
{
list_bamlist[[i]]<-read.table(paste0("02_info/",pop[i,1],"bam.filelist"))
dim_pop[i]<-dim(list_bamlist[[i]])[1]
}

n_ind<-min(dim_pop)

for (i in 1 : n_pop)
{
X_full<-list_bamlist[[i]]
X_rand<-X_full[sample(c(1:dim(X_full)[1]), n_ind),1]
write.table(X_rand, paste0("07_fst_by_pop_pair/",GROUP,"/", pop[i,1], "subsetbam.filelist"),row.names=F, col.names=F, quote=F)
}



