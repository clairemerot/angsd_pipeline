#this R script output bamfilelist by downsampling the different populations to compare population with the same sample size
#it now also output a file with the column to keep in the beagle file.


argv <- commandArgs(T)
GROUP <- argv[1] #list of subgroup on which we will loop
BAM_ALL <- argv[2] # file with all the bamfiles that was used to construct the whole beagle

library(dplyr)

BAM_beagle<-as.data.frame(read.table(BAM_ALL))
head(BAM_beagle)
colnames(BAM_beagle)<-c("file_name")

pop<-read.table(paste0("02_info/",GROUP,".txt"))
n_pop<-dim(pop)[1]

dim_pop<-vector(length=n_pop)
list_bamlist<-list()

#put in a list the bamfile for eahc subgroups
for (i in 1 : n_pop)
{
list_bamlist[[i]]<-read.table(paste0("02_info/",pop[i,1],"bam.filelist"))
dim_pop[i]<-dim(list_bamlist[[i]])[1]
print(dim_pop[i])
}

#calculate the size of the smaller group
n_ind<-min(dim_pop)

for (i in 1 : n_pop)
{
print(pop[i,1])
#for each subgroup take a random subset of individuals
X_full<-list_bamlist[[i]]
random_list<-sample(c(1:dim(X_full)[1]), n_ind,replace = FALSE)
X_rand<-X_full[random_list,1]

#write the random bamfilelist
write.table(X_rand, paste0("12_ngsLD/",GROUP,"/", pop[i,1], "subsetbam.filelist"),row.names=F, col.names=F, quote=F)

#merge random bam list with global bamlist to know the position of targeted sample in the big beagle file
X_rand<-as.data.frame(X_rand)
print(head(X_rand))
colnames(X_rand)<-c("file_name")
X_rand$keep<-rep("yes", dim(X_rand)[1])

BAM_beagle_X_rand<-left_join(BAM_beagle, X_rand)

#this is the index of our sample of interest
pos_vec<-which(BAM_beagle_X_rand$keep=="yes")
print(pos_vec)

#this will be the column of our sample of interest in the bamfile
#initialisation: 3 column with position, major minor
pos_beagle<-c(1,2,3)

for (j in 1 : length(pos_vec))
{
k<-pos_vec[j]
pos_beagle<-c(pos_beagle, (3*k)+1, (3*k)+2, (3*k)+3)
}

print(t(as.matrix(pos_beagle))[1:20])
#write the position of the target samples in a file to split the beagle
write.table(t(as.matrix(pos_beagle)), paste0("12_ngsLD/",GROUP,"/", pop[i,1], "_column.pos"),row.names=F, col.names=F, quote=F, sep=",")

}



