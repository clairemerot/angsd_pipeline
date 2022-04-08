#this R script output a file with the column to keep in the beagle file.


argv <- commandArgs(T)
GROUP <- argv[1] #file with the list of bamfile from the subgroup 
BAM_ALL <- argv[2] # file with all the bamfiles that was used to construct the whole beagle

library(dplyr)

BAM_beagle<-as.data.frame(read.table(BAM_ALL))
head(BAM_beagle)
colnames(BAM_beagle)<-c("file_name")

BAM_group<-as.data.frame(read.table(paste0("02_info/",GROUP,"bam.filelist")))
colnames(BAM_group)<-c("file_name")
BAM_group$keep<-"yes"


BAM_beagle_group<-left_join(BAM_beagle, BAM_group)

#this is the index of our sample of interest
pos_vec<-which(BAM_beagle_group$keep=="yes")
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
write.table(t(as.matrix(pos_beagle)), paste0("03_saf_maf_gl_all/subset_beagle/",GROUP,"_column.pos"),row.names=F, col.names=F, quote=F, sep=",")





