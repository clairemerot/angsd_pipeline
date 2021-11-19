#this R script takes as input the hwe info and output a matrix with average Hobs and average deviation from HW (F value) in sliding-windows.
#requires library windowscanr

#library(devtools)
#install_github('tavareshugo/windowscanr')

library(windowscanr)

argv <- commandArgs(T)
FILE <- argv[1]
WINDOW<- as.numeric(argv[2])
WINDOW_STEP<- as.numeric(argv[3])

print(paste("loading",FILE))
hwe<-read.table(FILE, header=T)

#calculate Hexp & Hobs
#the expected proportion of heterozygotes can be derived from the allelic frequency under HW (Hexp=2pq)
hwe$Hexp<-2*(hwe$hweFreq)*(1-hwe$hweFreq) #Hexp = 2*(hweFreq)(1-hweFreq)
#the observed proportion of heterozygotes can now be calculated using F, the departur from HW (F=1-Hobs/Hexp)
hwe$Hobs<-hwe$Hexp-(hwe$F*hwe$Hexp) #Hobs= Hexp - F* Hexp
head(hwe)

#get the propotion of homozygotes
hwe$Hminor<-hwe$Freq*hwe$Freq+hwe$Freq*(1-hwe$Freq)*hwe$F
#and the proportion of homozygotes at the major allele
hwe$Hmajor<-(1-hwe$Freq)*(1-hwe$Freq)+hwe$Freq*(1-hwe$Freq)*hwe$F

#output the same hwe matrix results but with estimated Hobs - big file - value are for each position
write.table(hwe, paste0(FILE,".Hobs"), row.names=F, quote=F, sep="\t")

#do sliding-window (time consuming)
print(paste("compute sliding window of size",WINDOW, "with step",WINDOW_STEP))
hwe_win <- winScan(x = hwe,groups = "Chromo", position = "Position",values = c("Hobs", "F"),win_size = WINDOW,win_step = WINDOW_STEP,funs = c("mean"))
head(hwe_win)

#output the sliding-windows file
write.table(hwe_win, paste0(FILE,".slidingwindows"), row.names=FALSE, quote=FALSE, sep="\t")

