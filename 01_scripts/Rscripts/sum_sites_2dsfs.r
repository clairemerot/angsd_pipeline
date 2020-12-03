# code R to take the sfs made on a subsample of sites in a usable format for subsequent analyses
argv <- commandArgs(T)
file<-argv[1]


sfs<-read.table (paste0(file))

sfs.sum<-colSums(sfs)
write.table(rbind(sfs.sum),  quote=F, col.names=F, row.names=F,paste0(file, ".2dsfs"))
