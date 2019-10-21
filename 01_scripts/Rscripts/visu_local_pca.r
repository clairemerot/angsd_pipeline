#this script is to visualise the results of pca done by scaffold, it will output a figure coloured by the variables 
# given - here it is plan for 4 variables)
#it will also look at correlation with PC axes done on the main global whole-genome pca


#import arguments for file naming
argv <- commandArgs(T)
MIN_MAF <- argv[1]
PERCENT_IND <- argv[2]
WINDOW_SIZE<- argv[3]
PATH_INFO_FILE<-argv[4]
MAX_DEPTH_FACTOR<-argv[5]
var1<-as.numeric(argv[6])
var2<-as.numeric(argv[7])
var3<-as.numeric(argv[8])
var4<-as.numeric(argv[9])


#upload the list of cov matrix 
window_file<-read.table(paste0("10_pca_by_window/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"/", WINDOW_SIZE, "/pdf_by_window/window_cov_files.txt"))
#head(window_file)

#upload an info file
info<-read.table(PATH_INFO_FILE, header=T)
print("plotting pca coloured by variable:")
colnames(info)[c(var1,var2,var3,var4)]
#head(info)
 
#upload the main pc coordinates
info_pc<-read.table(paste0("04_pca/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,".pca"))
#head (info_pc)

window_order<-matrix(nrow=dim(window_file)[1], ncol=5)
colnames(window_order)<-c("window", "corr_pc1", "corr_pc2", "corr_pc3", "var_pc1")

for (i in 1 : dim(window_file)[1])
{
#read the cov matrix
print(paste("doing pca and correlation analyses for ", window_file[i,1]))
cov.mati<-read.table(paste0("10_pca_by_window/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"/", WINDOW_SIZE, "/cov_by_window/",window_file[i,1]))
#make pca
pca.mati<-prcomp(cov.mati)
#calculate eigen value
eig1<-round(pca.mati$sdev[1]/sum(pca.mati$sdev),2)*100
eig2<-round(pca.mati$sdev[2]/sum(pca.mati$sdev),2)*100
eig3<-round(pca.mati$sdev[3]/sum(pca.mati$sdev),2)*100
eig4<-round(pca.mati$sdev[4]/sum(pca.mati$sdev),2)*100

#calculate correlation with main pca
window_order[i,1]<-window_file[i,1]
window_order[i,2]<-abs(cor.test(pca.mati$x[,1], info_pc[,1])$estimate)
window_order[i,3]<-abs(cor.test(pca.mati$x[,1], info_pc[,2])$estimate)
window_order[i,4]<-abs(cor.test(pca.mati$x[,1], info_pc[,3])$estimate)
window_order[i,5]<-eig1

#make a figure colouring by different factor
jpeg(paste0("10_pca_by_window/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"/", WINDOW_SIZE, "/pdf_by_window/", window_file[i,1], ".jpg"), height=1000, width=600)
par(mfrow=c(4,2))
#to colour the pca by one factor for pc1-2  & pc3-4
plot(pca.mati$x[,1],pca.mati$x[,2], pch=20, col=info[,var1], main=window_file[i,1], xlab=paste("PC1",eig1,"%"), ylab=paste("PC2",eig2,"%"))
plot(pca.mati$x[,3],pca.mati$x[,4], pch=20, col=info[,var1], main=window_file[i,1], xlab=paste("PC3",eig3,"%"), ylab=paste("PC4",eig4,"%"))
legend("topleft", legend=levels(info[,var1]), col=as.numeric(as.factor(levels(info[,var1]))), pch=20, cex=0.8)
plot(pca.mati$x[,1],pca.mati$x[,2], pch=20, col=info[,var2], main=window_file[i,1], xlab=paste("PC1",eig1,"%"), ylab=paste("PC2",eig2,"%"))
plot(pca.mati$x[,3],pca.mati$x[,4], pch=20, col=info[,var2], main=window_file[i,1], xlab=paste("PC3",eig3,"%"), ylab=paste("PC4",eig4,"%"))
legend("topleft", legend=levels(info[,var2]), col=as.numeric(as.factor(levels(info[,var2]))), pch=20, cex=0.8)
plot(pca.mati$x[,1],pca.mati$x[,2], pch=20, col=info[,var3], main=window_file[i,1], xlab=paste("PC1",eig1,"%"), ylab=paste("PC2",eig2,"%"))
plot(pca.mati$x[,3],pca.mati$x[,4], pch=20, col=info[,var3], main=window_file[i,1], xlab=paste("PC3",eig3,"%"), ylab=paste("PC4",eig4,"%"))
legend("topleft", legend=levels(info[,var3]), col=as.numeric(as.factor(levels(info[,var3]))), pch=20, cex=0.8)
plot(pca.mati$x[,1],pca.mati$x[,2], pch=20, col=info[,var4], main=window_file[i,1], xlab=paste("PC1",eig1,"%"), ylab=paste("PC2",eig2,"%"))
plot(pca.mati$x[,3],pca.mati$x[,4], pch=20, col=info[,var4], main=window_file[i,1], xlab=paste("PC3",eig3,"%"), ylab=paste("PC4",eig4,"%"))
legend("topleft", legend=levels(info[,var4]), col=as.numeric(as.factor(levels(info[,var4]))), pch=20, cex=0.8)

dev.off()
}

write.table(window_order, paste0("10_pca_by_window/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"/", WINDOW_SIZE, "/pdf_by_window/corr_mainpca_window.txt"), quote=F, row.names=F, sep="\t")