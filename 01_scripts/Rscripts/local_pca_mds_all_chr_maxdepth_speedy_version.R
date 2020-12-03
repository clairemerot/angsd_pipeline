#this code will use the covariance made along window of the genome
#CAUTION this version works on all the genome and not chromosome by chr 
#then it applies analysis from the losctruct packages as described in
#Li, H., & Ralph, P. (2019). Local PCA shows how the effect of population structure differs along the genome. 
#Genetics, 211(1), 289-304.

# And it also refers to the global pca to see which window have a PC1 that correlates with the first three PCs of the global PCA

# code R
library(lostruct)

#import arguments for file naming
argv <- commandArgs(T)
MIN_MAF <- argv[1]
PERCENT_IND <- argv[2]
WINDOW_SIZE<- argv[3]
MAX_DEPTH_FACTOR<-argv[4]
N_ind<-as.numeric(argv[5])
N_PC<-as.numeric(argv[6])
N_MDS<-as.numeric(argv[7])


#choose nb of PC kept for lostruct analysis and the nb of MDS
#N_PC<-2
#N_MDS<-50

#how many individuals? = dimension of the covariance matrix
#N_ind<-1446
#window_covmats_1<-as.matrix(read.table(paste0("10_pca_by_window/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"/",WINDOW_SIZE,"/cov_by_window/",window_order$file_name[1])))
#N_ind<-dim(window_covmats_1)[1]


#import the list of covariance matrix by window
window<-read.table(paste0("10_pca_by_window/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"/",WINDOW_SIZE,"/analyse_by_window/list_window.txt"),header=F)
colnames(window)<-"file_name"

#make a table to caracterize the different windows
window$LG<-vector(length=dim(window)[1])
window$pos<-vector(length=dim(window)[1])
window$n_snp<-vector(length=dim(window)[1])
window$start<-vector(length=dim(window)[1])
window$stop<-vector(length=dim(window)[1])

for (i in 1 :dim(window)[1] )
{
window$LG[i]<-unlist(strsplit(as.character(window$file_name[i]), split="_"))[2]
window$pos[i]<-unlist(strsplit(as.character(window$file_name[i]), split="_"))[3]
window$n_snp[i]<-as.numeric(unlist(strsplit(as.character(window$file_name[i]), split="_"))[4])
window$start[i]<-as.numeric(unlist(strsplit(as.character(window$pos[i]), split="-"))[1])
window$stop[i]<-as.numeric(unlist(strsplit(as.character(window$pos[i]), split="-"))[2])
}
#re-order the window correctly
window_order<-window[order(window$LG, window$start),]
head(window_order)


####work on all LG to run lostruct


#prepare an empty matrix with nrow= nb of window, and ncol= sumofsquares + nb of eigen value + dim of all eigen vectors (Nb individuals * N_Pc)
window_eigs<-matrix(nrow=dim(window_order)[1], ncol= (1+N_PC+ (N_ind*N_PC)))

#import each covariance matrix and make the eigen decomposition
for (i in 1: dim(window_order)[1])
{
	print(paste("import covariance matrix for window and get eigenvectors",i))
	window_covmats_i<-as.matrix(read.table(paste0("10_pca_by_window/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"/",WINDOW_SIZE,"/cov_by_window/",window_order$file_name[i])))
	window_eigs [i,]<- cov_pca(k=N_PC, covmat=window_covmats_i)
}

#save the big matrix just in case
print("saving matrix of eigenvectors")
write.table(window_eigs,paste0("10_pca_by_window/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"/",WINDOW_SIZE,"/analyse_by_window/window_eigenvector_all_LG.txt"), quote=F,row.names=F)  
 
#run lostruct fonctions to calculate  distances between pca made on each window
print(paste("calculating distance between window PCA with NPC=",N_PC))
window_dist<- pc_dist(window_eigs, npc=N_PC)
	
#run a MDS to visualise the distribution of window
print(paste("run MDS on distance matrix and will keep N_MDS=",N_MDS))
mds_axe<-cmdscale(window_dist, eig=TRUE, k=N_MDS)
window_order_mds<-cbind(window_order, mds_axe$points)

#plot the 12 first MDS
print("plotting MDS axis")
jpeg(paste0("10_pca_by_window/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"/",WINDOW_SIZE,"/analyse_by_window/mds_all_LG.jpg"))
par(mfrow=c(3,2))
plot(mds_axe$points[,1], mds_axe$points[,2], pch=20)
plot(mds_axe$points[,3], mds_axe$points[,4], pch=20)
plot(mds_axe$points[,5], mds_axe$points[,6], pch=20)
plot(mds_axe$points[,7], mds_axe$points[,8], pch=20)
plot(mds_axe$points[,9], mds_axe$points[,10], pch=20)
plot(mds_axe$points[,11], mds_axe$points[,12], pch=20)
dev.off()

print("writing mds results")
write.table(window_order_mds,paste0("10_pca_by_window/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"/",WINDOW_SIZE,"/analyse_by_window/window_mds_all_LG_values.txt"), quote=F,row.names=F)  
write.table(mds_axe$eig,paste0("10_pca_by_window/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"/",WINDOW_SIZE,"/analyse_by_window/window_mds_all_LG_eig.txt"), quote=F,row.names=F)  


## if wanted,  one can also look at correlation with the PCs observed in the global PCA on the whole genome
#import info on samples such as coordinates in the big PCA to see which portions of the genome are associated
# with the structure seen on the first 3 PCs of the PCa on the whole genome
info_pc<-read.table(paste0("04_pca/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,".pca"))[,1:20]


N_ind<-(dim(window_eigs)[1] - 1 - N_PC)/N_PC
for (i in 1: dim(window_order)[1])
{
	print(paste("calculating corr pca for window",i))
	PC_mat_i<- cbind(window_eigs[4:(3+N_ind),i],window_eigs[(4+N_ind):(3+2*N_ind),i])
	window_order$eig_pc1[i]<-window_eigs[1,i] #not sure which number is eigenvalue or sum of square and how to get variance from there
	window_order$eig_pc2[i]<-window_eigs[2,i]
	window_order$corr_pc1[i]<-abs(cor.test(PC_mat_i[,1], info_pc[,1])$estimate)
	window_order$corr_pc2[i]<-abs(cor.test(PC_mat_i[,1], info_pc[,2])$estimate)
	window_order$corr_pc3[i]<-abs(cor.test(PC_mat_i[,1], info_pc[,3])$estimate)
	window_order$corr_pc4[i]<-abs(cor.test(PC_mat_i[,1], info_pc[,4])$estimate)
	window_order$corr_pc5[i]<-abs(cor.test(PC_mat_i[,1], info_pc[,5])$estimate)
	window_order$corr_pc6[i]<-abs(cor.test(PC_mat_i[,1], info_pc[,6])$estimate)
	window_order$corr_pc7[i]<-abs(cor.test(PC_mat_i[,1], info_pc[,7])$estimate)
	window_order$corr_pc8[i]<-abs(cor.test(PC_mat_i[,1], info_pc[,8])$estimate)
    window_order$corr_pc9[i]<-abs(cor.test(PC_mat_i[,1], info_pc[,9])$estimate)
	window_order$corr_pc10[i]<-abs(cor.test(PC_mat_i[,1], info_pc[,10])$estimate)
	window_order$corr_pc11[i]<-abs(cor.test(PC_mat_i[,1], info_pc[,11])$estimate)
    window_order$corr_pc12[i]<-abs(cor.test(PC_mat_i[,1], info_pc[,12])$estimate)
}
write.table(window_order,paste0("10_pca_by_window/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"/",WINDOW_SIZE,"/analyse_by_window/window_corr_pc_all_LG.txt"), quote=F,row.names=F)  

