#this code will use the covariance made along window of the genome 
#then it applies analysis from the losctruct packages as described in
#Li, H., & Ralph, P. (2019). Local PCA shows how the effect of population structure differs along the genome. 
#Genetics, 211(1), 289-304.

# I have also added calculation of hopkins index to determine within which windows a trend of clustering is observed
# And it also refers to the global pca to see which window have a PC1 that correlates with the first three PCs of the global PCA

# code R
library(lostruct)
library(clustertend)

#import arguments for file naming
argv <- commandArgs(T)
MIN_MAF <- argv[1]
PERCENT_IND <- argv[2]
WINDOW_SIZE<- argv[3]
INFO<-argv[4]
MAX_DEPTH_FACTOR<-argv[5]
#change on which LG one want to work dependaing on the gneome
if (WINDOW_SIZE<=1000)
{LGs<-c("LG1","LG2", "LG3","LG4","LG5","LG6")}
if (WINDOW_SIZE>1000)
{LGs<-c("LG1","LG2", "LG3","LG4","LG5")}

#this alternativez allows looking at correlation with specific (continuous variables) givne in a x_info file whose path is given as arguments
#if FALSE it will be skipped
if (INFO!=FALSE)
{
	x_info<-read.table(INFO, header=T)
	n<-dim(x_info)[2]
}else{n=0}

#import info on samples such as coordinates in the big PCA to see which portions of the genome are associated
# with the structure seen on the first 3 PCs of the PCa on the whole genome
info_pc<-read.table(paste0("04_pca/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,".pca"))

#import the list of covariance matrix by window
window<-read.table(paste0("10_pca_by_window/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"/",WINDOW_SIZE,"/analyse_by_window/list_window.txt"), 
header=F)
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

#work on each LG to run lostruct, calculate hopkins and look at correlation within each window PC1 and given variables
#choose nb of PC kept for lostruct analysis
N_PC<-2
#initialise loop
window_result<-matrix(ncol=12)
window_result_bis<-matrix(ncol=12+n)

for (j in 1:length(LGs))
{
	print(paste("working on LG",j))
	
	window_order_j<-window_order[window_order$LG==LGs[j],]
	
	window_order_j$h_stat<-vector(length=dim(window_order_j)[1])
	window_order_j$corr_pc1<-vector(length=dim(window_order_j)[1])
	window_order_j$corr_pc2<-vector(length=dim(window_order_j)[1])
	window_order_j$corr_pc3<-vector(length=dim(window_order_j)[1])

	#import the cov matrix into a list
	window_covmats<-list()
	for (i in 1: dim(window_order_j)[1])
	{
		print(paste("import covariance matrix for LG",j,"window",i))
		window_covmats[[i]]<-as.matrix(read.table(paste0("10_pca_by_window/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"/",WINDOW_SIZE,"/cov_by_window/",window_order_j$file_name[i])))
	
	}
	
	#if window size is longer than the LG there is just one point and the analysis crash
	if(length(window_covmats)>1)
	{
		#run lostruct fonctions to calculate PC axis & then distances between pca made on each window
		window_eigs <- sapply(window_covmats, function (cm) cov_pca(k=N_PC, covmat=cm))
		window_dist<- pc_dist(t(window_eigs), npc=N_PC)
		
		#run a MDS to visualise the distribution of window
		mds_axe<-cmdscale(window_dist)
		window_order_j$mds1_byLG<-mds_axe[,1]
		window_order_j$mds2_byLG<-mds_axe[,2]
		jpeg(paste0("10_pca_by_window/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"/",WINDOW_SIZE,"/analyse_by_window/mds_",LGs[j],".jpg"))
		plot(mds_axe[,1], mds_axe[,2], pch=20, main=LGs[j])
		dev.off()
		
		#calculate for each window hopkins statistics
		#as well as correlation with the first 3 PCs observed in the global PCA on the whole genome
		
		N_ind<-(dim(window_eigs)[1] - 1 - N_PC)/N_PC
		for (i in 1: dim(window_order_j)[1])
		{
			print(paste("calculating hopkins for LG", j,"window",i))
			PC_mat_i<- cbind(window_eigs[4:(3+N_ind),i],window_eigs[(4+N_ind):(3+2*N_ind),i])
			window_order_j$h_stat[i]<-hopkins(PC_mat_i,N_ind/2, header=FALSE)
			window_order_j$corr_pc1[i]<-abs(cor.test(PC_mat_i[,1], info_pc[,1])$estimate)
			window_order_j$corr_pc2[i]<-abs(cor.test(PC_mat_i[,1], info_pc[,2])$estimate)
			window_order_j$corr_pc3[i]<-abs(cor.test(PC_mat_i[,1], info_pc[,3])$estimate)
		}
		
		#make a figure with MDS score, hopkins and correlation
		
		jpeg(paste0("10_pca_by_window/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"/",WINDOW_SIZE,"/analyse_by_window/mds_bypos_h_",LGs[j],".jpg"), height=1000, width=600)
		par(mfrow=c(3,2))
		plot((window_order_j$start+window_order_j$stop)/2, mds_axe[,1], main=LGs[j], pch=20)
		plot((window_order_j$start+window_order_j$stop)/2, window_order_j$corr_pc1, main="corr_PC1", pch=20, ylim=c(0,1))
		plot((window_order_j$start+window_order_j$stop)/2, mds_axe[,2], pch=20)
		plot((window_order_j$start+window_order_j$stop)/2, window_order_j$corr_pc2, main="corr_PC2", pch=20, ylim=c(0,1))
		plot((window_order_j$start+window_order_j$stop)/2, window_order_j$h_stat, main="hopkins", pch=20)
		abline(h=0.10)
		plot((window_order_j$start+window_order_j$stop)/2, window_order_j$corr_pc3, main="corr_PC3", pch=20,ylim=c(0,1))
		dev.off()
		
		#save all the information from LG j into a bigger table
		colnames(window_result)<-colnames(window_order_j)
		window_result<-rbind(window_result,window_order_j) 
		
		#if we give it an additional file to correlate with the analysis will run into that part
		window_order_j_bis<-window_order_j
		if (n>0)
		{
			print("looking at correlations with the variable given")
			jpeg(paste0("10_pca_by_window/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"/",WINDOW_SIZE,"/analyse_by_window/corr_with_info_",LGs[j],".jpg"), height=1200, width=600)
			par(mfrow=c(n/2,2))
			for (k in 1 : n)
			{
				corr_k<-vector(length= dim(window_order_j_bis)[1])
				for (i in 1: dim(window_order_j_bis)[1])
				{	
				corr_k[i]<-abs(cor.test(window_eigs[4:(3+N_ind),i], x_info[,k])$estimate)
				}
				window_order_j_bis<-cbind(window_order_j_bis, corr_k)
				colnames(window_order_j_bis)[k+12]<-paste0("corr_",colnames(x_info)[k])
				plot((window_order_j_bis$start+window_order_j_bis$stop)/2, corr_k, main=paste("correlation with",colnames(x_info)[k],"along", LGs[j]), pch=20, ylim=c(0,1))
			}
			dev.off()
		}
		colnames(window_result_bis)<-colnames(window_order_j_bis)
		#save all the information from LG j into a bigger table
		window_result_bis<-rbind(window_result_bis,window_order_j_bis)
	}
}

write.table(as.matrix(window_result_bis[2:(dim(window_result_bis)[1]),]),paste0("10_pca_by_window/all_maf",MIN_MAF,"_pctind",PERCENT_IND,"_maxdepth",MAX_DEPTH_FACTOR,"/",WINDOW_SIZE,"/analyse_by_window/window_result.txt"), quote=F,row.names=F)  

