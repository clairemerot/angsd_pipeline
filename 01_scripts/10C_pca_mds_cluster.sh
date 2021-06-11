#!/bin/bash
#SBATCH -J "10C_pca_cluster"
#SBATCH -o log_%j
#SBATCH -c 3 
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=7-00:00
#SBATCH --mem=100G

###this script will work on output from the lostruct analysis and implement the method form Huang et al 2020 to select significan tclusters of windows
#it requires two Rscripts/make_pca_simple and Rscripts/outliers_clusters_mds_allLG

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR
ulimit -S -n 2048

#prepare variables - avoid to modify
source 01_scripts/01_config.sh

#maybe edit
NB_CPU=3 #change accordingly in SLURM header
window_size=100
BAM_LIST=02_info/bam.filelist

#somewhere where the whole beagle file is uncompressed 
path_beagle_all=10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/pca_by_chosen_area
#if the all.beagle is not ready modify the path and uncomment the following line
#gunzip -c 03_saf_maf_gl_all/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".beagle.gz > "$path_beagle_all"/all.beagle

#parameters to group together into a cluster the outliers windows
max_MDS=11 #nb of mds on whcih to work
sd_lim=4 #limit for being outlier
x_between=20 #merge window with this number of windows between them
min_n_window=5 #remove cluster with less than this number of window



INPUT_FOLDER=10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$window_size"/

mkdir 10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$window_size"/outlier_mds
mkdir 10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$window_size"/outlier_mds/mds_plots
mkdir 10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$window_size"/outlier_mds/cluster_position

	

#this script will pull out the clusters of windows that are outliers on each MDS axis
echo "running Rscript on mds results to pull out clusters of outliers windows"
Rscript 01_scripts/Rscripts/outliers_clusters_mds_allLG.R "$INPUT_FOLDER" "$max_MDS" "$sd_lim" "$x_between" "$min_n_window"


mkdir 10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$window_size"/outlier_mds/pca_clusters"$x_between"_filter"$min_n_window"_sdLim"$sd_lim"
path=10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$window_size"/outlier_mds/pca_clusters"$x_between"_filter"$min_n_window"_sdLim"$sd_lim"

#list the file with position of clusters (a file with two column with the start and stop of the wanted regions in one cluster)
ls -1 10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/"$window_size"/outlier_mds/cluster_position/*"$x_between"_filter"$min_n_window"_sdLim"$sd_lim".pos | 
    while read j
    do
	POS_FILE="$j"
	echo "$POS_FILE"
	CLUSTER=$(basename "${j%.pos}")
	echo "working on $CLUSTER"		
	#count the nb of portions in the cluster
	num_pca=$(wc -l "$POS_FILE" | cut -d " " -f 1)
	echo "this cluster has $num_pca groups of windows"
	
	#make a header for the beagle of the joint cluster
	head -n 1 "$path_beagle_all"/all.beagle > "$path"/"$CLUSTER".beagle
	
	##split the beagle files and concatenate them below to group the diff portions into a single beagle		
	for i in $(seq $num_pca)
	do
	start=$(cat "$POS_FILE" | head -"$i" | tail -1 | cut -f 1)
	stop=$(cat "$POS_FILE" | head -"$i" | tail -1 | cut -f 2)
	
	echo "splitting beagle $start $stop"
	python2 01_scripts/utility_scripts/beagle_positions_extract_to_file.py "$path_beagle_all"/all.beagle $start $stop "$path"/"$start"-"$stop".beagle
	#remove header
	tail -n +2 "$path"/"$start"-"$stop".beagle > "$path"/"$start"-"$stop"_wo_header.beagle
	#concatenateq
	cat "$path"/"$CLUSTER".beagle "$path"/"$start"-"$stop"_wo_header.beagle > "$path"/"$CLUSTER"_intermediate.beagle
	#save into a header for the next concatenate
	cp "$path"/"$CLUSTER"_intermediate.beagle "$path"/"$CLUSTER".beagle
	
	done
	
	N_snp=$(wc -l "$path"/"$CLUSTER".beagle | cut -d " " -f 1)
	gzip "$path"/"$CLUSTER".beagle 
	rm "$path"/*.beagle
	
	#run pcangsd to get the covariance matrix
	echo "pcangsd on $CLUSTER"
	python2 $PCA_ANGSD_PATH/pcangsd.py -threads $NB_CPU -beagle "$path"/"$CLUSTER".beagle.gz -o "$path"/"$CLUSTER"_"$N_snp"snp

	#make the pca, a plot with kmeans group attribution
	INPUT="$path"/"$CLUSTER"_"$N_snp"snp.cov
	echo "pca with R on $INPUT"
	Rscript  01_scripts/Rscripts/make_pca_simple.r "$INPUT" "$BAM_LIST"
	
	gzip "$path"/"$CLUSTER"_"$N_snp"snp.cov
	
done

