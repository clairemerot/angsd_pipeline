#!/bin/bash
#SBATCH -J "10C_pca_chosen"
#SBATCH -o log_%j
#SBATCH -c 3 
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=X
#SBATCH --time=7-00:00
#SBATCH --mem=100G

###this script will work on all bamfiles and calculate saf, maf & genotype likelihood
#maybe edit
NB_CPU=3 #change accordingly in SLURM header

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

ulimit -S -n 2048

#prepare variables - avoid to modify
source 01_scripts/01_config.sh

##step1-split the beagle files
#this requires a file with two column with the start and stop of the wanted region
#eg LG1_125 LG1_456
CHOSEN_FILE=10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/pca_by_chosen_area/outlier_mds_100/inv_position_final.pos


num_pca=$(wc -l "$CHOSEN_FILE" | cut -d " " -f 1)
BAM_LIST=02_info/bam.filelist
path=10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/pca_by_chosen_area/outlier_mds_100
path_beagle_all=10_pca_by_window/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"/pca_by_chosen_area

#if the all.beagle is not ready
#gunzip -c 03_saf_maf_gl_all/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".beagle.gz > "$path_beagle_all"/all.beagle


for i in $(seq $num_pca)
	do
	start=$(cat "$CHOSEN_FILE" | head -"$i" | tail -1 | cut -f 1)
	stop=$(cat "$CHOSEN_FILE" | head -"$i" | tail -1 | cut -f 2)

	#echo "splitting beagle $start $stop"
	python2 01_scripts/utility_scripts/beagle_positions_extract_to_file.py "$path_beagle_all"/all.beagle $start $stop "$path"/"$start"-"$stop".beagle
	gzip "$path"/"$start"-"$stop".beagle

	echo "pcangsd"
	python2 $PCA_ANGSD_PATH/pcangsd.py -threads $NB_CPU -beagle "$path"/"$start"-"$stop".beagle.gz -o "$path"/"$start"-"$stop"

	
	INPUT="$path"/"$start"-"$stop".cov
	echo "pca with R on $INPUT"
	Rscript  01_scripts/Rscripts/make_pca_simple.r "$INPUT" "$BAM_LIST"
done

