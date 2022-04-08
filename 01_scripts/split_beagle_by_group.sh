#!/bin/bash
#SBATCH -J "split_beagle"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=XXX
#SBATCH --time=1-00:00
#SBATCH --mem=100G
# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR


###this script will take a beagle cut it for the list of individuals given

#to edit
GROUP=MyWantedGroup #the name of the subgroup on which we are making the subset 
#-> it should be associated with a file 02_info/"$GROUP"bam.filelist
BAM_GROUP=02_info/"$GROUP"bam.filelist

#if you want to make the bamlist for the group by joining several populations
#POP1=
#POP2=
#POP3=
#...
#cat 02_info/"$POP1"bam.filelist 02_info/"$POP2"bam.filelist 02_info/"$POP3"bam.filelist > 02_info/"$GROUP"bam.filelist

#standard variables - no need to change
source 01_scripts/01_config.sh
BAM_ALL=02_info/bam.filelist # where and what is the name of the big bam list on which the beagle has been built

#we will subset into a new folder
mkdir 03_saf_maf_gl_all/subset_beagle

#1 translate bamfilelist into beagle coordinates
#and make a file with the position of the individuals for each subgroup in the bigger beagle
Rscript 01_scripts/Rscripts/subset_ind_coordinates.r "$GROUP" "$BAM_ALL"


#2 subset beagle into a new folder

INFILE=03_saf_maf_gl_all/all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".beagle.gz
OUTFILE=03_saf_maf_gl_all/subset_beagle/"$GROUP"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".beagle.gz

k=$(cat 03_saf_maf_gl_all/subset_beagle/"$GROUP"_column.pos)

gunzip -c $INFILE | head | cut -f $k | gzip -c > $OUTFILE

