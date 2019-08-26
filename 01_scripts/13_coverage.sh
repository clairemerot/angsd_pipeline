#!/bin/bash
#SBATCH -J "13_coverage"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOURMAIL
#SBATCH --time=7-00:00
#SBATCH --mem=50G

###this script will work on all bamfiles and calculate coverage at each position
#maybe edit
NB_CPU=1 #change accordingly in SLURM header
BAM_FILELIST=02_info/bam.filelist

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

mkdir 13_coverage

module load samtools
ulimit -S -n 2048

####Calculate coverage for individual in bam.filelist
samtools depth -a -f $BAM_FILELIST | gzip - > 13_coverage/cov_trial.gz

