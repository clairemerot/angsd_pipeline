#!/bin/bash
#SBATCH -J "12C_simpl_LD"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p large
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=21-00:00
#SBATCH --mem=30G

###this script simplifies the Ld by taking a valu by block of Xkb

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#split the Ld into chr and inter chr
#If the file is big it can take long. You can split into several scripts that you will run in parallel (as this use only 1CPU)

LD_FILE=12_ngsLD/all.ld.gz
#edit with the name of your 2 chromosomes
CHR1=NC_027307.1
CHR2=NC_027328.1


gunzip -c "$LD_FILE" | awk '{if ($3 == "inf") {print}}' | gzip > "$LD_FILE"_inter.ld.gz

#this is for the 1st chr
gunzip -c "$LD_FILE" | grep ^$CHR1 | awk '{if ($3 != "inf") {print}}' | gzip > "$LD_FILE"_"$CHR1".ld.gz 

#this is for the 2nd chr
gunzip -c "$LD_FILE" | grep ^$CHR2 | awk '{if ($3 != "inf") {print}}' | gzip > "$LD_FILE"_"$CHR2".ld.gz 


#run Eric's script to keep quartile values per windows of 1000kb, 500 kb, 250 kb
cat 12_ngsLD/header.gz "$LD_FILE"_inter.ld.gz > "$LD_FILE"_inter_header.ld.gz

#python3 01_scripts/utility_scripts/ld_by_blocks_optimized_gzinput.py "$LD_FILE"_inter_header.ld.gz 1000 "$LD_FILE"_inter_by_1000.ld.gz
python3 01_scripts/utility_scripts/ld_by_blocks_optimized_gzinput.py "$LD_FILE"_inter_header.ld.gz 500 "$LD_FILE"_inter_by_500.ld.gz
#python3 01_scripts/utility_scripts/ld_by_blocks_optimized_gzinput.py "$LD_FILE"_inter_header.ld.gz 250 "$LD_FILE"_inter_by_250.ld.gz

#make space
rm "$LD_FILE"_inter.ld.gz


#run Eric's script to keep quartile values per windows of 1000kb, 500 kb, 250 kb
cat 12_ngsLD/header.gz "$LD_FILE"_"$CHR1".ld.gz > "$LD_FILE"_"$CHR1"_header.ld.gz

#python3 01_scripts/utility_scripts/ld_by_blocks_optimized.py "$LD_FILE"_"$CHR1"_header.ld 1000 "$LD_FILE"_"$CHR1"_by_1000.ld
python3 01_scripts/utility_scripts/ld_by_blocks_optimized.py "$LD_FILE"_"$CHR1"_header.ld 500 "$LD_FILE"_"$CHR1"_by_500.ld
#python3 01_scripts/utility_scripts/ld_by_blocks_optimized.py "$LD_FILE"_"$CHR1"_header.ld 250 "$LD_FILE"_"$CHR1"_by_250.ld

#make space
rm "$LD_FILE"_"$CHR1".ld.gz



#run Eric's script to keep quartile values per windows of 1000kb, 500 kb, 250 kb
cat 12_ngsLD/header.gz "$LD_FILE"_"$CHR2".ld.gz > "$LD_FILE"_"$CHR2"_header.ld

#python3 01_scripts/utility_scripts/ld_by_blocks_optimized.py "$LD_FILE"_"$CHR2"_header.ld 1000 "$LD_FILE"_"$CHR2"_by_1000.ld
python3 01_scripts/utility_scripts/ld_by_blocks_optimized.py "$LD_FILE"_"$CHR2"_header.ld 500 "$LD_FILE"_"$CHR2"_by_500.ld
#python3 01_scripts/utility_scripts/ld_by_blocks_optimized.py "$LD_FILE"_"$CHR2"_header.ld 250 "$LD_FILE"_"$CHR2"_by_250.ld

#make space
rm "$LD_FILE"_"$CHR2".ld.gz

