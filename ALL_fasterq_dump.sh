#!/bin/bash

#SBATCH --job-name=test_fasterq_dump
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --account=bisc028524

# Navigate to the project directory
cd /user/work/hg20812/project1

# Create the output directory
mkdir -p ./ALL_fastq_files_rnaseq

# Loop through each accession number in the list and process it
for i in $(cat test_accession_list.txt)
do
    fasterq-dump --outdir ./ALL_fastq_files_rnaseq --split-files $i
    gzip ./ALL_fastq_files_rnaseq/${i}*.fastq
done


