#!/bin/bash

#SBATCH --job-name=ALL_fastqc
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --account=bisc028524

## Load module
module load apps/fastqc/0.11.9

## Clarify scratch 
scratch=/user/work/hg20812/project1

mkdir ALL_fastqc_before_trim
for i in $(ls $scratch/ALL_fastq_files_rnaseq); do fastqc -t 24 $scratch/ALL_fastq_files_rnaseq/$i -o ALL_fastqc_before_trim; done