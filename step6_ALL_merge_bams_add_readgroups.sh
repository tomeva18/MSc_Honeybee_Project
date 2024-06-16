#!/bin/bash

#SBATCH --job-name=merge_bams
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --account=bisc028524
#SBATCH --output=merge_bams_%j.log
#SBATCH --error=merge_bams_%j.err

# Load modules
module load apps/picard/2.27.5
module load apps/samtools/1.16
module load lang/python/mamba/1.4.2

# Print the hostname and date for debugging purposes
echo "Running on $(hostname)"
date

# Run the Python script
python merge_bams_add_readgroups.py
