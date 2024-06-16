#!/bin/bash

#SBATCH --job-name=DBImport
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --account=bisc028524

# Specify directories
work_dir=/user/work/hg20812/project1
tmp=${work_dir}/tmp

# Create directory
mkdir -p ${tmp}

# Load module
module load apps/gatk/4.1.9

# Check if modules are loaded correctly
module list

# Run GATK GenomicsDBImport
gatk GenomicsDBImport \
    --genomicsdb-workspace-path ${work_dir}/my_database \
    -L ${work_dir}/intervals.list \
    --sample-name-map ${work_dir}/sample_map2.txt \
    --tmp-dir ${tmp} \
    --reader-threads 12
