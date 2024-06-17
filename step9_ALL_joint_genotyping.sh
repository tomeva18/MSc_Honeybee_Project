#!/bin/bash

#SBATCH --job-name=GenotypeGVCF
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --account=bisc028524

# Specify directories
work_dir=/user/work/hg20812/project1
output_dir=${work_dir}/joint_genotyping
reference=/user/work/hg20812/project/individual_test/STAR_mapping/Reference/Apis_mellifera.Amel_HAv3.1.dna.toplevel.fa

# Create output directory if it doesn't exist
mkdir -p ${output_dir}

# Load necessary modules
module load apps/gatk/4.1.9

# Run GATK GenotypeGVCFs
gatk GenotypeGVCFs \
    -R ${reference} \
    -V gendb://${work_dir}/my_database \
    -O ${output_dir}/joint_calls.vcf.gz
