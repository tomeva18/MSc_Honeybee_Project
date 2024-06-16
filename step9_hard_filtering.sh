#!/bin/bash

#SBATCH --job-name=HardFilterVariants
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --account=bisc028524

# Specify directories
work_dir=/user/work/hg20812/project1
output_dir=${work_dir}/filtered_variants
joint_genotyped_vcf=${work_dir}/joint_calls.vcf.gz

# Create output directory if it doesn't exist
mkdir -p ${output_dir}

# Load GATK module
module load apps/gatk/4.1.9

# Split into SNPs and indels
gatk SelectVariants \
    -V ${joint_genotyped_vcf} \
    -select-type SNP \
    -O ${output_dir}/snps.vcf.gz

gatk SelectVariants \
    -V ${joint_genotyped_vcf} \
    -select-type INDEL \
    -O ${output_dir}/indels.vcf.gz

# Filter SNPs
gatk VariantFiltration \
    -V ${output_dir}/snps.vcf.gz \
    --cluster-window-size 35 \
    --cluster-size 3 \
    -filter "QD < 2.0" -filter-name "QD2" \
    -filter "FS > 60.0" -filter-name "FS60" \
    -filter "MQ < 40.0" -filter-name "MQ40" \
    -filter "MQRankSum < -12.5" -filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" -filter-name "ReadPosRankSum-8" \
    -filter "SOR > 3.0" -filter-name "SOR3" \
    -O ${output_dir}/snps_filtered.vcf.gz

# Filter Indels
gatk VariantFiltration \
    -V ${output_dir}/indels.vcf.gz \
    --cluster-window-size 35 \
    --cluster-size 3 \
    -filter "QD < 2.0" -filter-name "QD2" \
    -filter "FS > 200.0" -filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" -filter-name "ReadPosRankSum-20" \
    -O ${output_dir}/indels_filtered.vcf.gz