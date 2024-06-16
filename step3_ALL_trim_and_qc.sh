#!/bin/bash

#SBATCH --job-name=test_fasterq_dump
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --account=bisc028524

# Load fastp module if available or ensure the path is correct
fastp_path="/user/work/hg20812/project1/fastp"  # Uncomment and set this if fastp is not a module but installed in a specific location

# Set directories
files_concat="/user/work/hg20812/project1/ALL_fastq_files_rnaseq"
work_dir="/user/work/hg20812/project1"

# Create output directories
mkdir -p ${work_dir}/ALL_trimmed_reads
mkdir -p ${work_dir}/ALL_fastp_reports
mkdir -p ${work_dir}/ALL_fastqc_after_trim

# Run fastp on .fastq files
for i in ${files_concat}/*_1.fastq; do
    r1=$(basename ${i})
    name=${r1%_1.fastq}
    r2="${name}_2.fastq"

    fastp -w 16 -l 40 \
    --cut_front --cut_tail --trim_poly_x \
    --html ${work_dir}/fastp_reports/${name}.html --json ${work_dir}/fastp_reports/${name}.json \
    -i ${files_concat}/${r1} -I ${files_concat}/${r2} \
    -o ${work_dir}/trimmed_reads/${r1} -O ${work_dir}/trimmed_reads/${r2}
done

# Run fastqc on trimmed read data
for i in ${work_dir}/trimmed_reads/*.fastq; do
    fastqc -t 8 ${i} -o ${work_dir}/fastqc_after_trim  # Adjust -t 8 as per available CPUs
done
