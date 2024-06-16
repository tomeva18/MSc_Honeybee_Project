#!/bin/bash

#SBATCH --job-name=ALL_star_mapping
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G
#SBATCH --account=bisc028524

# Load module
module load apps/star/2.7

# Number of permitted open files needs to be increased
ulimit -n 4096

# Set location of scratch space
scratch=/user/work/hg20812/project1

# Make STAR first pass output directory
mkdir -p $scratch/ALL_STAR_first_pass_files
# Make directory for STAR reports
mkdir -p ALL_star_reports

# Loop through each pair of trimmed reads
for i in ${scratch}/ALL_trimmed_reads/*_1.fastq.gz; do
    r1=$(basename $i)
    name=$(basename $i _1.fastq.gz)
    r2=${name}_2.fastq.gz

    # Run STAR with appropriate read group information
    STAR --runThreadN 12 \
         --genomeDir /user/work/hg20812/project/individual_test/STAR_mapping/Reference/StarRef1 \
         --sjdbGTFfile /user/work/hg20812/project/individual_test/STAR_mapping/Reference/Apis_mellifera.Amel_HAv3.1.56.gtf \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattrRGline ID:${name} LB:lib_${name} PL:ILLUMINA PU:Na SM:${name} \
         --outSAMunmapped Within \
         --readFilesCommand zcat \
         --readFilesIn ${scratch}/ALL_trimmed_reads/${r1} ${scratch}/ALL_trimmed_reads/${r2} \
         --outFileNamePrefix ${scratch}/ALL_STAR_first_pass_files/${name}_
done
