#!/bin/bash

#SBATCH --job-name=ALLHaploCaller
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --account=bisc028524

work_dir=/user/work/hg20812/project1
bam_files=/user/work/hg20812/project1/ALL_merged_bams  # Directory containing merged BAM files
ALL_called_variants_gvcf=/user/work/hg20812/project1/ALL_called_variants_gvcf

# Create the output directory if it doesn't exist
mkdir -p ${ALL_called_variants_gvcf}

# Load modules 
module load apps/gatk/4.1.9


# Call variants
for b in ${bam_files}/*.bam
do
    name=$(basename ${b} .bam) # Correct extension

    gatk HaplotypeCaller \
    -R /user/work/hg20812/project/individual_test/STAR_mapping/Reference/Apis_mellifera.Amel_HAv3.1.dna.toplevel.fa \
    --dont-use-soft-clipped-bases true \
    -I ${bam_files}/${name}.bam \
    -O ${ALL_called_variants_gvcf}/${name}.g.vcf \
    -ERC GVCF

    echo -e "${name}\t${ALL_called_variants_gvcf}/${name}.g.vcf" >> sample_map.txt

done
