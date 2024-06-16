#!/bin/bash

#SBATCH --job-name=PicardWorkflow
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --account=bisc028524

# Exit immediately if a command exits with a non-zero status
set -e  

# Define directories and reference file
bam_files=/user/work/hg20812/project1/ALL_STAR_first_pass_files
bam_proc=/user/work/hg20812/project1/bam_processing
bam_bqsr=/user/work/hg20812/project1/bam_bqsr
reference=/user/work/hg20812/project/individual_test/STAR_mapping/Reference/Apis_mellifera.Amel_HAv3.1.dna.toplevel.fa
vcf_file=/user/work/hg20812/project/individual_test/STAR_mapping/Reference/renamed_raw_QC.vcf

# Create necessary directories
mkdir -p ${bam_proc}
mkdir -p ${bam_bqsr}

# Load modules
module load apps/samtools/1.9
module load apps/picard/2.27.5
module load lang/java/1.8.0_201
module load apps/gatk/4.1.9

# Check if modules are loaded correctly
module list

# Path to picard.jar
PICARD_JAR=/sw/apps/picard/picard.jar

# Check if the picard.jar file exists
if [ ! -f "${PICARD_JAR}" ]; then
    echo "Picard JAR not found at ${PICARD_JAR}. Exiting."
    exit 1
fi

echo "Using Picard JAR at: ${PICARD_JAR}"

# Create FASTA index if it does not exist
if [ ! -f "${reference}.fai" ]; then
    echo "Creating FASTA index for ${reference}"
    samtools faidx ${reference} || { echo "Error creating FASTA index for ${reference}"; exit 1; }
fi

# Create VCF index if it does not exist
if [ ! -f "${vcf_file}.idx" ]; then
    echo "Creating VCF index for ${vcf_file}"
    gatk IndexFeatureFile -I ${vcf_file} || { echo "Error creating VCF index for ${vcf_file}"; exit 1; }
fi

# Set tags and process BAM files
for b in ${bam_files}/*.bam
do
    name=$(basename ${b} _Aligned.sortedByCoord.out.bam)

    # Set NM, MD, and UQ tags using Picard
    java -jar ${PICARD_JAR} SetNmMdAndUqTags \
    R=${reference} \
    I=${b} O=${bam_proc}/${name}_tags.bam || { echo "Error running SetNmMdAndUqTags on ${b}"; exit 1; }

    # Filter out non-primary and supplementary alignments
    samtools view -@ 8 -b -F 0x900 ${bam_proc}/${name}_tags.bam > ${bam_proc}/${name}_filt.bam || { echo "Error running Samtools on ${bam_proc}/${name}_tags.bam"; exit 1; }

    # Mark duplicates using Picard
    java -jar ${PICARD_JAR} MarkDuplicates \
    -I ${bam_proc}/${name}_filt.bam \
    -O ${bam_proc}/${name}_dedup.bam \
    -M ${bam_proc}/${name}_dup_metrics.txt \
    --REMOVE_DUPLICATES true || { echo "Error running MarkDuplicates on ${bam_proc}/${name}_filt.bam"; exit 1; }

    # Create sequence dictionary for reference genome (if it doesn't already exist)
    if [ ! -f ${reference%.*}.dict ]; then
        gatk CreateSequenceDictionary \
        -R ${reference} || { echo "Error creating sequence dictionary"; exit 1; }
    fi

    # Split Cigar Reads to enhance variant discovery
    gatk SplitNCigarReads \
    -R ${reference} \
    -I ${bam_proc}/${name}_dedup.bam \
    -O ${bam_proc}/${name}_cigar.bam || { echo "Error running SplitNCigarReads on ${bam_proc}/${name}_dedup.bam"; exit 1; }

    # Base recalibration
    gatk BaseRecalibrator \
    -R ${reference} \
    --known-sites ${vcf_file} \
    -I ${bam_proc}/${name}_cigar.bam \
    -O ${bam_bqsr}/${name}_recal_data.table || { echo "Error running BaseRecalibrator on ${bam_proc}/${name}_cigar.bam"; exit 1; }

    gatk ApplyBQSR \
    -R ${reference} \
    -I ${bam_proc}/${name}_cigar.bam \
    -bqsr-recal-file ${bam_bqsr}/${name}_recal_data.table \
    -O ${bam_bqsr}/${name}.bam || { echo "Error applying BQSR on ${bam_proc}/${name}_cigar.bam"; exit 1; }

    # Clean up intermediate files (optional, remove if you want to keep for debugging)
    rm ${bam_proc}/${name}_tags.bam ${bam_proc}/${name}_filt.bam ${bam_proc}/${name}_dedup.bam ${bam_proc}/${name}_cigar.bam

done
