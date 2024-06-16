import os
import subprocess

# Define the directory where the BAM files are located
bam_files_dir = "/user/work/hg20812/project1/bam_bqsr"

# Dictionary mapping individuals to their corresponding BAM files
bam_dict = {
    "211O10": ["SRR13889652.bam", "SRR13889705.bam"],
    "222G77": ["SRR13889719.bam", "SRR13889704.bam"],
    "231O48": ["SRR13889712.bam", "SRR13889696.bam"],
    "231O54": ["SRR13889718.bam", "SRR13889703.bam"],
    "241Y76": ["SRR13889717.bam", "SRR13889702.bam"],
    "242Y83": ["SRR13889716.bam", "SRR13889701.bam"],
    "251G2": ["SRR13889715.bam", "SRR13889699.bam"],
    "251G38": ["SRR13889714.bam", "SRR13889698.bam"],
    "321G50": ["SRR13889713.bam", "SRR13889697.bam"],
    "331O63": ["SRR13889710.bam", "SRR13889695.bam"],
    "341Y72": ["SRR13889709.bam", "SRR13889694.bam"],
    "351G40": ["SRR13889708.bam", "SRR13889693.bam"],
    "351G42": ["SRR13889707.bam", "SRR13889692.bam"],
    "361O81": ["SRR13889706.bam", "SRR13889687.bam"],
    "C1_1": ["SRR13889676.bam"],
    "C1_15": ["SRR13889686.bam", "SRR13889673.bam", "SRR13889659.bam"],
    "C1_2": ["SRR13889660.bam"],
    "C1_20": ["SRR13889683.bam"],
    "C1_7": ["SRR13889691.bam", "SRR13889688.bam", "SRR13889666.bam"],
    "C1_8": ["SRR13889674.bam"],
    "C2_1": ["SRR13889690.bam", "SRR13889682.bam", "SRR13889665.bam"],
    "C2_18": ["SRR13889671.bam", "SRR13889657.bam"],
    "C2_20": ["SRR13889675.bam", "SRR13889681.bam", "SRR13889663.bam"],
    "C2_21": ["SRR13889670.bam", "SRR13889656.bam"],
    "C2_5": ["SRR13889685.bam", "SRR13889672.bam", "SRR13889658.bam"],
    "C3_22": ["SRR13889664.bam", "SRR13889680.bam", "SRR13889662.bam"],
    "C3_26": ["SRR13889654.bam", "SRR13889679.bam", "SRR13889667.bam"],
    "C3_3": ["SRR13889689.bam"],
    "C3_32": ["SRR13889700.bam"],
    "C3_43": ["SRR13889678.bam", "SRR13889661.bam"],
    "C3_83": ["SRR13889669.bam", "SRR13889655.bam"],
    "C3_88": ["SRR13889684.bam", "SRR13889668.bam", "SRR13889653.bam"],
    "C3_97": ["SRR13889711.bam", "SRR13889677.bam"]
}

## Create the output directory
merged_bams_dir = "ALL_merged_bams"
os.makedirs(merged_bams_dir, exist_ok=True)

## Define the path to the Picard jar file
PICARD_JAR = "/sw/apps/picard/picard.jar"

## Create the function 'run_command' to run a shell command
def run_command(command):
    print(f"Running: {command}")
    subprocess.run(command, shell=True, check=True)

## Merge BAM files and add read groups
for individual, bams in bam_dict.items():
    merged_bam = os.path.join(merged_bams_dir, f"{individual}_merged.bam")
    input_bams = " ".join(os.path.join(bam_files_dir, bam) for bam in bams)
    
    # Merge BAM files
    merge_command = f"samtools merge {merged_bam} {input_bams}"
    run_command(merge_command)
    
    # Add or replace read groups
    output_bam = os.path.join(merged_bams_dir, f"{individual}_final.bam")
    add_rg_command = (
        f"java -jar {PICARD_JAR} AddOrReplaceReadGroups \\\n"
        f"       I={merged_bam} \\\n"
        f"       O={output_bam} \\\n"
        f"       RGID=4 \\\n"
        f"       RGLB=lib1 \\\n"
        f"       RGPL=ILLUMINA \\\n"
        f"       RGPU=unit1 \\\n"
        f"       RGSM={individual}"
    )
    run_command(add_rg_command)

    # Optionally remove the intermediate merged BAM file
    os.remove(merged_bam)

print("All BAM files have been processed.")
