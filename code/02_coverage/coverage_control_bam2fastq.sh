#!/bin/bash

# Record start time
start_time=$(date +%s)

# Define cancer type
cancer_type="control"

# Set directory paths
# input_dir_gbm="/path/to/data/control/GBM_Normal/bam"
# input_dir_brain="/path/to/data/control/brain/bam"
input_dir_testicle="/path/to/data/control/testicle_6/bam"

output_dir_base="/path/to/results/$cancer_type"

# Directly use a single cancer type value
cancer_type="testicle_6"
input_dir="$input_dir_testicle"
file_prefix="S"

fastq_output_dir="$output_dir_base/fastq_files/$cancer_type"
mkdir -p "$fastq_output_dir"

echo "Extracting FASTQ files for $cancer_type group..." >> "$output_dir_base/process_coverage_control.log"

# Iterate through specified BAM files and convert them to FASTQ
for bam_file in "$input_dir"/"$file_prefix"*.bam; do
    file_name=$(basename "$bam_file")
    fastq_file="$fastq_output_dir/${file_name%.bam}.fastq"
    
    # Check if the FASTQ file already exists
    if [ ! -s "$fastq_file" ]; then
        echo "Converting $bam_file to FASTQ..." >> "$output_dir_base/process_coverage_control.log"
        samtools fastq "$bam_file" > "$fastq_file" 2>> "$output_dir_base/process_coverage_control.log"
    else
        echo "FASTQ file $fastq_file already exists, skipping extraction." >> "$output_dir_base/process_coverage_control.log"
    fi
done

