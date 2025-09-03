#!/bin/bash

# Define the root directory
ROOT_DIR="/path/to/data/cancers"

# Define the failed samples record file
FAILED_SAMPLES_FILE="/path/to/data/cancers/failed_samples.txt"
> "$FAILED_SAMPLES_FILE"  # Clear the file

# Define the log file
LOG_FILE="/path/to/data/cancers/process_reads_ratio.log"
> "$LOG_FILE"  # Clear the log file

# Redirect stdout and stderr to the log file
exec > >(tee -a "$LOG_FILE") 2>&1

echo "Processing started at $(date)"

# Iterate through all cancer type directories
for cancer_dir in "$ROOT_DIR"/*; do
    # Check if it's a directory
    if [[ -d "$cancer_dir" ]]; then
        cancer_type=$(basename "$cancer_dir")
        echo "Processing cancer type: $cancer_type"
        
        # Define subdirectories (Tumor and Normal)
        sub_dirs=("Tumor" "Normal")
        
        # Iterate through Tumor and Normal subdirectories
        for sub_dir in "${sub_dirs[@]}"; do
            # Define the input, flagstat, and output directories
            input_dir="${cancer_dir}/${sub_dir}/bam"
            flagstat_dir="${cancer_dir}/${sub_dir}/flagstat"
            output_dir="${cancer_dir}/${sub_dir}/reads_ratio"

            # Create output directory if it does not exist
            mkdir -p "$output_dir"

            # Check if the BAM file directory exists
            if [[ ! -d "$input_dir" ]]; then
                echo "Warning: Directory $input_dir does not exist, skipping."
                continue
            fi

            # Iterate through all BAM files in the directory
            for bam_file in "$input_dir"/*_rna_output.pathseq.bam; do
                # Check if BAM files are present
                if [[ ! -f "$bam_file" ]]; then
                    echo "Warning: No BAM files found in $input_dir, skipping."
                    continue
                fi

                # Extract the sample name
                sample_name=$(basename "$bam_file" | sed 's/_rna_output.pathseq.bam//')

                # Define the path for the flagstat file
                flagstat_file="$flagstat_dir/${sample_name}_rna_hg38_align_sort.bam.flagstat.txt"

                # Check if flagstat file exists
                if [[ ! -f "$flagstat_file" ]]; then
                    echo "Warning: Missing flagstat file for $sample_name in $sub_dir, skipping."
                    continue
                fi

                # Retrieve total reads from the flagstat file
                total_reads=$(awk 'NR==1 {print $1}' "$flagstat_file")

                # Validate if total reads is a valid number
                if ! [[ "$total_reads" =~ ^[0-9]+$ ]]; then
                    echo "Error: Invalid total reads in $flagstat_file for $sample_name, skipping."
                    echo "$sample_name" >> "$FAILED_SAMPLES_FILE"
                    continue
                fi

                # Generate a temporary sorted BAM file
                temp_sorted_bam="${bam_file%.bam}.temp.sorted.bam"

                # Sort and index the BAM file
                samtools sort "$bam_file" -o "$temp_sorted_bam"
                samtools index "$temp_sorted_bam"

                # Calculate non-host reads
                non_host_reads=$(samtools idxstats "$temp_sorted_bam" | awk '$1 ~ /^N/ {sum += $3} END {print sum}')

                # Check if non-host reads are valid
                if [[ -z "$non_host_reads" ]]; then
                    echo "Error: Failed to calculate non-host reads for $sample_name, skipping."
                    rm -f "$temp_sorted_bam" "${temp_sorted_bam}.bai"
                    continue
                fi

                # Calculate the non-host read ratio
                non_host_ratio=$(echo "scale=6; $non_host_reads / $total_reads" | bc)

                # Output results to a file
                output_file="$output_dir/${sample_name}_non_host_ratio.txt"
                echo "Sample: $sample_name" > "$output_file"
                echo "Total reads: $total_reads" >> "$output_file"
                echo "Non-host reads: $non_host_reads" >> "$output_file"
                echo "Non-host ratio: $non_host_ratio" >> "$output_file"

                # Remove temporary files
                rm -f "$temp_sorted_bam" "${temp_sorted_bam}.bai"

                echo "Processed $sample_name in $sub_dir: Non-host ratio = $non_host_ratio"
            done
        done
    fi
done

echo "All cancer types processed. Results saved in respective directories."
echo "Processing finished at $(date)"
