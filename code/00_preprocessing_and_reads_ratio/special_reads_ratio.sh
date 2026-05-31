#!/bin/bash
# Compute the non-host (microbial) reads ratio for a special cohort
# (e.g. immunopeptidomics / responder cohorts). For each sample, the
# microbial read count from the PathSeq BAM is divided by the total reads
# reported in the corresponding flagstat file.

# Input / output directories (edit to point at your data)
input_dir="/path/to/project/data/special/<COHORT>/bam"        # PathSeq BAM files
flagstat_dir="/path/to/project/data/special/<COHORT>/flagstat" # flagstat outputs
output_dir="/path/to/project/data/special/<COHORT>/reads_ratio" # results

# Create the output directory if it does not exist
mkdir -p "$output_dir"

# Iterate over all PathSeq BAM files in the input directory
for bam_file in "$input_dir"/*_rna_output.pathseq.bam; do
    # Sample name
    sample_name=$(basename "$bam_file" | sed 's/_rna_output.pathseq.bam//')

    # Matching flagstat file
    flagstat_file="$flagstat_dir/${sample_name}_rna_hg38_align_sort.bam.flagstat.txt"

    # Skip if the flagstat file is missing
    if [[ ! -f "$flagstat_file" ]]; then
        echo "Warning: Missing flagstat file for $sample_name, skipping."
        continue
    fi

    # Total read count from the flagstat file
    total_reads=$(awk 'NR==6 {print $1}' "$flagstat_file")

    # Validate that the total read count is a number
    if ! [[ "$total_reads" =~ ^[0-9]+$ ]]; then
        echo "Error: Invalid total reads in $flagstat_file for $sample_name, skipping."
        continue
    fi

    # Temporary sorted BAM
    temp_sorted_bam="${bam_file%.bam}.temp.sorted.bam"

    # Sort and index the BAM
    samtools sort "$bam_file" -o "$temp_sorted_bam"
    samtools index "$temp_sorted_bam"

    # Count microbial reads (PathSeq contigs are prefixed with 'N')
    non_host_reads=$(samtools idxstats "$temp_sorted_bam" | awk '$1 ~ /^N/ {sum += $3} END {print sum}')

    # Validate microbial read count
    if [[ -z "$non_host_reads" ]]; then
        echo "Error: Failed to calculate non-host reads for $temp_sorted_bam, skipping."
        rm -f "$temp_sorted_bam" "${temp_sorted_bam}.bai"
        continue
    fi

    # Non-host (microbial) ratio
    non_host_ratio=$(echo "scale=6; $non_host_reads / $total_reads" | bc)

    # Write per-sample result
    output_file="$output_dir/${sample_name}_non_host_ratio.txt"
    echo "Sample: $sample_name" > "$output_file"
    echo "Total reads: $total_reads" >> "$output_file"
    echo "Non-host reads: $non_host_reads" >> "$output_file"
    echo "Non-host ratio: $non_host_ratio" >> "$output_file"

    # Clean up temporary files
    rm -f "$temp_sorted_bam" "${temp_sorted_bam}.bai"

    echo "Processed $sample_name: Non-host ratio = $non_host_ratio"
done

echo "All files processed. Results saved in $output_dir."
