#!/bin/bash

# Record the start time
start_time=$(date +%s)

# Set base directory paths
base_dir="/path/to/data"
input_dir_base="$base_dir/data/cancers_V1"
output_dir_base="$base_dir/results/cancers_V1"
reference_dir="$base_dir/data/reference"

# List of cancer types, more types can be added as needed
cancer_types=("CRC")

# Loop through each cancer type
for cancer_type in "${cancer_types[@]}"; do
    echo "Processing cancer type: $cancer_type"

    # Set paths for this cancer type
    cancer_type_reference="CRC_periodontion.csv"
    reference_list_dir="/path/to/data/CSVs_V1/$cancer_type_reference"

    input_dir="$input_dir_base/$cancer_type"
    output_dir_input="$output_dir_base/$cancer_type/03.coverage"
    fastq_output_dir="$output_dir_base/$cancer_type/03.coverage/fastq_files"
    mkdir -p "$fastq_output_dir"

    # Read the Taxon column (first column) from the CSV file and remove the header
    reference_list=$(cut -d',' -f1 "$reference_list_dir" | tail -n +2 | awk -v ref_dir="$reference_dir/" '{print ref_dir $1 ".fna"}')
    export reference_list

    # Define the parallel processing function
    process_reference_genome() {
        reference_genome=$1
        input_dir=$2
        output_dir_base=$3
        fastq_output_dir=$4
        data_type=$5

        reference_name=$(basename "$reference_genome" .fna)
        output_dir="$output_dir_base/$data_type/$reference_name"
        mkdir -p "$output_dir"
        
        {
            echo "Processing reference genome: $reference_name for $data_type"

            # Generate reference genome index
            if [ ! -s "${reference_genome}.bwt" ]; then
                echo "Indexing reference genome $reference_genome..."
                bwa index "$reference_genome"
            else
                echo "Reference genome $reference_genome already indexed, skipping."
            fi
            # Process BAM files
            echo "Processing non-control group data for $reference_name..."
            
            for bam_file in "$input_dir"/*.bam; do
                file_name=$(basename "$bam_file")
                fastq_file="$fastq_output_dir/${file_name%.bam}.fastq"
                sam_file="$output_dir/${file_name%.bam}.sam"
                new_bam_file="$output_dir/${file_name%.bam}_new.bam"
                sorted_bam_file="$output_dir/${file_name%.bam}_sorted.bam"
                bedgraph_file="$output_dir/${file_name%.bam}.bedgraph"

                # Perform alignment and save the result as a temporary SAM file
                bwa mem -t 32 "$reference_genome" "$fastq_file" > "$sam_file"
                
                # Use Python to process the SAM file
                python3 - <<END
import pysam
import re

# Get file paths passed from Bash
sam_file = "$sam_file"
new_bam_file = "$new_bam_file"

# Open SAM file
samfile = pysam.AlignmentFile(sam_file, "r")

# Calculate the maximum match length
max_match_length = 0
for read in samfile:
    if not read.is_unmapped:
        cigar_string = read.cigarstring
        match_length = 0
        for match in re.finditer(r'(\d+)M', cigar_string):
            match_length += int(match.group(1))
        max_match_length = max(max_match_length, match_length)

print(f"Maximum match length calculated: {max_match_length}")
samfile.close()

# Filter and sort
samfile = pysam.AlignmentFile(sam_file, "r")
filtered_reads = []

for read in samfile:
    if read.is_unmapped:
        continue
    cigar_string = read.cigarstring
    match_length = 0
    for match in re.finditer(r'(\d+)M', cigar_string):
        match_length += int(match.group(1))

    if match_length >= max_match_length * 0.95:
        filtered_reads.append(read)

if filtered_reads:
    with pysam.AlignmentFile(new_bam_file, "wb", header=samfile.header) as outfile:
        for read in filtered_reads:
            outfile.write(read)
    print(f"Successfully wrote filtered reads to {new_bam_file}")
else:
    print("No reads met the filter criteria, no BAM file created.")
samfile.close()
END
                rm "$sam_file"
                samtools sort -@ 16 "$new_bam_file" -o "$sorted_bam_file"
                samtools index "$sorted_bam_file"
                bedtools genomecov -ibam "$sorted_bam_file" -bga > "$bedgraph_file"
            done

            # Merge BAM files
            merged_bam="$output_dir/merged_sorted.bam"
            if [ ! -s "$merged_bam" ]; then
                samtools merge -@ 16 -f "$output_dir/merged.bam" "$output_dir"/*_sorted.bam
                samtools sort -@ 16 "$output_dir/merged.bam" -o "$merged_bam"
                samtools index -@ 16 "$merged_bam"
                rm "$output_dir/merged.bam"
            else
                echo "$merged_bam already exists, skipping."
            fi

            bedtools genomecov -ibam "$merged_bam" -bga > "$output_dir/coverage_${data_type}.bedgraph"
            
            rm "$output_dir"/*_sorted.bam "$output_dir"/*_sorted.bam.bai  "$output_dir"/*_new.bam
        } >> "$output_dir_base/process_coverage_main.log" 2>&1
    }

    export -f process_reference_genome

    # Loop to process Normal and Tumor data types
    for data_type in Normal; do
        input_dir="$input_dir_base/$cancer_type/$data_type/bam"
        output_dir="$output_dir_base/$cancer_type/03.coverage/$data_type"
        mkdir -p "$output_dir"

        fastq_output_dir="$output_dir/fastq_files"
        mkdir -p "$fastq_output_dir"

        # Convert BAM files to FASTQ
        for bam_file in "$input_dir"/*.bam; do
            file_name=$(basename "$bam_file")
            fastq_file="$fastq_output_dir/${file_name%.bam}.fastq"
            
            if [ ! -s "$fastq_file" ]; then
                echo "Converting $bam_file to FASTQ..." >> "$output_dir/process_coverage_main.log"
                samtools fastq "$bam_file" > "$fastq_file" 2>> "$output_dir/process_coverage_main.log"
            else
                echo "FASTQ file $fastq_file already exists, skipping extraction." >> "$output_dir/process_coverage_main.log"
            fi
        done

        # Use parallel processing to handle reference genomes
        echo "$reference_list" | parallel -j 5 process_reference_genome {} "$input_dir" "$output_dir_input" "$fastq_output_dir" "$data_type"
    done
done

# Record the end time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

echo "Processing complete."
echo "Total time: $hours hours, $minutes minutes, $seconds seconds."
