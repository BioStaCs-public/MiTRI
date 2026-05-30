#!/bin/bash

# Record the start time
start_time=$(date +%s)

# Define cancer type and reference file
cancer_type="control"
cancer_type_reference="CRC_R_NR_revised_reference.csv"

# Set directory paths
# input_dir_gbm="/path/to/data/control/GBM_Normal/bam"
# input_dir_brain="/path/to/data/control/brain/bam"
input_dir_testicle="/path/to/data/control/testicle_6/bam"
output_dir_base="/path/to/results/$cancer_type"
reference_dir="/path/to/data/reference"
reference_list_dir="/path/to/data/CSVs_special/$cancer_type_reference"

# Read the Taxon column (first column) from the CSV file and remove the header
reference_list=$(cut -d',' -f1 "$reference_list_dir" | tail -n +2 | awk -v ref_dir="$reference_dir/" '{print ref_dir $1 ".fna"}')
export reference_list

# Define parallel processing function to handle GBM and brain coverage separately and merge them for overall brain coverage
process_reference_genome() {
    reference_genome=$1
    input_dir=$2
    output_dir_base=$3
    fastq_output_dir=$4
    cancer_type=$5
    
    reference_name=$(basename "$reference_genome" .fna)
    
    output_dir="$output_dir_base/$reference_name/$cancer_type"
    
    mkdir -p "$output_dir"
    
    {
        echo "Processing reference genome: $reference_name for $cancer_type"

        final_result_file="$output_dir/coverage_${cancer_type}.bedgraph"
        if [ -s "$final_result_file" ]; then
            echo "Final result $final_result_file already exists, skipping."
            return 0
        fi
         
        if [ ! -s "${reference_genome}.bwt" ]; then
           
            # Create BWA index
            echo "Indexing reference genome $reference_genome..."
            bwa index "$reference_genome"
        else
            echo "Reference genome index already exists."
        fi

        # Handle individual coverage for GBM and brain, and merge them for overall brain coverage
        if [ "$cancer_type" == "brain" ]; then
            for subtype in "GBM_Normal" "brain"; do
                subtype_fastq_dir="$output_dir_base/fastq_files/$subtype"
                subtype_output_dir="$output_dir_base/$reference_name/$subtype"
                mkdir -p "$subtype_output_dir"
                rm -rf "$subtype_output_dir"/*
                echo "Processing $subtype data with reference genome $reference_name..."
                for bam_file in "$subtype_fastq_dir"/*.fastq; do
                    file_name=$(basename "$bam_file" .fastq)
                    sam_file="$subtype_output_dir/${file_name%.bam}.sam"
                    new_bam_file="$subtype_output_dir/${file_name}_new.bam"
                    sorted_bam_file="$subtype_output_dir/${file_name}_sorted.bam"
                    bedgraph_file="$subtype_output_dir/${file_name}.bedgraph"

                    # Align the BAM file and save the result as a temporary SAM file
                    bwa mem -t 64 "$reference_genome" "$bam_file" > "$sam_file"

                    # Use Python to process the SAM file
                    python3 - <<END
import pysam
import re

# Get the file paths from Bash
sam_file = "$sam_file"
new_bam_file = "$new_bam_file"

# Open SAM file
samfile = pysam.AlignmentFile(sam_file, "r")

# Calculate the maximum match length
max_match_length = 0
for read in samfile:
    if not read.is_unmapped:
        # Get match length from CIGAR string (only consider M)
        cigar_string = read.cigarstring
        match_length = 0
        for match in re.finditer(r'(\d+)M', cigar_string):
            match_length += int(match.group(1))  # Accumulate the match lengths
        max_match_length = max(max_match_length, match_length)

print(f"Maximum match length calculated: {max_match_length}")
samfile.close()

# Filter and sort
samfile = pysam.AlignmentFile(sam_file, "r")
filtered_reads = []

for read in samfile:
    if read.is_unmapped:
        continue
    
    # Get match length from CIGAR string (only consider M)
    cigar_string = read.cigarstring
    match_length = 0
    for match in re.finditer(r'(\d+)M', cigar_string):
        match_length += int(match.group(1))

    # Filter based on max match length 95%
    if match_length >= max_match_length * 0.95:
        filtered_reads.append(read)

# Output filtered reads to new BAM file
if filtered_reads:
    try:
        with pysam.AlignmentFile(new_bam_file, "wb", header=samfile.header) as outfile:
            for read in filtered_reads:
                outfile.write(read)
        print(f"Successfully wrote filtered reads to {new_bam_file}")
    except Exception as e:
        print(f"Error writing BAM file: {e}")
else:
    print("No reads met the filter criteria, no BAM file created.")

samfile.close()
END
                    rm "$sam_file"
                    samtools sort -@ 16 "$new_bam_file" -o "$sorted_bam_file"
                    echo "Indexing $sorted_bam_file..."
                    samtools index "$sorted_bam_file"

                    echo "Calculating coverage for $sorted_bam_file..."
                    bedtools genomecov -ibam "$sorted_bam_file" -bga > "$bedgraph_file"
                done

                # Merge BAM files and calculate coverage for the subtype
                merged_bam="$subtype_output_dir/merged_${subtype}_sorted.bam"

                if [ ! -s "$merged_bam" ]; then
                    echo "Merging and processing $subtype group BAM files..."
                    samtools merge -@ 16 -f "$subtype_output_dir/merged_${subtype}.bam" "$subtype_output_dir"/*_sorted.bam
                    samtools sort -@ 16 "$subtype_output_dir/merged_${subtype}.bam" -o "$merged_bam"
                    samtools index -@ 16 "$merged_bam"
                    rm "$subtype_output_dir/merged_${subtype}.bam"
                else
                    echo "$merged_bam already exists, skipping."
                fi

                echo "Calculating coverage for merged $subtype group..."
                bedtools genomecov -ibam "$merged_bam" -bga > "$subtype_output_dir/coverage_${subtype}.bedgraph"
            done

            # Merge GBM and brain BAM files for overall brain coverage
            merged_bam="$output_dir/merged_brain_combined_sorted.bam"

            if [ ! -s "$merged_bam" ]; then
                echo "Merging and processing combined brain (GBM + brain) BAM files..."
                samtools merge -@ 16 -f "$output_dir/merged_brain_combined.bam" "$output_dir_base/$reference_name/GBM_Normal/merged_GBM_Normal_sorted.bam" "$output_dir_base/$reference_name/brain/merged_brain_sorted.bam"
                samtools sort -@ 16 "$output_dir/merged_brain_combined.bam" -o "$merged_bam"
                samtools index -@ 16 "$merged_bam"
                rm "$output_dir/merged_brain_combined.bam"
            else
                echo "$merged_bam already exists, skipping."
            fi

            echo "Calculating coverage for combined brain (GBM + brain)..."
            bedtools genomecov -ibam "$merged_bam" -bga > "$output_dir/coverage_brain_combined.bedgraph"
        
        elif [ "$cancer_type" == "testicle_6" ]; then
            # Process testicle sample separately
            echo "Processing testicle data with reference genome $reference_name..."
            for bam_file in "$fastq_output_dir"/*.fastq; do
                file_name=$(basename "$bam_file" .fastq)
                sam_file="$output_dir/${file_name%.bam}.sam"
                new_bam_file="$output_dir/${file_name}_new.bam"
                sorted_bam_file="$output_dir/${file_name}_sorted.bam"
                bedgraph_file="$output_dir/${file_name}.bedgraph"

                # Align and process the testicle data
                bwa mem -t 64 "$reference_genome" "$bam_file" > "$sam_file"
                python3 - <<END
# Python filtering as done in the previous section
END
                rm "$sam_file"
                samtools sort -@ 16 "$new_bam_file" -o "$sorted_bam_file"
                samtools index "$sorted_bam_file"

                # Calculate coverage for the testicle data
                bedtools genomecov -ibam "$sorted_bam_file" -bga > "$bedgraph_file"
            done

            # Merge BAM files for testicle samples and calculate coverage
            merged_bam="$output_dir/merged_testicle_sorted.bam"

            if [ ! -s "$merged_bam" ]; then
                echo "Merging and processing testicle group BAM files..."
                samtools merge -@ 16 -f "$output_dir/merged_testicle.bam" "$output_dir"/*_sorted.bam
                samtools sort -@ 16 "$output_dir/merged_testicle.bam" -o "$merged_bam"
                samtools index -@ 16 "$merged_bam"
                rm "$output_dir/merged_testicle.bam"
            else
                echo "$merged_bam already exists, skipping."
            fi

            echo "Calculating coverage for merged testicle group..."
            bedtools genomecov -ibam "$merged_bam" -bga > "$output_dir/coverage_testicle.bedgraph"
        fi
    rm "$output_dir"/*_sorted.bam "$output_dir"/*_sorted.bam.bai "$output_dir"/*_new.bam
    } >> "$output_dir_base/process_coverage_control.log" 2>&1
}

export -f process_reference_genome

# Use parallel to process brain (GBM and brain combined) and testicle reference genomes
echo "$reference_list" | parallel -j 10 process_reference_genome {} "$input_dir_brain" "$output_dir_base" "$output_dir_base/fastq_files/brain" "brain"
echo "$reference_list" | parallel -j 10 process_reference_genome {} "$input_dir_testicle" "$output_dir_base" "$output_dir_base/fastq_files/testicle_6" "testicle_6"

# Record the end time
end_time=$(date +%s)

# Calculate and display total elapsed time
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

echo "Processing complete." | tee -a "$output_dir_base/process_coverage_control.log"
echo "Total time: $hours hours, $minutes minutes, $seconds seconds." | tee -a "$output_dir_base/process_coverage_control.log"
