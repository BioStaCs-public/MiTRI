#!/bin/bash
# Map per-base coverage (bedgraph) onto fixed-size genome windows for each
# candidate microbe, separately for the responder (R) and non-responder (NR)
# groups. The window-averaged coverage is then compared in bedgraph_significance.R.
#
# Prerequisite: build fixed-size BED windows from each genome's .fai index, e.g.
#
#   INPUT_DIR="/path/to/project/data/reference"        # .fna genomes
#   OUTPUT_DIR="/path/to/project/data/reference_bed"   # .bed windows
#   WINDOW_SIZE=500                                     # requires samtools + bedtools
#   for GENOME in Fusobacterium_animalis Fusobacterium_nucleatum ... ; do
#       FNA_FILE="${INPUT_DIR}/${GENOME}.fna"
#       [ ! -f "${FNA_FILE}.fai" ] && samtools faidx "$FNA_FILE"
#       cut -f1,2 "${FNA_FILE}.fai" > "${FNA_FILE}.genome"
#       bedtools makewindows -g "${FNA_FILE}.genome" -w "$WINDOW_SIZE" > "${OUTPUT_DIR}/${GENOME}.bed"
#   done

# Microbial species to process
MICROORGANISMS=("Fusobacterium_nucleatum" "Fusobacterium_animalis" "Fusobacterium_canifelinum" "Fusobacterium_hwasookii" "Fusobacterium_massiliense" "Fusobacterium_polymorphum" "Fusobacterium_pseudoperiodonticum" "Fusobacterium_vincentii" "Fusobacterium_periodonticum")

# Input / output directories
BED_DIR="/path/to/project/data/reference_bed"
R_DIR="/path/to/project/results/special/CRC/coverage/R"
NR_DIR="/path/to/project/results/special/CRC/coverage/NR"
OUTPUT_DIR="/path/to/project/results/special/CRC/bedgraph_signif"
mkdir -p "$OUTPUT_DIR"

# Process each microbe
for MICROBE in "${MICROORGANISMS[@]}"
do
    # Responder group
    echo "Processing ${MICROBE} R data..."
    sort -k1,1 -k2,2n "${R_DIR}/${MICROBE}/coverage_R.bedgraph" > "${R_DIR}/${MICROBE}/coverage_R.sorted.bedgraph"
    bedtools map -a "${BED_DIR}/${MICROBE}.bed" -b "${R_DIR}/${MICROBE}/coverage_R.sorted.bedgraph" -c 4 -o mean > "${OUTPUT_DIR}/${MICROBE}_R_mapped.bed"

    # Non-responder group
    echo "Processing ${MICROBE} NR data..."
    sort -k1,1 -k2,2n "${NR_DIR}/${MICROBE}/coverage_NR.bedgraph" > "${NR_DIR}/${MICROBE}/coverage_NR.sorted.bedgraph"
    bedtools map -a "${BED_DIR}/${MICROBE}.bed" -b "${NR_DIR}/${MICROBE}/coverage_NR.sorted.bedgraph" -c 4 -o mean > "${OUTPUT_DIR}/${MICROBE}_NR_mapped.bed"
done

echo "All processes completed."

# The same mapping is used for the tumour-vs-normal cohorts by pointing R_DIR/NR_DIR
# at the Tumor/Normal coverage directories and the corresponding species list, e.g.
# BRCA: MICROORGANISMS=("Bifidobacterium_longum" "Cutibacterium_acnes")
# CRC : MICROORGANISMS=("Fusobacterium_nucleatum")
