# ============================================================================
# Build the presence/absence matrix for an UpSet plot from cached peptide lists
# Reads <CANCER>_tumor_peptides.txt or <CANCER>_normal_peptides.txt directly
# ============================================================================

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================================
# Configuration
# ============================================================================

# Cancer types
# cancer_types <- c("THCA", "BRCA", "CRC",
#                   "LUNG", "KIDNEY", "BLCA", "STAD", "LIHC")
cancer_types <- c("THCA", "BRCA", "CRC",
                  "LUNG", "KIDNEY", "BLCA", "STAD", "LIHC")

# Cached peptide-list directory
peptide_list_dir <- "/path/to/project/results/summary_V2/heatmap_peptides/data"

# Output directory
output_dir <- "/path/to/project/results/summary_V2/unpset_peptides"

# Peptide category to analyse:
# - "tumor_specific": peptides unique to tumour
# - "normal_specific": peptides unique to normal
# - "shared": peptides shared by tumour and normal
status_to_analyze <- "tumor_specific"

# ============================================================================
# Function definitions
# ============================================================================

#' Load a peptide list from file
#'
#' @param cancer_type cancer-type name
#' @param status peptide category ("tumor_specific", "normal_specific", "shared")
#' @param base_dir directory containing the peptide lists
#' @return character vector of peptides, or NULL if the file is missing
load_peptide_list <- function(cancer_type, status, base_dir) {
  # file_name <- paste0(cancer_type, "_", status, "_peptides.txt")
  # file_path <- file.path(base_dir, file_name)
  
  file_name <- paste0(cancer_type, "_", status, "_peptides.txt")
  file_path <- file.path(base_dir, cancer_type, file_name)
  
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }

  # Read the file, one peptide per line
  peptides <- readLines(file_path, warn = FALSE)

  # Drop blank lines and surrounding whitespace
  peptides <- trimws(peptides)
  peptides <- peptides[peptides != ""]

  cat(sprintf("%s (%s): loaded %d peptides\n",
              cancer_type, status, length(peptides)))

  return(peptides)
}


#' Build the presence/absence matrix
#'
#' @param peptide_list_by_cancer peptide lists named by cancer type
#' @return data.frame with peptides as rows and 0/1 columns per cancer type
generate_presence_matrix <- function(peptide_list_by_cancer) {

  # Convert to a long data frame
  epitope_long <- stack(peptide_list_by_cancer) %>%
    rename(Epitope = values, Cancer = ind) %>%
    mutate(Presence = 1)

  # Pivot to wide format (presence/absence matrix)
  epitope_matrix <- epitope_long %>%
    pivot_wider(names_from = Cancer,
                values_from = Presence,
                values_fill = 0)

  # Ensure all cancer columns are integer
  epitope_matrix <- epitope_matrix %>%
    mutate(across(-Epitope, ~ as.integer(.)))

  # Convert to data.frame
  epitope_matrix <- as.data.frame(epitope_matrix)

  return(epitope_matrix)
}


#' Print matrix statistics
#'
#' @param matrix presence/absence matrix
print_matrix_statistics <- function(matrix) {
  cat("\n" , rep("=", 80), "\n", sep = "")
  cat("Matrix statistics\n")
  cat(rep("=", 80), "\n", sep = "")

  cat(sprintf("Total peptides: %d\n", nrow(matrix)))
  cat(sprintf("Cancer types: %d\n", ncol(matrix) - 1))  # minus the Epitope column

  cat("\nPeptides per cancer type:\n")
  for (cancer in colnames(matrix)[-1]) {  # skip the Epitope column
    count <- sum(matrix[[cancer]])
    cat(sprintf("  %s: %d\n", cancer, count))
  }

  # Count peptides shared by all cancer types
  cancer_cols <- colnames(matrix)[-1]
  if (length(cancer_cols) > 1) {
    shared_all <- sum(rowSums(matrix[, cancer_cols, drop = FALSE]) == length(cancer_cols))
    cat(sprintf("\nPeptides shared by all cancer types: %d\n", shared_all))
  }

  cat(rep("=", 80), "\n", sep = "")
}

# ============================================================================
# Main
# ============================================================================

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("Build UpSet-plot matrix from peptide lists\n")
cat(rep("=", 80), "\n", sep = "")
cat(sprintf("Category: %s\n", toupper(status_to_analyze)))
cat(sprintf("Cancer types: %d\n", length(cancer_types)))
cat(sprintf("Peptide-list directory: %s\n", peptide_list_dir))
cat(rep("=", 80), "\n", sep = "")

# 1. Load peptide lists for all cancer types
cat("\nLoading peptide lists...\n")
cancer_peptide_list <- list()

for (cancer_type in cancer_types) {
  peptides <- load_peptide_list(cancer_type, status_to_analyze, peptide_list_dir)

  if (!is.null(peptides) && length(peptides) > 0) {
    cancer_peptide_list[[cancer_type]] <- peptides
  } else {
    warning(sprintf("Skipping %s: no peptide data", cancer_type))
  }
}

# Check that some data were loaded
if (length(cancer_peptide_list) == 0) {
  stop("Error: no peptide data were loaded.")
}

cat(sprintf("\nSuccessfully loaded data for %d cancer types\n", length(cancer_peptide_list)))

# 2. Build the presence/absence matrix
cat("\nBuilding presence/absence matrix...\n")
epitope_matrix <- generate_presence_matrix(cancer_peptide_list)

# 3. Print statistics
print_matrix_statistics(epitope_matrix)

# 4. Save results
cat("\nSaving results...\n")

# CSV file
csv_file <- file.path(output_dir,
                      sprintf("epitope_presence_matrix_%dcancers_%s.csv",
                              length(cancer_peptide_list),
                              status_to_analyze))
write.csv(epitope_matrix, csv_file, row.names = FALSE)
cat(sprintf("CSV file saved: %s\n", csv_file))

# RDS file (for fast loading in R)
rds_file <- file.path(output_dir,
                      sprintf("epitope_presence_matrix_%dcancers_%s.rds",
                              length(cancer_peptide_list),
                              status_to_analyze))
saveRDS(epitope_matrix, rds_file)
cat(sprintf("RDS file saved: %s\n", rds_file))

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("Done.\n")
cat(rep("=", 80), "\n", sep = "")

# ============================================================================
# Optional: generate the UpSet plot
# ============================================================================

# To generate the UpSet plot directly, uncomment the following:

# library(ComplexUpset)
#
# # Convert the matrix to a boolean format suitable for upset
# upset_data <- epitope_matrix %>%
#   mutate(across(-Epitope, ~ as.logical(.)))
#
# # Build the UpSet plot
# p <- upset(
#   upset_data,
#   colnames(upset_data)[-1],  # all cancer types except the Epitope column
#   name = 'Cancers',
#   width_ratio = 0.1,
#   min_size = 10,  # only show intersections with at least 10 peptides
#   sort_sets = 'descending',
#   sort_intersections = 'descending'
# )
#
# # Save the figure
# upset_file <- file.path(save_path,
#                         sprintf("upset_plot_%dcancers_%s.pdf",
#                                 length(cancer_peptide_list),
#                                 status_to_analyze))
# ggsave(upset_file, p, width = 12, height = 8)
# cat(sprintf("UpSet plot saved: %s\n", upset_file))