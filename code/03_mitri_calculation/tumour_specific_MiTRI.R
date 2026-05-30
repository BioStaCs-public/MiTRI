setwd("/path/to/project/")
# Load required packages
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(data.table)
library(parallel)
library(openxlsx)  # Excel I/O
library(readxl)

# Tumour-specific MiTRI: restrict the activity signal to genomic positions
# covered in the tumour but NOT in the matched normal, then normalise by the
# non-host read ratio and the length of the retained (tumour-restricted) region.

# Cancer types to process (parameterised; the full panel is listed below)
cancer_types <- c("CRC", "BRCA", "CESC", "LIHC", "OSCC", "LUNG", "PAAD", "STAD", "THCA")
for (cancer_type in cancer_types) {
  cat("cancer_type", cancer_type, ":\n")
  reflist <- paste0(cancer_type, "_V4.csv")

  output_base <- file.path("/path/to/project/results/cancers", cancer_type, "03.coverage")

  output_dir <- file.path("/path/to/project/results/cancers", cancer_type, "04.activity")
  dir.create(output_dir, recursive = TRUE)

  flagstat_base = file.path("data/cancers", cancer_type, "Tumor/flagstat/")
  non_host_ratio_base <- file.path("data/cancers", cancer_type, "Tumor/reads_ratio/")

  reference_base <- file.path("/path/to/project/data/reference")
  Reference_length_dir <- file.path("/path/to/project/data/CSVs", reflist)
  reference_lengths <- fread(Reference_length_dir, colClasses = list(numeric = "Flag"))

  # Record the start time
  start_time <- Sys.time()
  cat("Script start time:", format(start_time), "\n")

  generate_boolean_vector <- function(coverage) {
    seq_names <- unique(seqnames(coverage))
    genome_length <- sum(sapply(split(coverage, seqnames(coverage)), function(gr) max(end(gr))))
    boolean_vector <- rep(0L, genome_length)
    offset <- 0

    for (chr in seq_names) {
      chr_coverage <- coverage[seqnames(coverage) == chr]
      chr_length <- max(end(chr_coverage))
      # Keep only ranges with score > 0
      valid_ranges <- chr_coverage[chr_coverage$score > 0]

      # Vectorised assignment
      starts <- start(valid_ranges) + offset
      ends <- end(valid_ranges) + offset

      for (i in seq_along(starts)) {
        boolean_vector[starts[i]:ends[i]] <- 1L
      }

      # Advance the offset
      offset <- offset + chr_length

    }
    return(boolean_vector)
  }

  # Remove positions where tumour and normal coverage overlap
  remove_overlap_regions <- function(t_boolean, other_boolean) {
    if (length(t_boolean) != length(other_boolean)) {
      stop("t_boolean and other_boolean have different lengths")
    }
    # Initialise the result vector
    result <- integer(length(t_boolean))

    # Classify each position
    overlap <- (t_boolean == 1) & (other_boolean == 1)   # covered in both
    only_red <- (t_boolean == 1) & (other_boolean == 0)  # tumour only
    only_blue <- (t_boolean == 0) & (other_boolean == 1) # normal only
    neither <- (t_boolean == 0) & (other_boolean == 0)   # covered in neither

    # Encode into the result vector
    result[overlap] <- -1
    result[only_red] <- 1
    result[only_blue] <- 2
    result[neither] <- 0

    return(result)
  }
  calculate_remaining_coverage <- function(coverage, boolean_vector, non_host_ratio, activity_region_length) {
    # Per-base coverage weighted by "score"
    coverage_values <- coverage(coverage, weight = "score")
    total_coverage <- 0
    global_pos_index <- 1

    # Iterate over per-chromosome coverage
    for (i in seq_along(coverage_values)) {
      coverage_rle <- coverage_values[[i]]  # per-chromosome coverage RLE
      run_lengths <- runLength(coverage_rle)  # RLE run lengths
      run_values <- runValue(coverage_rle)    # RLE run values
      chr_length <- sum(run_lengths)          # chromosome length

      # Global position indices for this chromosome
      start_positions <- global_pos_index + cumsum(c(0, run_lengths[-length(run_lengths)]))
      end_positions <- pmin(start_positions + run_lengths - 1, length(boolean_vector))

      # Sum coverage over active (tumour-restricted) positions
      for (j in seq_along(run_lengths)) {
        boolean_slice <- boolean_vector[start_positions[j]:end_positions[j]]
        active_bases <- sum(boolean_slice == 1)  # active bases in this run

        if (active_bases > 0) {
          total_coverage <- total_coverage + run_values[j] * active_bases
        }
      }

      # Advance the global position index
      global_pos_index <- global_pos_index + chr_length
    }

    # Coverage ratio
    ratio <- (total_coverage / non_host_ratio) / activity_region_length
    return(list(ratio = ratio))
  }

  # Read the non-host (microbial) read ratio for a sample
  read_non_host_ratio <- function(sample_name) {
    # Result file path
    reads_ratio_file <- paste0(non_host_ratio_base, sample_name, "_non_host_ratio.txt")

    # Check the file exists
    if (file.exists(reads_ratio_file)) {
      lines <- readLines(reads_ratio_file)

      # Find and extract the "Non-host ratio" value
      ratio_line <- grep("Non-host ratio:", lines, value = TRUE)

      if (length(ratio_line) == 1) {
        reads_ratio <- as.numeric(sub("Non-host ratio: ", "", ratio_line))
        return(reads_ratio)
      } else {
        stop(paste("Failed to find Non-host ratio in file:", reads_ratio_file))
      }
    } else {
      stop(paste("No reads_ratio file found for sample:", sample_name))
    }
  }

  results_list <- list()

  reference_dir <- file.path(reference_base)

  tumor_dir_base <- file.path(output_base, "Tumor")
  normal_dir_base <- file.path(output_base, "Normal")

  fna_files <- list.files(reference_dir, pattern = "\\.fna$", full.names = TRUE)
  reference_names <- basename(fna_files)
  reference_names <- gsub("\\.fna$", "", reference_names)

  valid_reference_names <- reference_names[reference_names %in% reference_lengths$Taxon]

  results_dt <- data.table(
    Reference = character(),
    Activity = numeric()
  )

  results_dt <- rbindlist(mclapply(valid_reference_names, function(reference_name) {
    tryCatch({
      tumor_dir <- file.path(tumor_dir_base, reference_name)
      normal_dir <- file.path(normal_dir_base, reference_name)

      # Import tumour and matched-normal coverage files
      t_files <- list.files(tumor_dir, pattern = ".*_rna_.*\\.bedgraph$", full.names = TRUE)
      t_coverages <- lapply(t_files, import)

      n_file <- file.path(normal_dir, "coverage_Normal.bedgraph")
      n_coverage <- import(n_file)

      t_boolean_vectors <- lapply(t_coverages, generate_boolean_vector)
      n_boolean_vectors <- generate_boolean_vector(n_coverage)

      # reference_length <- reference_lengths[Taxon == reference_name, Reference_length] / 10^3

      calculate_coverage_for_control <- function(control_vector) {
        lapply(seq_along(t_boolean_vectors), function(i) {
          file_name <- basename(t_files[i])
          sample_name <- strsplit(file_name, "_rna")[[1]][1]

          # total_library_size <- read_total_library_value(sample_name) / 10^6
          non_host_ratio <- read_non_host_ratio(sample_name)
          cleaned_vector <- remove_overlap_regions(t_boolean_vectors[[i]], control_vector)

          activity_region_length <- sum(cleaned_vector != 0)
          remaining_coverage <- calculate_remaining_coverage(t_coverages[[i]], cleaned_vector, non_host_ratio, activity_region_length)

          list(sample_name = sample_name, ratio = remaining_coverage$ratio)

        })
      }

      # Compute coverage against the matched-normal control
      normal_results <- calculate_coverage_for_control(n_boolean_vectors)

      # Extract the tumour-restricted ratios
      Signal_normal_ratio_values <- sapply(normal_results, function(x) x$ratio)
      sample_names <- sapply(normal_results, function(x) x$sample_name)

      # Return the result
      data.table(
        Reference = reference_name,
        Sample_Name = sample_names,
        Activity = Signal_normal_ratio_values
      )

    }, error = function(e) {
      # Catch and log errors, returning a flagged result
      message(paste("Error processing reference:", reference_name, "\n", e))
      return(data.table(
        Reference = reference_name,
        Sample_Name = sample_names,
        Activity = NA
      ))
    })
  }, mc.cores = 5), use.names = TRUE, fill = TRUE)

  results_list[] <- results_dt

  # Save the results to an Excel file
  output_xlsx_file <- file.path(output_dir,  paste0(cancer_type,"_tumor_specific_MiTRI.xlsx"))
  write.xlsx(results_dt, file = output_xlsx_file)

  cat("Results saved to:", output_xlsx_file, "\n")
}
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "secs")
elapsed_time_formatted <- sprintf("%02d:%02d:%02d", as.integer(elapsed_time) %/% 3600, (as.integer(elapsed_time) %% 3600) %/% 60, as.integer(elapsed_time) %% 60)

cat("Total run time:", elapsed_time_formatted, "\n")
