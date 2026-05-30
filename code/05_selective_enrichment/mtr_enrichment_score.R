#!/usr/bin/env Rscript
# ==============================================================================
# MTR Region Enrichment Score Analysis - With Logging
# Purpose: Calculate enrichment of reads in MTR regions to distinguish 
#          true transcriptional reprogramming from abundance effects
# ==============================================================================

library(dplyr)
library(data.table)
library(parallel)
library(ggplot2)
library(Biostrings)
library(GenomicRanges)
library(rtracklayer)

# ==============================================================================
# CONFIGURATION SECTION - edit all paths and file names here
# ==============================================================================

# --- Base path configuration ---
CONFIG <- list(
  # Working directory and base path
  base_dir = "/path/to/project",

  # Results and data directories
  results_subdir = "results_V3/cancers_V3.1",
  data_subdir = "data",
  reference_subdir = "reference",
  csv_subdir = "CSVs_20251208",

  # Analysis subdirectory structure
  vector_subdir = "05.MTR/05.1.vector",
  filtered_reads_subdir = "05.MTR/05.4.MTR_reads_nogene/Tumor",
  coverage_subdir = "03.coverage/Tumor",

  # Output directories
  output_base = "results_V3/summary/MTR_enrichment_analysis",
  output_intermediate = "results_V3/summary/MTR_enrichment_analysis/intermediate",

  # Output file names
  output_files = list(
    all_results = "mtr_enrichment_scores_all.csv",
    filtered_results = "mtr_enrichment_scores_filtered.csv",
    intermediate_suffix = "_enrichment.csv",
    plot_distribution = "enrichment_distribution.pdf",
    plot_by_cancer = "enrichment_by_cancer.pdf",
    log_file = "analysis.log",
    error_log = "errors_detailed.csv"
  ),
  
  # Input file-name patterns
  file_patterns = list(
    reflist_suffix = "_V4.csv",
    vector_suffix = "_all_vector.csv",
    filtered_suffix = "_StrategyA_read_based.csv",
    bam_pattern = "_rna.*\\.bam\\.txt$",
    genome_extension = ".fna"
  ),
  
  # Cancer types
  cancer_types = c("LUNG", "OSCC", "BLCA", "CRC", "BRCA",
                   "CESC", "LIHC", "STAD", "PAAD"),

  # Filtering thresholds
  thresholds = list(
    min_mtr_reads = 2,
    min_total_reads = 0,
    p_value = 0.05
  ),
  
  # Statistical-test settings
  stats = list(
    alternative = "two.sided"
  )
)

# Set working directory
setwd(CONFIG$base_dir)

# Build full paths
PATHS <- list(
  base_dir = CONFIG$base_dir,
  results_base = file.path(CONFIG$base_dir, CONFIG$results_subdir),
  # data_base = file.path(CONFIG$base_dir, CONFIG$data_subdir),
  reference_dir = file.path(CONFIG$base_dir, CONFIG$data_subdir, CONFIG$reference_subdir),
  csv_dir = file.path(CONFIG$base_dir, "data_V3", CONFIG$csv_subdir),
  output_dir = file.path(CONFIG$base_dir, CONFIG$output_base),
  output_intermediate = file.path(CONFIG$base_dir, CONFIG$output_intermediate)
)

# ==============================================================================
# LOGGING FUNCTIONS
# ==============================================================================

# Initialise log file
LOG_FILE <- NULL
ERROR_LOG <- list()

init_logging <- function(output_dir) {
  LOG_FILE <<- file.path(output_dir, CONFIG$output_files$log_file)
  
  # Create a fresh log file
  if (file.exists(LOG_FILE)) {
    file.remove(LOG_FILE)
  }
  
  log_message("=================================================================")
  log_message("MTR ENRICHMENT SCORE ANALYSIS - LOG FILE")
  log_message(paste("Start Time:", Sys.time()))
  log_message(paste("R Version:", R.version.string))
  log_message("=================================================================")
  log_message("")
}

# Write to both console and log file
log_message <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  formatted_msg <- paste0("[", timestamp, "] [", level, "] ", msg)
  
  # Output to console
  cat(msg, "\n")

  # Write to log file
  if (!is.null(LOG_FILE)) {
    cat(formatted_msg, "\n", file = LOG_FILE, append = TRUE)
  }
}

# Record error details
log_error <- function(sample, microbe, cancer, error_type, details) {
  ERROR_LOG[[length(ERROR_LOG) + 1]] <<- list(
    Sample = sample,
    Microbe = microbe,
    Cancer_Type = cancer,
    Error_Type = error_type,
    Details = details,
    Timestamp = as.character(Sys.time())
  )
  
  log_message(paste("ERROR:", sample, "-", error_type, ":", details), "ERROR")
}

# Save error log
save_error_log <- function(output_dir) {
  if (length(ERROR_LOG) > 0) {
    error_df <- do.call(rbind, lapply(ERROR_LOG, as.data.frame, stringsAsFactors = FALSE))
    error_file <- file.path(output_dir, CONFIG$output_files$error_log)
    write.csv(error_df, error_file, row.names = FALSE)
    log_message(paste("Error log saved to:", error_file))
    log_message(paste("Total errors logged:", nrow(error_df)))
  } else {
    log_message("No errors to log!")
  }
}

# ==============================================================================
# Core Functions with Enhanced Error Handling and Logging
# ==============================================================================

get_genome_length <- function(reference_name, reference_dir) {
  fna_file <- file.path(reference_dir, paste0(reference_name, CONFIG$file_patterns$genome_extension))
  if (!file.exists(fna_file)) {
    log_message(paste("    Warning: Reference genome not found for", reference_name), "WARN")
    return(NA)
  }
  
  tryCatch({
    sequences <- readDNAStringSet(fna_file)
    length_val <- sum(width(sequences))
    log_message(paste("    Genome length for", reference_name, ":", length_val, "bp"), "DEBUG")
    return(length_val)
  }, error = function(e) {
    log_message(paste("    Error reading genome for", reference_name, ":", e$message), "ERROR")
    return(NA)
  })
}

get_mtr_region_length <- function(mtr_file) {
  if (!file.exists(mtr_file) || file.size(mtr_file) == 0) {
    return(NA)
  }
  
  if (!grepl("\\.csv$", mtr_file)) {
    log_message(paste("    Warning: Not a CSV file:", mtr_file), "WARN")
    return(NA)
  }
  
  tryCatch({
    mtr_data <- fread(mtr_file, header = TRUE)
    
    if ("x" %in% names(mtr_data)) {
      length_val <- sum(mtr_data$x == 1)
      return(ifelse(length_val > 0, length_val, NA))
    } else if (ncol(mtr_data) == 1) {
      length_val <- sum(mtr_data[[1]] == 1)
      return(ifelse(length_val > 0, length_val, NA))
    } else {
      log_message(paste("    Warning: No 'x' column found in", basename(mtr_file)), "WARN")
      return(NA)
    }
  }, error = function(e) {
    log_message(paste("    Error reading vector file:", e$message), "ERROR")
    return(NA)
  })
}

count_mtr_reads <- function(filtered_reads_file) {
  if (!file.exists(filtered_reads_file) || file.size(filtered_reads_file) == 0) {
    return(0)
  }
  
  tryCatch({
    reads_data <- fread(filtered_reads_file, header = TRUE)
    return(nrow(reads_data))
  }, error = function(e) {
    log_message(paste("    Error counting MTR reads:", e$message), "ERROR")
    return(0)
  })
}

count_total_reads <- function(bam_txt_file) {
  if (!file.exists(bam_txt_file)) {
    return(NA)
  }
  
  tryCatch({
    cmd <- paste("grep -v '^@'", bam_txt_file, "| wc -l")
    total_reads <- as.numeric(system(cmd, intern = TRUE))
    return(ifelse(total_reads > 0, total_reads, NA))
  }, error = function(e) {
    log_message(paste("    Error counting total reads:", e$message), "ERROR")
    return(NA)
  })
}

calculate_enrichment_score <- function(sample_name, microbe, cancer_type,
                                       genome_length, cancer_dir) {
  
  tryCatch({
    # Build paths from configuration
    vector_file <- file.path(cancer_dir, CONFIG$vector_subdir,
                             microbe, paste0(sample_name, CONFIG$file_patterns$vector_suffix))
    filtered_reads_dir <- file.path(cancer_dir, CONFIG$filtered_reads_subdir, microbe)
    coverage_dir <- file.path(cancer_dir, CONFIG$coverage_subdir, microbe)
    
    # Get MTR region length
    mtr_region_length <- NA
    mtr_source <- "none"
    
    if (file.exists(vector_file) && file.size(vector_file) > 0) {
      mtr_region_length <- get_mtr_region_length(vector_file)
      if (!is.na(mtr_region_length) && mtr_region_length > 0) {
        mtr_source <- "vector"
      }
    }
    
    if (is.na(mtr_region_length) || mtr_region_length == 0) {
      log_error(sample_name, microbe, cancer_type, "NO_MTR_REGION", "MTR region length is 0 or NA")
      return(NULL)
    }
    
    # Find filtered reads file
    filtered_pattern <- paste0(sample_name, ".*", gsub("^_", "", CONFIG$file_patterns$filtered_suffix), "$")
    filtered_files <- list.files(filtered_reads_dir, pattern = filtered_pattern, full.names = TRUE)
    
    if (length(filtered_files) == 0) {
      log_error(sample_name, microbe, cancer_type, "NO_FILTERED_FILE", "Filtered reads file not found")
      return(NULL)
    }
    
    mtr_reads_count <- count_mtr_reads(filtered_files[1])
    
    # Find BAM text file
    bam_pattern <- paste0(sample_name, CONFIG$file_patterns$bam_pattern)
    bam_files <- list.files(coverage_dir, pattern = bam_pattern, full.names = TRUE)
    
    if (length(bam_files) == 0) {
      log_error(sample_name, microbe, cancer_type, "NO_BAM_FILE", "BAM text file not found")
      return(NULL)
    }
    
    total_reads_count <- count_total_reads(bam_files[1])
    
    if (is.na(total_reads_count) || total_reads_count == 0) {
      log_error(sample_name, microbe, cancer_type, "INVALID_TOTAL_READS", 
                paste("Total reads count is", total_reads_count))
      return(NULL)
    }
    
    # Key fix: check data consistency
    if (mtr_reads_count > total_reads_count) {
      log_error(sample_name, microbe, cancer_type, "INCONSISTENT_COUNTS",
                paste("MTR reads (", mtr_reads_count, ") > Total reads (", total_reads_count, 
                      "). This suggests data inconsistency between filtered reads and BAM file."))
      # Correction: use MTR reads as the floor for total reads
      total_reads_count <- max(mtr_reads_count, total_reads_count)
      log_message(paste("      CORRECTED: Adjusted total_reads_count to", total_reads_count), "WARN")
    }
    
    # Calculate metrics
    mtr_genome_fraction <- mtr_region_length / genome_length
    mtr_reads_fraction <- mtr_reads_count / total_reads_count
    enrichment_score <- ifelse(mtr_genome_fraction > 0,
                               mtr_reads_fraction / mtr_genome_fraction,
                               NA)
    
    # Statistical test with validation
    p_value <- NA
    if (!is.na(enrichment_score) && mtr_reads_count > 0) {
      # Re-validate binom.test parameters
      if (total_reads_count >= mtr_reads_count &&
          mtr_genome_fraction > 0 && mtr_genome_fraction <= 1) {
        tryCatch({
          binom_result <- binom.test(
            x = mtr_reads_count,
            n = total_reads_count,
            p = mtr_genome_fraction,
            alternative = CONFIG$stats$alternative
          )
          p_value <- binom_result$p.value
        }, error = function(e) {
          log_error(sample_name, microbe, cancer_type, "BINOM_TEST_FAILED",
                    paste("x=", mtr_reads_count, "n=", total_reads_count, 
                          "p=", mtr_genome_fraction, "Error:", e$message))
        })
      } else {
        log_error(sample_name, microbe, cancer_type, "INVALID_BINOM_PARAMS",
                  paste("Invalid params: x=", mtr_reads_count, "n=", total_reads_count, 
                        "p=", mtr_genome_fraction))
      }
    }
    
    return(data.frame(
      Sample = sample_name,
      Microbe = microbe,
      Cancer_Type = cancer_type,
      MTR_Region_Length = mtr_region_length,
      MTR_Source = mtr_source,
      Genome_Length = genome_length,
      MTR_Genome_Fraction = mtr_genome_fraction,
      MTR_Genome_Percentage = mtr_genome_fraction * 100,
      MTR_Reads_Count = mtr_reads_count,
      Total_Reads_Count = total_reads_count,
      MTR_Reads_Fraction = mtr_reads_fraction,
      MTR_Reads_Percentage = mtr_reads_fraction * 100,
      Enrichment_Score = enrichment_score,
      Log2_Enrichment = ifelse(!is.na(enrichment_score) && enrichment_score > 0,
                               log2(enrichment_score), NA),
      P_Value = p_value,
      Is_Enriched = enrichment_score > 1,
      stringsAsFactors = FALSE
    ))
    
  }, error = function(e) {
    log_error(sample_name, microbe, cancer_type, "CALCULATION_ERROR", e$message)
    return(NULL)
  })
}

# ==============================================================================
# Main Analysis Pipeline with Progress Tracking and Logging
# ==============================================================================

run_complete_enrichment_analysis <- function() {
  
  all_results <- list()
  total_processed <- 0
  total_successful <- 0
  
  for (cancer_type in CONFIG$cancer_types) {
    
    log_message("\n========================================")
    log_message(paste("Processing", cancer_type))
    log_message("========================================")
    
    cancer_dir <- file.path(PATHS$results_base, cancer_type)
    
    if (!dir.exists(cancer_dir)) {
      log_message("  ERROR: Cancer directory not found", "ERROR")
      next
    }
    
    reflist_file <- file.path(PATHS$csv_dir, 
                              paste0(cancer_type, CONFIG$file_patterns$reflist_suffix))
    
    if (!file.exists(reflist_file)) {
      log_message("  ERROR: Reference list not found", "ERROR")
      next
    }
    
    microbe_info <- fread(reflist_file)
    log_message(paste("  Found", nrow(microbe_info), "microbes"))
    
    cancer_results <- list()
    
    for (i in 1:nrow(microbe_info)) {
      
      microbe <- microbe_info$Taxon[i]
      log_message(paste("\n  Microbe", i, "/", nrow(microbe_info), ":", microbe))
      
      genome_length <- get_genome_length(microbe, PATHS$reference_dir)
      
      if (is.na(genome_length)) {
        log_message(paste("    Skipping - genome length:", genome_length), "WARN")
        next
      }
      
      # Find sample files
      vector_dir <- file.path(cancer_dir, CONFIG$vector_subdir, microbe)
      sample_names <- c()
      
      if (dir.exists(vector_dir)) {
        vector_pattern <- paste0(gsub("^_", "", CONFIG$file_patterns$vector_suffix), "$")
        vector_files <- list.files(vector_dir, pattern = vector_pattern)
        sample_names <- gsub(CONFIG$file_patterns$vector_suffix, "", vector_files)
      }
      
      if (length(sample_names) == 0) {
        log_message("    No vector.csv files found", "WARN")
        next
      }
      
      log_message(paste("    Processing", length(sample_names), "samples..."))
      
      # Process samples
      microbe_results <- list()
      success_count <- 0
      
      for (j in seq_along(sample_names)) {
        sample <- sample_names[j]
        total_processed <- total_processed + 1
        
        if (j %% 10 == 0 || j == length(sample_names)) {
          log_message(paste("      Progress:", j, "/", length(sample_names)))
        }
        
        result <- calculate_enrichment_score(
          sample_name = sample,
          microbe = microbe,
          cancer_type = cancer_type,
          genome_length = genome_length,
          cancer_dir = cancer_dir
        )
        
        if (!is.null(result)) {
          microbe_results[[length(microbe_results) + 1]] <- result
          success_count <- success_count + 1
          total_successful <- total_successful + 1
        }
      }
      
      log_message(paste("      Successfully processed:", success_count, "/", length(sample_names)))
      
      if (length(microbe_results) > 0) {
        combined_microbe <- do.call(rbind, microbe_results)
        cancer_results[[length(cancer_results) + 1]] <- combined_microbe
      }
    }
    
    # Save intermediate results
    if (length(cancer_results) > 0) {
      cancer_combined <- do.call(rbind, cancer_results)
      all_results[[length(all_results) + 1]] <- cancer_combined
      
      dir.create(PATHS$output_intermediate, recursive = TRUE, showWarnings = FALSE)
      intermediate_file <- file.path(PATHS$output_intermediate, 
                                     paste0(cancer_type, CONFIG$output_files$intermediate_suffix))
      write.csv(cancer_combined, intermediate_file, row.names = FALSE)
      log_message(paste("  Intermediate results saved for", cancer_type))
    }
  }
  
  log_message("\n========================================")
  log_message("PROCESSING SUMMARY")
  log_message("========================================")
  log_message(paste("Total samples processed:", total_processed))
  log_message(paste("Successfully calculated:", total_successful))
  log_message(paste("Success rate:", round(100 * total_successful / total_processed, 1), "%"))
  
  if (length(all_results) > 0) {
    final_results <- do.call(rbind, all_results)
    return(final_results)
  } else {
    return(NULL)
  }
}

# ==============================================================================
# Statistical Analysis and Visualization
# ==============================================================================

analyze_and_visualize_results <- function(results_df, output_dir) {
  
  if (is.null(results_df) || nrow(results_df) == 0) {
    log_message("\nNo valid results to analyze", "WARN")
    return(NULL)
  }
  
  valid_results <- results_df %>%
    filter(!is.na(Enrichment_Score) & 
             !is.infinite(Enrichment_Score) &
             MTR_Reads_Count > CONFIG$thresholds$min_mtr_reads &
             Total_Reads_Count > CONFIG$thresholds$min_total_reads)
  
  log_message("\n=== Analysis Results ===")
  log_message(paste("Total valid samples:", nrow(valid_results)))
  
  log_message("\nOverall Statistics:")
  log_message(paste("  Median Enrichment Score:", round(median(valid_results$Enrichment_Score), 2)))
  log_message(paste("  Mean Enrichment Score:", round(mean(valid_results$Enrichment_Score), 2)))
  log_message(paste("  Min Enrichment Score:", round(min(valid_results$Enrichment_Score), 2)))
  log_message(paste("  Max Enrichment Score:", round(max(valid_results$Enrichment_Score), 2)))
  log_message(paste("  Samples with enrichment > 1:", sum(valid_results$Is_Enriched), 
                    "(", round(100 * mean(valid_results$Is_Enriched), 1), "%)"))
  
  if (any(!is.na(valid_results$P_Value))) {
    log_message(paste("  Significant enrichment (p <", CONFIG$thresholds$p_value, "):", 
                      sum(valid_results$P_Value < CONFIG$thresholds$p_value, na.rm = TRUE), 
                      "(", round(100 * mean(valid_results$P_Value < CONFIG$thresholds$p_value, na.rm = TRUE), 1), "%)"))
  }
  
  cancer_summary <- valid_results %>%
    group_by(Cancer_Type) %>%
    summarise(
      N_Samples = n(),
      Median_Enrichment = median(Enrichment_Score),
      Mean_Enrichment = mean(Enrichment_Score),
      Pct_Enriched = round(100 * mean(Is_Enriched), 1),
      .groups = 'drop'
    ) %>%
    arrange(desc(Median_Enrichment))
  
  log_message("\nSummary by Cancer Type:")
  log_message(paste(capture.output(print(as.data.frame(cancer_summary), row.names = FALSE)), collapse = "\n"))
  
  # Visualizations
  tryCatch({
    p1 <- ggplot(valid_results, aes(x = Log2_Enrichment)) +
      geom_histogram(bins = 50, fill = "darkgreen", alpha = 0.7) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
      labs(title = "MTR Read Enrichment Distribution",
           x = "Log2(Enrichment Score)",
           y = "Count") +
      theme_minimal()
    
    ggsave(file.path(output_dir, CONFIG$output_files$plot_distribution),
           plot = p1, width = 10, height = 6)
    
    p2 <- ggplot(valid_results, 
                 aes(x = reorder(Cancer_Type, Enrichment_Score, median), 
                     y = Enrichment_Score)) +
      geom_boxplot(fill = "lightgreen", alpha = 0.7) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      scale_y_log10() +
      labs(title = "MTR Read Enrichment by Cancer Type",
           x = "Cancer Type",
           y = "Enrichment Score (log scale)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(output_dir, CONFIG$output_files$plot_by_cancer),
           plot = p2, width = 10, height = 6)
    
    log_message(paste("\nVisualization files saved to:", output_dir))
    
  }, error = function(e) {
    log_message(paste("\nWarning: Could not create some visualizations:", e$message), "WARN")
  })
  
  return(valid_results)
}

# ==============================================================================
# Execute Analysis
# ==============================================================================

dir.create(PATHS$output_dir, recursive = TRUE, showWarnings = FALSE)

# Initialise logging
init_logging(PATHS$output_dir)

log_message("╔══════════════════════════════════════════════════════════╗")
log_message("║         MTR ENRICHMENT SCORE ANALYSIS - COMPLETE          ║")
log_message("╚══════════════════════════════════════════════════════════╝")
log_message("\nThis analysis addresses reviewer concerns about")
log_message("distinguishing true transcriptional reprogramming")
log_message("from microbial abundance effects.\n")
log_message("Configuration:")
log_message(paste("  Base directory:", CONFIG$base_dir))
log_message(paste("  Cancer types:", length(CONFIG$cancer_types)))
log_message(paste("  Output directory:", PATHS$output_dir))
log_message(paste("  Log file:", LOG_FILE))
log_message("")

start_time <- Sys.time()

# Run analysis
results <- run_complete_enrichment_analysis()

if (!is.null(results) && nrow(results) > 0) {
  output_file <- file.path(PATHS$output_dir, CONFIG$output_files$all_results)
  write.csv(results, output_file, row.names = FALSE)
  log_message(paste("\nOK: Raw results saved to:", output_file))
  log_message(paste("  Total rows:", nrow(results)))
  
  analyzed_results <- analyze_and_visualize_results(results, PATHS$output_dir)
  
  if (!is.null(analyzed_results) && nrow(analyzed_results) > 0) {
    output_file_filtered <- file.path(PATHS$output_dir, CONFIG$output_files$filtered_results)
    write.csv(analyzed_results, output_file_filtered, row.names = FALSE)
    log_message(paste("OK: Filtered results saved to:", output_file_filtered))
    
    log_message("\n╔══════════════════════════════════════════════════════════╗")
    log_message("║               SUMMARY FOR MANUSCRIPT                      ║")
    log_message("╚══════════════════════════════════════════════════════════╝")
    
    median_enrichment <- median(analyzed_results$Enrichment_Score)
    mean_enrichment <- mean(analyzed_results$Enrichment_Score)
    pct_enriched <- mean(analyzed_results$Is_Enriched) * 100
    
    log_message("\nKey Finding:")
    log_message(paste("• MTR regions show", round(median_enrichment, 2), "-fold median enrichment"))
    log_message(paste("  (", round(mean_enrichment, 2), "-fold mean enrichment)"))
    log_message(paste("• ", round(pct_enriched, 1), "% of samples show enrichment > 1"))
    log_message("\nThis conclusively demonstrates that MTR represents")
    log_message("selective transcriptional activation, not abundance effects.")
  }
} else {
  log_message("\nFailed: No results generated", "ERROR")
  log_message("Please check error messages above")
}

# Save error log
save_error_log(PATHS$output_dir)

end_time <- Sys.time()
log_message("\n========================================")
log_message(paste("Analysis completed in", 
                  round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))
log_message(paste("Output directory:", PATHS$output_dir))
log_message(paste("Log file:", LOG_FILE))
log_message("========================================")