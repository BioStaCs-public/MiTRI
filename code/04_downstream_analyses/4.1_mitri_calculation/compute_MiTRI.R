setwd("/path/to/your/project/")  # Set working directory

# Load necessary packages
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(data.table)
library(parallel)
library(openxlsx)
library(readxl)

# Define cancer types
cancer_types <- c("THCA")

for (cancer_type in cancer_types) {
  cat("Processing cancer type:", cancer_type, ":\n")
  
  # Define paths
  reflist <- paste0("THCA_TvsN.csv")  # Reference list file
  output_base <- file.path("/path/to/output", cancer_type, "03.coverage")
  output_dir <- file.path("/path/to/output", cancer_type, "04.activity")
  dir.create(output_dir, recursive = TRUE)  # Create directory if it does not exist
  
  reference_base <- file.path("/path/to/reference")
  Reference_length_dir <- file.path("/path/to/reference", reflist)
  reference_lengths <- fread(Reference_length_dir, colClasses = list(numeric = "Reference_length"))
  
  control_base <- "/path/to/control"
  
  # Define sample scenarios
  sample_scenarios <- list(
    "Tumor" = list(
      tumor_dir_base = file.path(output_base, "Tumor"),
      flagstat_base = file.path("data", cancer_type, "Tumor/flagstat/"),
      non_host_ratio_base = file.path("data", cancer_type, "Tumor/reads_ratio/")
    ),
    "Normal" = list(
      tumor_dir_base = file.path(output_base, "Normal"),
      flagstat_base = file.path("data", cancer_type, "Normal/flagstat/"),
      non_host_ratio_base = file.path("data", cancer_type, "Normal/reads_ratio/")
    )
  )
  
  # Record script start time
  start_time <- Sys.time()
  cat("Script started at:", format(start_time), "\n")
  
  generate_boolean_vector <- function(coverage) {
    seq_names <- unique(seqnames(coverage))
    genome_length <- sum(sapply(split(coverage, seqnames(coverage)), function(gr) max(end(gr))))
    boolean_vector <- rep(0L, genome_length)
    offset <- 0
    
    for (chr in seq_names) {
      chr_coverage <- coverage[seqnames(coverage) == chr]
      chr_length <- max(end(chr_coverage))
      valid_ranges <- chr_coverage[chr_coverage$score > 0]
      
      starts <- start(valid_ranges) + offset
      ends <- end(valid_ranges) + offset
      
      for (i in seq_along(starts)) {
        boolean_vector[starts[i]:ends[i]] <- 1L
      }
      
      offset <- offset + chr_length
    }
    return(boolean_vector)
  }
  
  remove_overlap_regions <- function(t_boolean, other_boolean) {
    if (length(t_boolean) != length(other_boolean)) {
      stop("Boolean vectors do not have the same length")
    }
    result <- integer(length(t_boolean))
    
    overlap <- (t_boolean == 1) & (other_boolean == 1)  # Overlap between two coverage vectors
    only_red <- (t_boolean == 1) & (other_boolean == 0) # Only red coverage
    only_blue <- (t_boolean == 0) & (other_boolean == 1) # Only blue coverage
    neither <- (t_boolean == 0) & (other_boolean == 0)  # Neither coverage
    
    result[overlap] <- -1
    result[only_red] <- 1
    result[only_blue] <- 2
    result[neither] <- 0
    
    return(result)
  }
  
  calculate_remaining_coverage <- function(coverage, boolean_vector, non_host_ratio, activity_region_length) {
    coverage_values <- coverage(coverage, weight = "score")
    total_coverage <- 0
    global_pos_index <- 1
    
    for (i in seq_along(coverage_values)) {
      coverage_rle <- coverage_values[[i]]
      run_lengths <- runLength(coverage_rle)
      run_values <- runValue(coverage_rle)
      chr_length <- sum(run_lengths)
      
      start_positions <- global_pos_index + cumsum(c(0, run_lengths[-length(run_lengths)]))
      end_positions <- pmin(start_positions + run_lengths - 1, length(boolean_vector))
      
      for (j in seq_along(run_lengths)) {
        boolean_slice <- boolean_vector[start_positions[j]:end_positions[j]]
        active_bases <- sum(boolean_slice == 1)
        
        if (active_bases > 0) {
          total_coverage <- total_coverage + run_values[j] * active_bases
        }
      }
      
      global_pos_index <- global_pos_index + chr_length
    }
    
    ratio <- (total_coverage / non_host_ratio) / activity_region_length
    return(list(ratio = ratio))
  }
  
  read_non_host_ratio <- function(sample_name) {
    reads_ratio_file <- paste0(non_host_ratio_base, sample_name, "_non_host_ratio.txt")
    
    if (file.exists(reads_ratio_file)) {
      lines <- readLines(reads_ratio_file)
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
  
  # Main loop: iterate over each sample scenario
  results_list <- list()
  for (scenario_name in names(sample_scenarios)) {
    sample_scenario <- sample_scenarios[[scenario_name]]
    
    flagstat_base <- sample_scenario$flagstat_base
    tumor_dir_base <- sample_scenario$tumor_dir_base
    non_host_ratio_base <- sample_scenario$non_host_ratio_base
    
    reference_dir <- file.path(reference_base)
    control_dir_base <- file.path(control_base)
    
    fna_files <- list.files(reference_dir, pattern = "\\.fna$", full.names = TRUE)
    reference_names <- basename(fna_files)
    reference_names <- gsub("\\.fna$", "", reference_names)
    
    valid_reference_names <- reference_names[reference_names %in% reference_lengths$Taxon]
    
    results_dt <- rbindlist(mclapply(valid_reference_names, function(reference_name) {
      tryCatch({
        tumor_dir <- file.path(tumor_dir_base, reference_name)
        brain_dir <- file.path(control_dir_base, reference_name, "brain")
        gbm_dir <- file.path(control_dir_base, reference_name, "GBM_Normal")
        brain_combined_dir <- brain_dir
        
        brain_file <- file.path(brain_dir, "coverage_brain.bedgraph")
        brain_combined_file <- file.path(brain_combined_dir, "coverage_brain_combined.bedgraph")
        gbm_file <- file.path(gbm_dir, "coverage_GBM_Normal.bedgraph")
        testicle_file <- file.path(control_dir_base, reference_name, "testicle_6", "coverage_testicle.bedgraph")
        
        t_files <- list.files(tumor_dir, pattern = ".*_rna_.*\\.bedgraph$", full.names = TRUE)
        
        # Handle empty files
        zero_size_files <- t_files[file.info(t_files)$size == 0]
        if (length(zero_size_files) > 0) {
          cat("Reference", reference_name, ":\n")
          print(basename(zero_size_files))
          t_files <- t_files[file.info(t_files)$size > 0]  # Keep only non-zero files
        }
        
        t_coverages <- lapply(t_files, import)
        
        brain_coverage <- import(brain_file)
        gbm_coverage <- import(gbm_file)
        brain_combined_coverage <- import(brain_combined_file)
        testicle_coverage <- import(testicle_file)
        
        t_boolean_vectors <- lapply(t_coverages, generate_boolean_vector)
        brain_boolean_vector <- generate_boolean_vector(brain_coverage)
        gbm_boolean_vector <- generate_boolean_vector(gbm_coverage)
        brain_combined_boolean_vector <- generate_boolean_vector(brain_combined_coverage)
        testicle_boolean_vector <- generate_boolean_vector(testicle_coverage)
        
        calculate_coverage_for_control <- function(control_vector) {
          lapply(seq_along(t_boolean_vectors), function(i) {
            file_name <- basename(t_files[i])
            sample_name <- strsplit(file_name, "_rna")[[1]][1]
            
            cleaned_vector <- remove_overlap_regions(t_boolean_vectors[[i]], control_vector)
            non_host_ratio <- read_non_host_ratio(sample_name)
            activity_region_length <- sum(cleaned_vector != 0)
            remaining_coverage <- calculate_remaining_coverage(t_coverages[[i]], cleaned_vector, non_host_ratio, activity_region_length)
            
            list(sample_name = sample_name, ratio = remaining_coverage$ratio)
          })
        }
        
        brain_results <- calculate_coverage_for_control(brain_boolean_vector)
        gbm_results <- calculate_coverage_for_control(gbm_boolean_vector)
        brain_combined_results <- calculate_coverage_for_control(brain_combined_boolean_vector)
        testicle_results <- calculate_coverage_for_control(testicle_boolean_vector)
        
        sample_names <- sapply(brain_results, function(x) x$sample_name)
        Signal_Brain_ratio_values <- sapply(brain_results, function(x) x$ratio)
        Signal_GBM_ratio_values <- sapply(gbm_results, function(x) x$ratio)
        Signal_Brain_Combined_ratio_values <- sapply(brain_combined_results, function(x) x$ratio)
        Signal_Testicle_ratio_values <- sapply(testicle_results, function(x) x$ratio)
        
        data.table(
          Reference = reference_name,
          Scenario = scenario_name,
          Sample_Name = sample_names,
          'Brain 1(n=17)' = Signal_Brain_ratio_values,
          'Brain 2(n=9)'  = Signal_GBM_ratio_values,
          'Brain 1\u222A2(n=26)' = Signal_Brain_Combined_ratio_values,
          'Testis(n=6)' = Signal_Testicle_ratio_values
        )
      }, error = function(e) {
        cat("Error processing reference:", reference_name, "\n", e$message, "\n")
        return(data.table(
          Reference = reference_name,
          Scenario = scenario_name,
          Sample_Name = NA,
          Error = e$message
        ))
      })
    }, mc.cores = 20), use.names = TRUE, fill = TRUE)
    
    results_list[[paste(scenario_name)]] <- results_dt
  }
  
  # Save results to Excel
  output_xlsx_file <- file.path(output_dir, paste0(cancer_type, "_Tumor_Normal_activity_ratiolist_spp2.xlsx"))
  wb <- createWorkbook()
  for (sheet_name in names(results_list)) {
    truncated_name <- substr(sheet_name, 1, 31)
    addWorksheet(wb, truncated_name)
    writeDataTable(wb, truncated_name, results_list[[sheet_name]])
  }
  
  tryCatch({
    saveWorkbook(wb, output_xlsx_file, overwrite = TRUE)
    cat("Results saved to:", output_xlsx_file, "\n")
  }, error = function(e) {
    cat("Error saving Excel file:", e$message, "\n")
  })
  
  # Clean up memory
  rm(results_list, wb, output_xlsx_file)
  gc()  # Garbage collection
}

end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "secs")
elapsed_time_formatted <- sprintf("%02d:%02d:%02d", as.integer(elapsed_time) %/% 3600, (as.integer(elapsed_time) %% 3600) %/% 60, as.integer(elapsed_time) %% 60)

cat("Total execution time:", elapsed_time_formatted, "\n")
