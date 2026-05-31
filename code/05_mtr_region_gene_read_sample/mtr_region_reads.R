# Load required libraries
library(dplyr)
library(readr)
library(stringr)
library(Biostrings)
library(data.table)
library(parallel)

# ============================================
# Core function: CIGAR parsing
# ============================================
calculate_reference_length <- function(cigar_string) {
  if (is.na(cigar_string) || cigar_string == "*") {
    return(NA_integer_)
  }

  operations <- str_extract_all(cigar_string, "\\d+[MIDNSHPX=]")[[1]]

  if (length(operations) == 0) {
    return(NA_integer_)
  }

  ref_length <- 0

  for (op in operations) {
    len <- as.integer(str_extract(op, "\\d+"))
    type <- str_extract(op, "[MIDNSHPX=]")

    if (type %in% c("M", "D", "N", "=", "X")) {
      ref_length <- ref_length + len
    }
  }

  return(ref_length)
}

calculate_reference_length_vectorized <- function(cigar_vector) {
  sapply(cigar_vector, calculate_reference_length)
}

# ============================================
# FNA parser (new format - includes gene annotations)
# ============================================
parse_fna_with_annotations <- function(fna_file) {
  sequences <- readDNAStringSet(fna_file)
  headers <- names(sequences)

  # Parse each header
  regions_list <- lapply(headers, function(header) {
    # Separate position info from annotation info
    # Format: seqid:start-end gene_id=xxx gene_name=xxx product=xxx strand=x length=xxxbp
    parts <- strsplit(header, " ", fixed = TRUE)[[1]]

    # Extract position info (first field)
    position_part <- parts[1]

    # Parse seqid:start-end
    if (grepl(":", position_part)) {
      coord_parts <- strsplit(position_part, ":", fixed = TRUE)[[1]]
      if (length(coord_parts) == 2) {
        seqid <- coord_parts[1]
        range_part <- coord_parts[2]

        # Parse start-end
        if (grepl("-", range_part)) {
          range_coords <- strsplit(range_part, "-", fixed = TRUE)[[1]]
          if (length(range_coords) == 2) {
            start_pos <- as.integer(range_coords[1])
            end_pos <- as.integer(range_coords[2])

            # Initialise annotation fields
            gene_id <- NA_character_
            gene_name <- NA_character_
            locus_tag <- NA_character_
            product <- NA_character_
            biotype <- NA_character_
            strand <- NA_character_
            gene_length <- NA_integer_

            # Parse annotation info (remaining fields)
            if (length(parts) > 1) {
              annotations <- parts[2:length(parts)]

              for (anno in annotations) {
                if (grepl("=", anno)) {
                  kv <- strsplit(anno, "=", fixed = TRUE)[[1]]
                  if (length(kv) == 2) {
                    key <- kv[1]
                    value <- kv[2]

                    if (key == "gene_id") {
                      gene_id <- value
                    } else if (key == "gene_name" || key == "gene") {
                      gene_name <- value
                    } else if (key == "locus_tag") {
                      locus_tag <- value
                    } else if (key == "product") {
                      product <- value
                    } else if (key == "biotype") {
                      biotype <- value
                    } else if (key == "strand") {
                      strand <- value
                    } else if (key == "length") {
                      gene_length <- as.integer(gsub("bp", "", value))
                    }
                  }
                }
              }
            }

            return(data.frame(
              Reference = seqid,
              Start = start_pos,
              End = end_pos,
              gene_id = gene_id,
              gene_name = gene_name,
              locus_tag = locus_tag,
              product = product,
              biotype = biotype,
              strand = strand,
              gene_length = gene_length,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }

    # Return NULL if parsing failed
    return(NULL)
  })

  # Drop NULLs and combine
  regions_list <- regions_list[!sapply(regions_list, is.null)]

  if (length(regions_list) == 0) {
    stop("Could not parse any FNA header")
  }

  regions <- do.call(rbind, regions_list)
  return(regions)
}

# ============================================
# FLAG-parsing helpers
# ============================================
is_proper_pair <- function(flag) {
  bitwAnd(as.integer(flag), 2) > 0
}

is_first_in_pair <- function(flag) {
  bitwAnd(as.integer(flag), 64) > 0
}

is_second_in_pair <- function(flag) {
  bitwAnd(as.integer(flag), 128) > 0
}

is_reverse_strand <- function(flag) {
  bitwAnd(as.integer(flag), 16) > 0
}

# ============================================
# Per-sample processing (includes gene annotations)
# ============================================
process_sample <- function(microbe, fna_sample, bam_txt_file, fna_file, output_dir) {

  # Initialise the return value
  result <- list(
    microbe = microbe,
    sample = fna_sample,
    success = FALSE,
    fully_contained = 0L,
    boundary = 0L,
    total = 0L,
    error_message = NULL
  )

  tryCatch({
    # -------------------- Read and prepare data --------------------

    # Parse regions (includes gene annotations)
    regions <- parse_fna_with_annotations(fna_file)

    if (is.null(regions) || nrow(regions) == 0L) {
      result$error_message <- "regions is empty"
      return(result)
    }

    message("      Parsed ", nrow(regions), " gene regions (with annotations)")

    # Read BAM data
    bam_data <- read.table(bam_txt_file, sep = "\t", header = FALSE,
                           comment.char = "@", fill = TRUE, stringsAsFactors = FALSE)

    if (is.null(bam_data) || nrow(bam_data) == 0L) {
      result$error_message <- "bam data is empty"
      return(result)
    }

    # Standardise column names
    nkeep <- min(11L, ncol(bam_data))
    bam_data <- bam_data[, seq_len(nkeep), drop = FALSE]
    if (ncol(bam_data) < 11L) {
      for (i in (ncol(bam_data) + 1):11) {
        bam_data[[i]] <- NA
      }
    }
    colnames(bam_data)[1:11] <- c("QNAME","FLAG","RNAME","POS","MAPQ",
                                  "CIGAR","RNEXT","PNEXT","TLEN","SEQ","QUAL")

    # Convert to data.table and clean
    bam_dt <- as.data.table(bam_data)
    bam_dt <- bam_dt[!is.na(POS) & !is.na(CIGAR) & CIGAR != "*"]

    if (nrow(bam_dt) == 0L) {
      result$error_message <- "no valid data"
      return(result)
    }

    # Type conversion
    bam_dt[, `:=`(
      FLAG = as.integer(FLAG),
      POS = as.integer(POS),
      PNEXT = as.integer(PNEXT),
      TLEN = as.integer(TLEN),
      CIGAR = as.character(CIGAR)
    )]

    # -------------------- Compute ref_length --------------------
    message("[", microbe, "][", fna_sample, "] computing CIGAR lengths...")
    bam_dt[, ref_length := calculate_reference_length_vectorized(CIGAR)]
    bam_dt <- bam_dt[!is.na(ref_length) & ref_length > 0]

    if (nrow(bam_dt) == 0L) {
      result$error_message <- "no valid CIGAR"
      return(result)
    }

    # End position of each read
    bam_dt[, End_POS := POS + ref_length - 1L]

    # Add strand info
    bam_dt[, is_reverse := is_reverse_strand(FLAG)]

    # Prepare regions
    regions_dt <- as.data.table(regions)
    regions_dt[, `:=`(Start = as.integer(Start), End = as.integer(End))]
    regions_dt <- regions_dt[!is.na(Start) & !is.na(End) & End >= Start]

    # ============================================
    # Fully-contained + boundary statistics
    # ============================================
    message("[", microbe, "][", fna_sample, "] computing fully-contained matches + boundary stats...")

    # Coarse filter
    regions_bounds <- regions_dt[, .(
      minS = min(Start, na.rm = TRUE),
      maxE = max(End, na.rm = TRUE)
    ), by = Reference]

    bam_filtered <- merge(bam_dt, regions_bounds,
                          by.x = "RNAME", by.y = "Reference",
                          all.x = FALSE, all.y = FALSE)

    bam_filtered <- bam_filtered[End_POS >= minS & POS <= maxE]

    if (nrow(bam_filtered) == 0) {
      result$error_message <- "no data after coarse filter"
      return(result)
    }

    # Fine filter (with gene annotation)
    setkey(regions_dt, Reference)
    setkey(bam_filtered, RNAME)

    all_overlaps <- bam_filtered[regions_dt, allow.cartesian = TRUE, nomatch = 0]

    # Classify overlap type
    all_overlaps[, `:=`(
      is_fully_contained = (POS >= Start & End_POS <= End),
      is_any_overlap = (
        (POS >= Start & POS <= End) |
          (End_POS >= Start & End_POS <= End) |
          (POS <= Start & End_POS >= End)
      )
    )]

    all_overlaps[, is_boundary := is_any_overlap & !is_fully_contained]

    # -------------------- Extract results (with gene annotation) --------------------

    result_fully_contained <- all_overlaps[is_fully_contained == TRUE,
                                           .(QNAME, FLAG, RNAME, POS, End_POS, CIGAR, PNEXT, TLEN,
                                             ref_length, is_reverse, SEQ,
                                             Region_Start = Start, Region_End = End,
                                             gene_id, gene_name, locus_tag, product, biotype, strand, gene_length,
                                             overlap_type = "fully_contained")
    ]

    result_boundary <- all_overlaps[is_boundary == TRUE,
                                    .(QNAME, FLAG, RNAME, POS, End_POS, CIGAR, PNEXT, TLEN,
                                      ref_length, is_reverse, SEQ,
                                      Region_Start = Start, Region_End = End,
                                      gene_id, gene_name, locus_tag, product, biotype, strand, gene_length,
                                      overlap_type = "boundary")
    ]

    # -------------------- Write results --------------------
    base_name <- sub("\\.bam\\.txt$", "", basename(bam_txt_file))
    sample_name <- sub("\\..*", "", basename(base_name))

    # Write fully-contained reads (with gene annotation)
    if (nrow(result_fully_contained) > 0) {
      out_contained <- file.path(output_dir,
                                 paste0(sample_name, "_fully_contained_reads.csv"))
      fwrite(result_fully_contained, out_contained)
      message("      -> fully contained: ", nrow(result_fully_contained), " reads")
    } else {
      message("      -> fully contained: 0 reads")
    }

    # Write boundary reads (with gene annotation)
    if (nrow(result_boundary) > 0) {
      out_boundary <- file.path(output_dir,
                                paste0(sample_name, "_boundary_reads.csv"))
      fwrite(result_boundary, out_boundary)
      message("      -> boundary reads: ", nrow(result_boundary), " reads")
    } else {
      message("      -> boundary reads: 0 reads")
    }

    # Write summary statistics (grouped by gene annotation)
    summary_file <- file.path(output_dir, paste0(sample_name, "_summary.txt"))

    # Per-gene-region statistics (with gene annotation)
    gene_stats <- all_overlaps[, .(
      fully_contained = sum(is_fully_contained),
      boundary = sum(is_boundary),
      total_overlapping = sum(is_any_overlap)
    ), by = .(RNAME, Start, End, gene_id, gene_name, locus_tag, product)]

    # Rename for output consistency
    setnames(gene_stats, c("Start", "End"), c("Region_Start", "Region_End"))

    gene_stats_file <- file.path(output_dir,
                                 paste0(sample_name, "_gene_stats.csv"))
    fwrite(gene_stats, gene_stats_file)

    # Overall statistics
    sink(summary_file)
    cat("Sample:", fna_sample, "\n")
    cat("Microbe:", microbe, "\n")
    cat("Gene regions:", nrow(regions_dt), "\n\n")
    cat("=" , rep("=", 50), "\n", sep = "")
    cat("Read statistics:\n")
    cat("  fully-contained reads:", nrow(result_fully_contained), "\n")
    cat("  boundary reads:", nrow(result_boundary), "\n")
    cat("  total overlapping reads:", nrow(result_fully_contained) + nrow(result_boundary), "\n\n")
    cat("=" , rep("=", 50), "\n", sep = "")
    cat("\nPer-gene-region statistics (top 10, sorted by fully-contained count):\n")
    print(head(gene_stats[order(-fully_contained)], 10))
    sink()

    message("      -> summary file: ", basename(summary_file))
    message("      -> gene statistics: ", basename(gene_stats_file))

    # Update the return value
    result$success <- TRUE
    result$fully_contained <- nrow(result_fully_contained)
    result$boundary <- nrow(result_boundary)
    result$total <- nrow(result_fully_contained) + nrow(result_boundary)

  }, error = function(e) {
    result$error_message <<- as.character(e$message)
    message("      Error: ", e$message)
  })

  return(result)
}

# ============================================
# Function to process a single microbe
# ============================================
process_microbe <- function(microbe, base_dirs) {
  message("\n========== Processing microbe: ", microbe, " ==========")

  coverage_dir <- file.path(base_dirs$coverage, microbe)
  fna_dir <- file.path(base_dirs$fna, microbe)
  output_dir <- file.path(base_dirs$output, microbe)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  bam_txt_files <- list.files(coverage_dir, pattern = "\\.bam\\.txt$", full.names = TRUE)
  fna_files <- list.files(fna_dir, pattern = "_target_genes_gff\\.fna$", full.names = TRUE)

  # Filter out empty files
  if (length(fna_files)) {
    finfo <- file.info(fna_files)
    fna_files <- fna_files[!is.na(finfo$size) & finfo$size > 0]
  }

  if (length(fna_files) == 0L) {
    message("[", microbe, "] no valid fna file")
    return(data.frame(
      microbe = microbe,
      sample = NA,
      success = FALSE,
      fully_contained = 0L,
      boundary = 0L,
      total = 0L,
      error_message = "no valid fna file",
      stringsAsFactors = FALSE
    ))
  }

  # Match sample names
  fna_samples <- sub("_target_genes_.*", "", basename(fna_files))
  bam_samples <- sub("_.*", "", basename(bam_txt_files))
  bam_txt_files <- bam_txt_files[bam_samples %in% fna_samples]

  if (length(bam_txt_files) == 0L) {
    message("[", microbe, "] no matching bam file")
    return(data.frame(
      microbe = microbe,
      sample = NA,
      success = FALSE,
      fully_contained = 0L,
      boundary = 0L,
      total = 0L,
      error_message = "no matching bam file",
      stringsAsFactors = FALSE
    ))
  }

  message("  FNA file count: ", length(fna_files))
  message("  matching BAM count: ", length(bam_txt_files))

  # Process each sample
  sample_results <- list()
  for (fna_file in fna_files) {
    fna_sample <- sub("_target_genes_.*", "", basename(fna_file))
    bam_samples <- sub("_.*", "", basename(bam_txt_files))
    cand <- bam_txt_files[bam_samples == fna_sample]

    if (length(cand) == 1L) {
      message("\n  Processing sample: [", fna_sample, "]")
      message("    FNA: ", basename(fna_file))
      message("    BAM: ", basename(cand))

      res <- process_sample(microbe, fna_sample, cand, fna_file, output_dir)
      sample_results[[length(sample_results) + 1]] <- res
    }
  }

  # Convert to data.frame
  if (length(sample_results) > 0) {
    results_df <- do.call(rbind, lapply(sample_results, function(x) {
      data.frame(
        microbe = x$microbe,
        sample = x$sample,
        success = x$success,
        fully_contained = x$fully_contained,
        boundary = x$boundary,
        total = x$total,
        error_message = ifelse(is.null(x$error_message), "", x$error_message),
        stringsAsFactors = FALSE
      )
    }))
    return(results_df)
  } else {
    return(data.frame(
      microbe = microbe,
      sample = NA,
      success = FALSE,
      fully_contained = 0L,
      boundary = 0L,
      total = 0L,
      error_message = "no sample processed",
      stringsAsFactors = FALSE
    ))
  }
}

# ============================================
# Main program - multiple cancers + R/NR loop
# ============================================
#
# # List of cancer types to process
# cancer_list <- c("CRC", "BRCA", "CESC", "LIHC", "OSCC", "LUNG", "PAAD", "STAD", "THCA")
# List of cancer types to process
cancer_list <- c("BLCA")

# Response types
response_types <- c("Tumor")

# Base paths
base_path <- "/path/to/project/results_V3/cancers_V3.1"
csv_path <- "/path/to/project/data_V3/CSVs_20251208"

# Loop over each cancer type
for (cancer in cancer_list) {

  message("\n")
  message("================================================================================")
  message("Start processing cancer type: ", cancer)
  message("================================================================================")

  # Loop over response types (R and NR)
  for (response_type in response_types) {

    message("\n")
    message("========================================")
    message("Processing: ", cancer, " - ", response_type)
    message("========================================")

    # Read the corresponding microbe-list CSV file
    csv_file <- file.path(csv_path, paste0(cancer, "_V4.csv"))

    if (!file.exists(csv_file)) {
      message("[", cancer, "-", response_type, "] CSV file does not exist: ", csv_file)
      next
    }

    # Read the microbe list
    tryCatch({
      microbe_df <- read.csv(csv_file)
      if (!"Taxon" %in% colnames(microbe_df)) {
        message("[", cancer, "-", response_type, "] CSV file is missing the Taxon column")
        next
      }
      microbe_list <- microbe_df$Taxon

      if (length(microbe_list) == 0) {
        message("[", cancer, "-", response_type, "] microbe list is empty")
        next
      }

      message("Read ", length(microbe_list), " microbes from CSV")

    }, error = function(e) {
      message("[", cancer, "-", response_type, "] failed to read CSV file: ", e$message)
      next
    })

    # Set paths for this cancer type + response type
    base_dirs <- list(
      coverage = file.path(base_path, cancer, "03.coverage", response_type),
      fna = file.path(base_path, cancer, "05.MTR/05.3.MTR_genes/2.filtered_genes_no_merge"),
      output = file.path(base_path, cancer, "05.MTR/05.4.MTR_reads",
                         paste0(response_type, "_fully_contained"))
    )

    # Check that the paths exist
    if (!dir.exists(base_dirs$fna)) {
      message("[", cancer, "-", response_type, "] FNA path does not exist: ", base_dirs$fna)
      next
    }

    if (!dir.exists(base_dirs$coverage)) {
      message("[", cancer, "-", response_type, "] Coverage path does not exist: ", base_dirs$coverage)
      next
    }

    # Create the output directories
    if (!dir.exists(dirname(base_dirs$output))) {
      dir.create(dirname(base_dirs$output), recursive = TRUE)
    }
    if (!dir.exists(base_dirs$output)) {
      dir.create(base_dirs$output, recursive = TRUE)
    }

    message("========== Start parallel processing ==========")
    message("Cancer type: ", cancer)
    message("Response type: ", response_type)
    message("Total microbes: ", length(microbe_list))
    message("Parallel cores: 32")

    # Process each microbe in parallel
    result_list <- mclapply(microbe_list, function(microbe) {
      process_microbe(microbe, base_dirs)
    }, mc.cores = 32)

    # ============================================
    # Aggregate the statistics
    # ============================================
    message("\n========== Aggregate the statistics ==========")

    # Combine all results
    total_stats <- do.call(rbind, result_list)

    # Filter out failed samples
    successful_stats <- total_stats[total_stats$success == TRUE, ]
    failed_stats <- total_stats[total_stats$success == FALSE, ]

    # Save all results
    output_summary_dir <- file.path(base_path, cancer, "05.MTR/05.4.MTR_reads", response_type)
    if (!dir.exists(output_summary_dir)) {
      dir.create(output_summary_dir, recursive = TRUE)
    }

    # Save the complete statistics
    write.csv(total_stats,
              file.path(output_summary_dir,
                        paste0(cancer, "_", response_type, "_overall_summary_complete.csv")),
              row.names = FALSE)

    # Save the successful statistics
    if (nrow(successful_stats) > 0) {
      write.csv(successful_stats,
                file.path(output_summary_dir,
                          paste0(cancer, "_", response_type, "_overall_summary_success.csv")),
                row.names = FALSE)
    }

    # Save the failed statistics
    if (nrow(failed_stats) > 0) {
      write.csv(failed_stats,
                file.path(output_summary_dir,
                          paste0(cancer, "_", response_type, "_overall_summary_failed.csv")),
                row.names = FALSE)
    }

    # Print the summary
    message("\n" , rep("=", 60))
    message("[", cancer, "-", response_type, "] Overall statistics:")
    message("  microbes processed: ", length(unique(total_stats$microbe)))
    message("  total samples processed: ", nrow(total_stats))
    message("  successful samples: ", nrow(successful_stats))
    message("  failed samples: ", nrow(failed_stats))
    message(rep("=", 60))

    if (nrow(successful_stats) > 0) {
      message("\n[", cancer, "-", response_type, "] Successful-sample statistics:")
      message("  total fully-contained reads: ", sum(successful_stats$fully_contained, na.rm = TRUE))
      message("  total boundary reads: ", sum(successful_stats$boundary, na.rm = TRUE))
      message("  total overlapping reads: ", sum(successful_stats$total, na.rm = TRUE))
    }

    if (nrow(failed_stats) > 0) {
      message("\n[", cancer, "-", response_type, "] Failed-sample details:")
      print(failed_stats[, c("microbe", "sample", "error_message")])
    }

    message("\n[", cancer, "-", response_type, "] Files saved:")
    message("  - ", paste0(cancer, "_", response_type, "_overall_summary_complete.csv"), " (all samples)")
    message("  - ", paste0(cancer, "_", response_type, "_overall_summary_success.csv"), " (successful samples)")
    message("  - ", paste0(cancer, "_", response_type, "_overall_summary_failed.csv"), " (failed samples)")
    message(rep("=", 60))

    message("\n[", cancer, "-", response_type, "] ========== Done ==========\n")
  }

  message("\n[", cancer, "] ========== All response types for this cancer type are done ==========\n")
}

message("\n")
message("================================================================================")
message("All cancer types and response types are done")
message("================================================================================")
