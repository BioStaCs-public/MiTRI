# ---- Setup ----
setwd("/project")

# ---- Dependencies (kept minimal, no unused packages) ----
library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(parallel)
library(Biostrings)

# ---- Config (anonymized) ----
cancer_type <- c("PROJECT_TYPE")
# reference list CSV
reflist <- "project_reference.csv"

output_base <- file.path("/project/results/special", cancer_type, "06.coverage_target")
output_dir  <- file.path("/project/results/special", cancer_type, "04.activity_region_T_Nback")
dir.create(output_dir, recursive = TRUE)

reference_base <- file.path("/project/data/reference")
reference_length_csv <- file.path("/project/data/csvs_special", reflist)
reference_lengths <- fread(reference_length_csv, colClasses = list(numeric = "Reference_length"))

control_base <- "/project/results/control"

# Sample scenarios (kept fields for compatibility, even if not used elsewhere)
sample_scenarios <- list(
  "T" = list(
    tumor_dir_base      = file.path(output_base, "T"),
    flagstat_base       = file.path("data/special", cancer_type, "T/flagstat/"),
    non_host_ratio_base = file.path("data/special", cancer_type, "T/reads_ratio/")
  )
)

normal_dir_base <- "/project/results/cancers/PROJECT/03.coverage/Normal"

# ---- Logging ----
start_time <- Sys.time()
cat("Script start time:", format(start_time), "\n")

# ---- Helpers (logic unchanged) ----
generate_boolean_vector <- function(coverage) {
  # Build a flattened boolean vector across contigs based on score > 0
  seq_names <- sort(unique(seqnames(coverage)))
  genome_length <- sum(sapply(split(coverage, seqnames(coverage)), function(gr) max(end(gr))))
  boolean_vector <- rep(0L, genome_length)
  offset <- 0

  for (chr in seq_names) {
    chr_coverage <- coverage[seqnames(coverage) == chr]
    chr_length <- max(end(chr_coverage))
    valid_ranges <- chr_coverage[chr_coverage$score > 0]
    starts <- start(valid_ranges) + offset
    ends   <- end(valid_ranges)   + offset

    for (i in seq_along(starts)) {
      boolean_vector[starts[i]:ends[i]] <- 1L
    }
    offset <- offset + chr_length
  }
  return(boolean_vector)
}

remove_overlap_regions <- function(t_boolean, other_boolean) {
  if (length(t_boolean) != length(other_boolean)) {
    stop("Lengths of t_boolean and other_boolean do not match.")
  }
  result <- integer(length(t_boolean))
  overlap   <- (t_boolean == 1) & (other_boolean == 1)
  only_red  <- (t_boolean == 1) & (other_boolean == 0)
  only_blue <- (t_boolean == 0) & (other_boolean == 1)
  neither   <- (t_boolean == 0) & (other_boolean == 0)

  result[overlap]   <- 0
  result[only_red]  <- 1
  result[only_blue] <- 0
  result[neither]   <- 0
  return(result)
}

# ---- Main ----
for (scenario_name in names(sample_scenarios)) {
  sample_scenario <- sample_scenarios[[scenario_name]]

  flagstat_base       <- sample_scenario$flagstat_base
  tumor_dir_base      <- sample_scenario$tumor_dir_base
  non_host_ratio_base <- sample_scenario$non_host_ratio_base

  reference_dir <- file.path(reference_base)
  control_dir_base <- file.path(control_base)

  fna_files <- list.files(reference_dir, pattern = "\\.fna$", full.names = TRUE)
  reference_names <- basename(fna_files)
  reference_names <- gsub("\\.fna$", "", reference_names)

  valid_reference_names <- reference_names[reference_names %in% reference_lengths$Taxon]

  result <- mclapply(valid_reference_names, function(reference_name) {
    tryCatch({
      tumor_dir  <- file.path(tumor_dir_base, reference_name)
      normal_dir <- file.path(normal_dir_base, reference_name)

      n_file <- file.path(normal_dir, "coverage_Normal.bedgraph")
      n_coverage <- import(n_file)

      t_files <- list.files(tumor_dir, pattern = ".*_rna_.*\\.bedgraph$", full.names = TRUE)
      t_coverages <- lapply(t_files, import)

      t_boolean_vectors <- lapply(t_coverages, generate_boolean_vector)
      n_boolean_vectors <- generate_boolean_vector(n_coverage)

      lapply(seq_along(t_boolean_vectors), function(i) {
        file_name   <- basename(t_files[i])
        sample_name <- strsplit(file_name, "_rna")[[1]][1]

        cleaned_vector <- remove_overlap_regions(t_boolean_vectors[[i]], n_boolean_vectors)

        # Save boolean vector
        output_path_vector <- file.path(output_dir, "04.1.vector", reference_name)
        dir.create(output_path_vector, recursive = TRUE)
        output_csv <- file.path(output_path_vector, paste0(sample_name, "_all_vector.csv"))
        write.csv(cleaned_vector, output_csv, row.names = FALSE)

        # Extract sequences according to boolean runs (logic kept exactly as provided)
        fna_file <- file.path(reference_dir, paste0(reference_name, ".fna"))
        reference_dna <- readDNAStringSet(fna_file)

        seq_names <- names(reference_dna)
        seq_names <- sort(seq_names)
        offset <- 0
        all_extracted_sequences <- DNAStringSet()

        for (chr in seq_names) {
          chr_dna <- reference_dna[[chr]]

          # NOTE: kept as in your original script (no change to indexing expression)
          chr_boolean_vector <- cleaned_vector[offset + 1 : length(chr_dna)]

          runs <- rle(chr_boolean_vector)
          segments <- which(runs$values == 1)

          for (segment in segments) {
            start_pos <- sum(runs$lengths[1:(segment - 1)]) + 1
            end_pos   <- start_pos + runs$lengths[segment] - 1

            extracted_sequence <- subseq(chr_dna, start = start_pos, end = end_pos)

            id <- strsplit(chr, " ")[[1]][1]
            new_names <- paste0(id, "_", start_pos, "_", end_pos)

            extracted_sequence_set <- DNAStringSet(extracted_sequence)
            names(extracted_sequence_set) <- new_names
            all_extracted_sequences <- c(all_extracted_sequences, extracted_sequence_set)
          }
          offset <- offset + length(chr_dna)
        }

        output_dir_region <- file.path(output_dir, "04.2.region", reference_name)
        dir.create(output_dir_region, recursive = TRUE)
        output_seq_file <- file.path(output_dir_region, paste0(sample_name, "_all_extracted_sequences.fna"))
        writeXStringSet(all_extracted_sequences, output_seq_file, width = 100000)
      })
    }, error = function(e) {
      cat("Error processing reference:", reference_name, "\n", e$message, "\n")
    })
  }, mc.cores = 30)

  cat("Completed mclapply for scenario:", scenario_name, "\n")
}

end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "secs")
elapsed_time_formatted <- sprintf(
  "%02d:%02d:%02d",
  as.integer(elapsed_time) %/% 3600,
  (as.integer(elapsed_time) %% 3600) %/% 60,
  as.integer(elapsed_time) %% 60
)
cat("Total runtime:", elapsed_time_formatted, "\n")
