# ---- Dependencies ----
library(dplyr)
library(stringr)
library(Biostrings)
library(data.table)

# ---- Config (anonymized) ----
PROJECT_ROOT <- "/project"
PROJECT_TYPE <- "PROJECT_TYPE"
MICROBE_CSV  <- file.path(PROJECT_ROOT, "data/csvs_special", "project_reference.csv")

COVERAGE_DIR_ROOT <- file.path(PROJECT_ROOT, "results/special", PROJECT_TYPE, "06.coverage_target", "T")
REGION_DIR_ROOT   <- file.path(PROJECT_ROOT, "results/special", PROJECT_TYPE, "04.activity_region_T_Nback_share", "04.2.region")
OUTPUT_DIR_ROOT   <- file.path(PROJECT_ROOT, "results/special", PROJECT_TYPE, "04.activity_region_reads_T_Nback_share")

# ---- Input list of target microbes ----
microbe_list <- read.csv(MICROBE_CSV)$Taxon

# ---- Parse .fna headers to regions (expects headers like: <refA>_<refB>_<start>_<end>) ----
parse_fna <- function(fna_file) {
  sequences <- readDNAStringSet(fna_file)
  headers <- names(sequences)
  parts <- data.table::tstrsplit(headers, "_", fixed = TRUE)

  regions <- data.frame(
    Reference = paste(parts[[1]], parts[[2]], sep = "_"),
    Start     = as.integer(parts[[3]]),
    End       = as.integer(parts[[4]]),
    stringsAsFactors = FALSE
  )
  regions
}

# ---- Main loop ----
for (microbe in microbe_list) {
  coverage_dir <- file.path(COVERAGE_DIR_ROOT, microbe)
  fna_dir      <- file.path(REGION_DIR_ROOT,   microbe)
  output_dir   <- file.path(OUTPUT_DIR_ROOT,   microbe)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  bam_txt_files <- list.files(coverage_dir, pattern = "\\.bam\\.txt$", full.names = TRUE)
  bam_txt_files <- bam_txt_files[!grepl("tmp", bam_txt_files)]

  fna_files <- list.files(fna_dir, pattern = "\\.fna$", full.names = TRUE)

  # Filter empty .fna files
  if (length(fna_files)) {
    finfo <- file.info(fna_files)
    fna_files <- fna_files[!is.na(finfo$size) & finfo$size > 0]
  }
  if (length(fna_files) == 0L) {
    message("[", microbe, "] no usable .fna files; skipping.")
    next
  }

  # Match samples between fna and bam.txt
  fna_samples <- sub("_.*", "", basename(fna_files))
  bam_samples <- sub("_.*", "", basename(bam_txt_files))
  bam_txt_files <- bam_txt_files[bam_samples %in% fna_samples]

  if (length(bam_txt_files) == 0L) {
    message("[", microbe, "] no bam txt files matched to fna samples; skipping.")
    next
  }

  for (fna_file in fna_files) {
    fna_sample <- sub("_.*", "", basename(fna_file))

    bam_samples <- sub("_.*", "", basename(bam_txt_files))
    cand <- bam_txt_files[bam_samples == fna_sample]

    if (length(cand) == 0L) {
      message("[", microbe, "][", fna_sample, "] no matched bam txt; skipping.")
      next
    }
    if (length(cand) > 1L) {
      message("[", microbe, "][", fna_sample, "] multiple bam txt matched; using the first.")
    }
    bam_txt_file <- cand[1]

    # Parse regions from .fna (must contain Reference/Start/End)
    regions <- tryCatch(parse_fna(fna_file), error = function(e) NULL)
    if (is.null(regions) || nrow(regions) == 0L) {
      message("[", microbe, "][", fna_sample, "] regions empty/failed; skipping.")
      next
    }
    if (!all(c("Reference","Start","End") %in% names(regions))) {
      message("[", microbe, "][", fna_sample, "] regions missing Reference/Start/End; skipping.")
      next
    }

    # Read bam txt (TAB-separated)
    bam_data <- tryCatch(
      read.table(bam_txt_file, sep = "\t", header = FALSE,
                 comment.char = "@", fill = TRUE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(bam_data) || nrow(bam_data) == 0L) {
      message("[", microbe, "][", fna_sample, "] empty/failed bam data; skipping.")
      next
    }

    # Keep first up to 10 columns; DO NOT pad if fewer than 10 (per your request)
    nkeep <- min(10L, ncol(bam_data))
    bam_data <- bam_data[, seq_len(nkeep), drop = FALSE]
    colnames(bam_data) <- c("QNAME","FLAG","RNAME","POS","MAPQ","CIGAR","RNEXT","PNEXT","TLEN","SEQ")[seq_len(nkeep)]

    # Drop rows with missing POS
    # (use qualified namespace to avoid loading tidyr)
    bam_data <- bam_data %>% tidyr::drop_na(POS)
    if (nrow(bam_data) == 0L) {
      message("[", microbe, "][", fna_sample, "] no valid bam rows after POS filtering; skipping.")
      next
    }

    # Largest 'M' in CIGAR (numeric before 'M')
    bam_data <- bam_data %>%
      mutate(
        POS   = suppressWarnings(as.integer(POS)),
        CIGAR = as.character(CIGAR),
        SEQ   = as.character(if ("SEQ" %in% names(.)) SEQ else NA_character_),
        maxM  = suppressWarnings(
          as.numeric(str_remove(str_extract(CIGAR, "\\d+M"), "M"))
        )
      )

    # Determine read length: prefer maxM; fallback to max nchar(SEQ) excluding "*"
    RL1 <- suppressWarnings(max(bam_data$maxM, na.rm = TRUE))
    if (is.infinite(RL1)) RL1 <- NA_real_
    if (is.na(RL1) || RL1 <= 0) {
      seq_len_vec <- if ("SEQ" %in% names(bam_data))
        nchar(bam_data$SEQ[!is.na(bam_data$SEQ) & bam_data$SEQ != "*"]) else integer(0)
      ReadLength <- if (length(seq_len_vec)) max(seq_len_vec, na.rm = TRUE) else NA_real_
    } else {
      ReadLength <- RL1
    }

    if (is.na(ReadLength) || ReadLength <= 0) {
      message("[", microbe, "][", fna_sample, "] cannot determine ReadLength; skipping.")
      next
    }

    # Convert to data.table
    bam_dt <- as.data.table(bam_data)
    regions_dt <- as.data.table(regions)

    # Sanitize region coordinates
    regions_dt[, `:=`(
      Start = suppressWarnings(as.integer(Start)),
      End   = suppressWarnings(as.integer(End))
    )]
    regions_dt <- regions_dt[!is.na(Start) & !is.na(End) & End >= Start]
    if (nrow(regions_dt) == 0L) {
      message("[", microbe, "][", fna_sample, "] no valid regions; skipping.")
      next
    }

    # Remove reads with NA POS and compute End_POS
    bam_dt <- bam_dt[!is.na(POS)]
    if (nrow(bam_dt) == 0L) {
      message("[", microbe, "][", fna_sample, "] all POS are NA; skipping.")
      next
    }
    bam_dt[, End_POS := POS + as.integer(ReadLength) - 1L]

    # Coarse filter by global bounds per Reference
    regions_bounds <- regions_dt[, .(minS = min(Start, na.rm = TRUE),
                                     maxE = max(End,   na.rm = TRUE)),
                                 by = Reference]

    bam_dt <- merge(
      bam_dt, regions_bounds,
      by.x = "RNAME", by.y = "Reference",
      all.x = FALSE, all.y = FALSE
    )
    if (nrow(bam_dt) == 0L) {
      message("[", microbe, "][", fna_sample, "] reference mismatch or no candidate reads; skipping.")
      next
    }

    bam_dt <- bam_dt[End_POS >= minS & POS <= maxE]
    if (nrow(bam_dt) == 0L) {
      message("[", microbe, "][", fna_sample, "] no candidates after coarse filter; skipping.")
      next
    }

    # Fine overlap filter (allow cartesian)
    setkey(regions_dt, Reference)
    setkey(bam_dt, RNAME)

    res <- bam_dt[regions_dt, allow.cartesian = TRUE][
      (POS >= Start & POS <= End) |
        (End_POS >= Start & End_POS <= End) |
        (POS <= Start & End_POS >= End),
      .(QNAME, RNAME, POS, SEQ, maxM)
    ]

    # Write result
    if (!is.null(res) && nrow(res) > 0L) {
      out_name <- paste0(basename(bam_txt_file), "_filtered.csv")
      output_file <- file.path(output_dir, out_name)
      fwrite(res, output_file)
      message("[", microbe, "][", fna_sample, "] wrote: ", output_file, " (", nrow(res), " rows)")
    } else {
      message("[", microbe, "][", fna_sample, "] no overlapping reads; no file written.")
    }
  }
}
