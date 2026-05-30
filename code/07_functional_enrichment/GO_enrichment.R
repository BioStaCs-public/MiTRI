# ============================================================
# Multi-cancer GO enrichment analysis, RefSeq version (integrated output)
# Based on RefSeq ID to GO ID mapping
# Features: one CSV per cancer type, visualizations combined into a single figure
# ============================================================

library(data.table)
library(clusterProfiler)
library(GO.db)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(cowplot)

# ------------------------------------------------
# 1. Path configuration
# ------------------------------------------------

# Mapping file path (RefSeq -> GO ID)
mapping_file <- "/path/to/project/data_V3/Uniport/mapping_results/protein_GO_mapping.tsv"

# Gene-list directory
gene_list_base_dir <- "/path/to/project/results_V3/summary/sharing_rate_gene/with_speciesINFO"

# Output directory
output_base_dir <- "/path/to/project/results_V3/summary/GO_enrichment"
dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)

# Log files
log_file <- file.path(output_base_dir, "enrichment_log.txt")
progress_file <- file.path(output_base_dir, "progress_checkpoint.rds")

# Define the 12 analysis targets (11 cancer types + PanCancer)
cancer_types <- c("BLCA", "BRCA", "CESC", "CRC", "KIDNEY",
                  "LIHC", "LUNG", "OSCC", "PAAD", "STAD",
                  "THCA", "PanCancer")

# Define gene groups (5 files)
gene_groups <- c("Normal_genes", "Normal_unique_genes",
                 "Shared_genes", "Tumor_genes", "Tumor_unique_genes")

# ------------------------------------------------
# 2. Logging function
# ------------------------------------------------

log_message <- function(msg, file = log_file) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_msg <- paste0("[", timestamp, "] ", msg, "\n")
  cat(full_msg)
  cat(full_msg, file = file, append = TRUE)
}

# ------------------------------------------------
# 3. Read and process mapping file
# ------------------------------------------------

log_message("=" %>% rep(60) %>% paste(collapse = ""))
log_message("Start reading the RefSeq-GO mapping file...")
log_message("=" %>% rep(60) %>% paste(collapse = ""))

start_time <- Sys.time()

# Read only the header first to detect column names.
header <- fread(mapping_file, sep = "\t", nrows = 0, header = TRUE)
log_message("Detected columns:")
log_message(paste("  ", paste(colnames(header), collapse = ", ")))

# Detect RefSeq and GO columns, case-insensitive.
refseq_col <- grep("^refseq$|matched_id", colnames(header), 
                   ignore.case = TRUE, value = TRUE)[1]
go_col <- grep("^go$|go_id|go_term|gene.*ontology", colnames(header), 
               ignore.case = TRUE, value = TRUE)[1]

if (is.na(refseq_col) || is.na(go_col)) {
  log_message("Error: unable to find RefSeq or GO columns.")
  log_message("Please check the mapping file column names.")
  stop("Column detection failed")
}

log_message(sprintf("Using columns: RefSeq='%s', GO='%s'", refseq_col, go_col))

# Read the selected data columns.
mapping_dt <- fread(mapping_file, sep = "\t", header = TRUE, 
                    stringsAsFactors = FALSE,
                    select = c(refseq_col, go_col),
                    showProgress = TRUE)

# Rename to standard column names.
setnames(mapping_dt, old = c(refseq_col, go_col), new = c("RefSeq", "GO"))

log_message(sprintf("Mapping file loaded, elapsed time: %.2f seconds", 
                    difftime(Sys.time(), start_time, units = "secs")))
log_message(sprintf("Total records: %d", nrow(mapping_dt)))
log_message(sprintf("Unique RefSeq IDs: %d", uniqueN(mapping_dt$RefSeq)))

# Remove missing values and invalid GO entries.
mapping_dt <- mapping_dt[!is.na(GO) & GO != "" & !is.na(RefSeq) & RefSeq != ""]
log_message(sprintf("Records after cleanup: %d", nrow(mapping_dt)))

# ------------------------------------------------
# 4. Retrieve GO term information with GO.db
# ------------------------------------------------

log_message("\nExtracting GO term information from GO.db...")

# Extract unique GO IDs.
unique_go_ids <- unique(mapping_dt$GO)
log_message(sprintf("Unique GO IDs: %d", length(unique_go_ids)))

# Query GO information in batches.
get_go_info <- function(go_ids) {
  # Keep only IDs present in GO.db.
  valid_go_ids <- go_ids[go_ids %in% keys(GO.db, keytype = "GOID")]
  
  if (length(valid_go_ids) == 0) {
    return(data.frame(GOID = character(0), 
                      TERM = character(0), 
                      ONTOLOGY = character(0),
                      stringsAsFactors = FALSE))
  }
  
  # Query GO term names.
  terms <- AnnotationDbi::select(GO.db, 
                                 keys = valid_go_ids,
                                 columns = c("TERM", "ONTOLOGY"),
                                 keytype = "GOID")
  
  # Remove duplicate records for the same GO ID.
  terms <- unique(terms)
  
  return(terms)
}

go_info <- get_go_info(unique_go_ids)
log_message(sprintf("Retrieved GO term information for %d GO IDs", nrow(go_info)))

# Merge GO information into the mapping table.
mapping_dt <- merge(mapping_dt, go_info, 
                    by.x = "GO", by.y = "GOID", 
                    all.x = TRUE)

# Remove records without matched GO information.
mapping_dt <- mapping_dt[!is.na(TERM) & !is.na(ONTOLOGY)]
log_message(sprintf("Valid records after merging GO information: %d", nrow(mapping_dt)))

# Summarize the ontology distribution.
ontology_stats <- mapping_dt[, .N, by = ONTOLOGY]
log_message("\nGO ontology distribution:")
for (i in 1:nrow(ontology_stats)) {
  log_message(sprintf("  %s: %d relationships", 
                      ontology_stats$ONTOLOGY[i], 
                      ontology_stats$N[i]))
}

# ------------------------------------------------
# 5. Build TERM2GENE tables by ontology
# ------------------------------------------------

log_message("\nBuilding TERM2GENE data structures...")

build_term2gene <- function(dt, ontology) {
  sub_dt <- dt[ONTOLOGY == ontology, .(GO, RefSeq, TERM)]
  
  # Use GO ID plus term name as TERM, and RefSeq as GENE.
  term2gene <- sub_dt[, .(TERM = paste(GO, TERM, sep = "|"), 
                          GENE = RefSeq)]
  term2gene <- unique(term2gene)
  
  # Build TERM2NAME for later display.
  term2name <- unique(sub_dt[, .(TERM = paste(GO, TERM, sep = "|"), 
                                 NAME = TERM)])
  
  return(list(term2gene = as.data.frame(term2gene),
              term2name = as.data.frame(term2name)))
}

term2gene_list <- list(
  BP = build_term2gene(mapping_dt, "BP"),
  CC = build_term2gene(mapping_dt, "CC"),
  MF = build_term2gene(mapping_dt, "MF")
)

for (ont in c("BP", "CC", "MF")) {
  log_message(sprintf("%s: %d relationships, %d genes", 
                      ont,
                      nrow(term2gene_list[[ont]]$term2gene),
                      uniqueN(term2gene_list[[ont]]$term2gene$GENE)))
}

# Clear memory.
rm(mapping_dt, go_info)
gc()

# ------------------------------------------------
# 6. Fast gene-list reader
# ------------------------------------------------

read_gene_list <- function(file_path) {
  if (!file.exists(file_path)) {
    return(character(0))
  }
  
  genes <- fread(file_path, header = FALSE, stringsAsFactors = FALSE)[[1]]
  genes <- unique(as.character(genes))
  genes <- genes[genes != "" & !is.na(genes)]
  
  return(genes)
}

# ------------------------------------------------
# 7. Enrichment analysis function with tuned parameters
# ------------------------------------------------

run_enrichment_analysis <- function(gene_list, term2gene, term2name, 
                                    background, ontology,
                                    pval_cutoff = 0.05, 
                                    qval_cutoff = 0.2,
                                    min_gs_size = 5,
                                    max_gs_size = 3000) {
  
  if (length(gene_list) == 0) {
    return(NULL)
  }
  
  # Filter TERM2GENE to genes present in the background set.
  term2gene_filtered <- term2gene[term2gene$GENE %in% background, ]
  
  # Check gene annotation coverage.
  genes_with_annotation <- intersect(gene_list, term2gene_filtered$GENE)
  coverage <- length(genes_with_annotation) / length(gene_list) * 100
  
  log_message(sprintf("    Input genes: %d, annotated genes: %d (%.1f%%)",
                      length(gene_list), 
                      length(genes_with_annotation),
                      coverage))
  
  if (length(genes_with_annotation) < 5) {
    log_message("    Warning: fewer than 5 annotated genes, skipping analysis")
    return(NULL)
  }
  
  # Run enrichment analysis.
  result <- tryCatch({
    enricher(
      gene = gene_list,
      TERM2GENE = term2gene_filtered,
      TERM2NAME = term2name,
      universe = background,
      pvalueCutoff = pval_cutoff,
      qvalueCutoff = qval_cutoff,
      pAdjustMethod = "BH",
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size
    )
  }, error = function(e) {
    log_message(sprintf("    Error: %s", e$message))
    return(NULL)
  })
  
  if (!is.null(result) && nrow(result@result) > 0) {
    log_message(sprintf("    Found %d significant enriched terms", nrow(result@result)))
  } else {
    log_message("    No significant enrichment results")
  }
  
  return(result)
}

# ------------------------------------------------
# 8. Integrated visualization function for three ontologies
# ------------------------------------------------

create_integrated_plots <- function(all_results, cancer, output_dir, n_show_per_ont = 10) {
  
  # Collect data for all result-bearing analyses.
  plot_data_list <- list()
  
  for (ontology in c("BP", "CC", "MF")) {
    for (group in gene_groups) {
      
      if (!is.null(all_results[[cancer]][[ontology]][[group]])) {
        # Extract the data component from the stored result object.
        df <- all_results[[cancer]][[ontology]][[group]]$data
        
        if (!is.null(df) && nrow(df) > 0) {
          # Select the top N terms.
          df_top <- head(df[order(df$p.adjust), ], n_show_per_ont)
          df_top$Ontology <- ontology
          df_top$Group <- group
          
          # Convert GeneRatio to a numeric value for sorting and plotting.
          df_top$GeneRatio_val <- sapply(strsplit(as.character(df_top$GeneRatio), "/"), 
                                         function(x) as.numeric(x[1]) / as.numeric(x[2]))
          
          # Shorten long descriptions while keeping non-empty labels.
          df_top$Description_short <- sapply(df_top$Description, function(desc) {
            if (is.na(desc) || desc == "") {
              return("Unknown")
            } else if (nchar(desc) > 50) {
              return(paste0(substr(desc, 1, 47), "..."))
            } else {
              return(desc)
            }
          })
          
          plot_data_list[[paste(ontology, group, sep = "_")]] <- df_top
        }
      }
    }
  }
  
  if (length(plot_data_list) == 0) {
    log_message(sprintf("  %s: no data available for visualization", cancer))
    return()
  }
  
  # Combine all plotting data.
  combined_data <- rbindlist(plot_data_list, fill = TRUE)
  
  # Create integrated plots for each group.
  for (group in unique(combined_data$Group)) {
    
    group_data <- combined_data[combined_data$Group == group, ]
    
    if (nrow(group_data) == 0) next
    
    # Sort by ontology and adjusted p-value, then assign term order within each ontology.
    group_data <- group_data %>%
      group_by(Ontology) %>%
      arrange(p.adjust) %>%
      mutate(term_order = paste0(Ontology, "_", row_number())) %>%
      ungroup() %>%
      arrange(Ontology, p.adjust)
    
    # ========== Dotplot ==========
    tryCatch({
      # Ensure required columns are present and non-empty.
      if (all(c("GeneRatio_val", "Description_short", "Count", "p.adjust") %in% names(group_data)) &&
          nrow(group_data) > 0 &&
          !all(is.na(group_data$Description_short))) {
        
        p1 <- ggplot(group_data, 
                     aes(x = GeneRatio_val, 
                         y = reorder(Description_short, GeneRatio_val),
                         size = Count, 
                         color = p.adjust)) +
          geom_point() +
          scale_color_gradient(low = "red", high = "blue", 
                               name = "Adjusted\nP-value") +
          scale_size_continuous(name = "Gene Count", range = c(2, 8)) +
          facet_wrap(~ Ontology, scales = "free_y", ncol = 1) +
          labs(title = paste(cancer, "-", group, "- GO Enrichment"),
               x = "Gene Ratio",
               y = NULL) +
          theme_bw() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.text.y = element_text(size = 9),
            axis.text.x = element_text(size = 10),
            strip.text = element_text(size = 11, face = "bold"),
            strip.background = element_rect(fill = "lightgray"),
            legend.position = "right",
            panel.grid.minor = element_blank()
          )
        
        out_file <- file.path(output_dir, 
                              paste0(cancer, "_", group, "_integrated_dotplot.pdf"))
        
        # Adjust plot height dynamically.
        n_terms <- nrow(group_data)
        plot_height <- max(10, n_terms * 0.25)
        
        ggsave(out_file, p1, width = 7, height = plot_height, limitsize = FALSE)
        log_message(sprintf("  Saved: %s", basename(out_file)))
      } else {
        log_message("  Dotplot skipped: incomplete data")
      }
      
    }, error = function(e) {
      log_message(sprintf("  Dotplot failed: %s", e$message))
    })
    
    # ========== Barplot ==========
    tryCatch({
      # Ensure required columns are present and non-empty.
      if (all(c("Count", "Description_short", "p.adjust") %in% names(group_data)) &&
          nrow(group_data) > 0 &&
          !all(is.na(group_data$Description_short))) {
        
        p2 <- ggplot(group_data, 
                     aes(x = Count, 
                         y = reorder(Description_short, Count),
                         fill = p.adjust)) +
          geom_bar(stat = "identity") +
          scale_fill_gradient(low = "red", high = "blue", 
                              name = "Adjusted\nP-value") +
          facet_wrap(~ Ontology, scales = "free_y", ncol = 1) +
          labs(title = paste(cancer, "-", group, "- GO Enrichment"),
               x = "Gene Count",
               y = NULL) +
          theme_bw() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.text.y = element_text(size = 9),
            axis.text.x = element_text(size = 10),
            strip.text = element_text(size = 11, face = "bold"),
            strip.background = element_rect(fill = "lightgray"),
            legend.position = "right",
            panel.grid.minor = element_blank()
          )
        
        out_file <- file.path(output_dir, 
                              paste0(cancer, "_", group, "_integrated_barplot.pdf"))
        
        n_terms <- nrow(group_data)
        plot_height <- max(10, n_terms * 0.25)
        
        ggsave(out_file, p2, width = 7, height = plot_height, limitsize = FALSE)
        log_message(sprintf("  Saved: %s", basename(out_file)))
      } else {
        log_message("  Barplot skipped: incomplete data")
      }
      
    }, error = function(e) {
      log_message(sprintf("  Barplot failed: %s", e$message))
    })
  }
}

# ------------------------------------------------
# 9. Main analysis workflow with checkpoint resume support
# ------------------------------------------------

# Read progress.
if (file.exists(progress_file)) {
  progress <- readRDS(progress_file)
  log_message("\nExisting progress detected, resuming from checkpoint...")
  log_message(sprintf("Completed tasks: %d", length(progress$completed)))
} else {
  progress <- list(completed = character(0), 
                   results = list())
}

# Store all results.
all_results <- progress$results

# Iterate over all cancer types.
for (cancer in cancer_types) {
  
  log_message("\n")
  log_message("=" %>% rep(60) %>% paste(collapse = ""))
  log_message(sprintf("Processing: %s", cancer))
  log_message("=" %>% rep(60) %>% paste(collapse = ""))
  
  # Create output directory.
  cancer_output_dir <- file.path(output_base_dir, cancer)
  dir.create(cancer_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Read gene lists.
  cancer_gene_dir <- file.path(gene_list_base_dir, cancer)
  
  if (!dir.exists(cancer_gene_dir)) {
    log_message(sprintf("Warning: directory not found, skipping: %s", cancer))
    next
  }
  
  gene_lists <- list()
  for (group in gene_groups) {
    file_path <- file.path(cancer_gene_dir, paste0(group, ".txt"))
    gene_lists[[group]] <- read_gene_list(file_path)
    log_message(sprintf("  %s: %d genes", group, length(gene_lists[[group]])))
  }
  
  # Build the background set as the union of all genes.
  background_genes <- unique(unlist(gene_lists))
  log_message(sprintf("  Background: %d genes", length(background_genes)))
  
  if (length(background_genes) == 0) {
    log_message("  Warning: no valid genes, skipping this cancer type")
    next
  }
  
  # Initialize result storage for this cancer type.
  if (is.null(all_results[[cancer]])) {
    all_results[[cancer]] <- list()
  }
  
  # Store data frames for integrated output.
  integrated_results <- list()
  
  # Run enrichment analysis for each ontology and gene group.
  for (ontology in c("BP", "CC", "MF")) {
    
    if (is.null(all_results[[cancer]][[ontology]])) {
      all_results[[cancer]][[ontology]] <- list()
    }
    
    for (group in gene_groups) {
      
      task_id <- paste(cancer, ontology, group, sep = "_")
      
      # Skip completed tasks.
      if (task_id %in% progress$completed) {
        log_message(sprintf("\n  Skipping completed task: %s - %s - %s", 
                            cancer, ontology, group))
        next
      }
      
      log_message(sprintf("\n  Analysis: %s - %s - %s", cancer, ontology, group))
      
      gene_list <- gene_lists[[group]]
      
      if (length(gene_list) == 0) {
        log_message("    Gene list is empty, skipping")
        next
      }
      
      # Run enrichment analysis.
      enrich_result <- run_enrichment_analysis(
        gene_list = gene_list,
        term2gene = term2gene_list[[ontology]]$term2gene,
        term2name = term2gene_list[[ontology]]$term2name,
        background = background_genes,
        ontology = ontology
      )
      
      # Save results.
      if (!is.null(enrich_result) && nrow(enrich_result@result) > 0) {
        
        # Use the @result slot directly for stable output.
        result_df <- as.data.frame(enrich_result@result)
        
        # Continue only when the result data frame has rows.
        if (nrow(result_df) > 0) {
          # Add ontology and group columns.
          result_df$Ontology <- ontology
          result_df$Group <- group
          
          # Store both the enrichResult object and the data frame.
          all_results[[cancer]][[ontology]][[group]] <- list(
            enrichResult = enrich_result,
            data = result_df
          )
          integrated_results[[task_id]] <- result_df
          
          # Save the single-analysis RDS file.
          output_prefix <- file.path(cancer_output_dir, 
                                     paste(cancer, ontology, group, sep = "_"))
          saveRDS(enrich_result, paste0(output_prefix, "_enrichResult.rds"))
        } else {
          log_message("    Warning: result data frame is empty after conversion, skipping save")
        }
      }
      
      # Mark the task as completed.
      progress$completed <- c(progress$completed, task_id)
      progress$results <- all_results
      
      # Save progress.
      saveRDS(progress, progress_file)
    }
  }
  
  # ========== Integrated CSV output ==========
  if (length(integrated_results) > 0) {
    
    # Combine all results.
    combined_df <- rbindlist(integrated_results, fill = TRUE)
    
    # Put Ontology and Group at the front.
    col_order <- c("Ontology", "Group", 
                   setdiff(names(combined_df), c("Ontology", "Group")))
    combined_df <- combined_df[, ..col_order]
    
    # Save the integrated CSV.
    integrated_csv_file <- file.path(cancer_output_dir, 
                                     paste0(cancer, "_all_enrichment_results.csv"))
    fwrite(combined_df, integrated_csv_file)
    
    log_message(sprintf("\nOK: %s integrated results saved: %s", cancer, basename(integrated_csv_file)))
    log_message(sprintf("  Total GO terms: %d", nrow(combined_df)))
    
    # Save the integrated enrichResult object list for later batch operations.
    integrated_rds_file <- file.path(cancer_output_dir, 
                                     paste0(cancer, "_all_enrichResults.rds"))
    
    # Extract all enrichResult objects.
    enrich_obj_list <- list()
    for (ont in c("BP", "CC", "MF")) {
      for (grp in gene_groups) {
        if (!is.null(all_results[[cancer]][[ont]][[grp]])) {
          key <- paste(ont, grp, sep = "_")
          enrich_obj_list[[key]] <- all_results[[cancer]][[ont]][[grp]]$enrichResult
        }
      }
    }
    
    if (length(enrich_obj_list) > 0) {
      saveRDS(enrich_obj_list, integrated_rds_file)
      log_message(sprintf("  Integrated RDS object saved: %s", basename(integrated_rds_file)))
      log_message(sprintf("  Includes %d enrichResult objects", length(enrich_obj_list)))
    }
  }
  
  # ========== Generate integrated visualizations ==========
  log_message(sprintf("\nGenerating integrated visualizations for %s...", cancer))
  create_integrated_plots(all_results, cancer, cancer_output_dir, n_show_per_ont = 10)
  
  # ========== Generate per-cancer summary report ==========
  summary_data <- list()
  
  for (ontology in c("BP", "CC", "MF")) {
    for (group in gene_groups) {
      if (!is.null(all_results[[cancer]][[ontology]][[group]])) {
        # Extract the data component.
        df <- all_results[[cancer]][[ontology]][[group]]$data
        if (!is.null(df) && nrow(df) > 0) {
          summary_data[[length(summary_data) + 1]] <- data.frame(
            Cancer = cancer,
            Ontology = ontology,
            Group = group,
            Significant_Terms = nrow(df),
            Top_Pvalue = min(df$pvalue),
            Top_Term = df$Description[1],
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  if (length(summary_data) > 0) {
    summary_df <- rbindlist(summary_data)
    summary_file <- file.path(cancer_output_dir, 
                              paste0(cancer, "_enrichment_summary.csv"))
    fwrite(summary_df, summary_file)
    log_message(sprintf("\n%s summary report saved", cancer))
  }
}

# ------------------------------------------------
# 10. Cross-cancer commonality analysis
# ------------------------------------------------

log_message("\n")
log_message("=" %>% rep(60) %>% paste(collapse = ""))
log_message("Starting cross-cancer commonality analysis...")
log_message("=" %>% rep(60) %>% paste(collapse = ""))

cross_cancer_dir <- file.path(output_base_dir, "Cross_Cancer_Analysis")
dir.create(cross_cancer_dir, recursive = TRUE, showWarnings = FALSE)

for (ontology in c("BP", "CC", "MF")) {
  for (group in gene_groups) {
    
    log_message(sprintf("\nAnalyzing %s - %s", ontology, group))
    
    # Collect enriched terms from all cancer types.
    all_terms <- list()
    
    for (cancer in cancer_types) {
      if (!is.null(all_results[[cancer]][[ontology]][[group]])) {
        # Extract the data component from the stored list.
        result_data <- all_results[[cancer]][[ontology]][[group]]$data
        if (!is.null(result_data) && nrow(result_data) > 0) {
          all_terms[[cancer]] <- result_data$Description
        }
      }
    }
    
    if (length(all_terms) == 0) {
      log_message("  No data available, skipping")
      next
    }
    
    # Calculate term frequencies.
    all_unique_terms <- unique(unlist(all_terms))
    term_freq <- sapply(all_unique_terms, function(term) {
      sum(sapply(all_terms, function(x) term %in% x))
    })
    
    term_freq_df <- data.frame(
      Term = all_unique_terms,
      Frequency = term_freq,
      Percentage = round(term_freq / length(cancer_types) * 100, 1),
      Cancer_List = sapply(all_unique_terms, function(term) {
        cancers <- names(all_terms)[sapply(all_terms, function(x) term %in% x)]
        paste(cancers, collapse = ";")
      }),
      stringsAsFactors = FALSE
    )
    
    term_freq_df <- term_freq_df[order(-term_freq_df$Frequency), ]
    
    # Save results.
    out_file <- file.path(cross_cancer_dir, 
                          paste0("Common_Terms_", ontology, "_", group, ".csv"))
    fwrite(term_freq_df, out_file)
    
    log_message(sprintf("  Unique terms: %d", nrow(term_freq_df)))
    log_message(sprintf("  Present in >=5 cancers: %d", sum(term_freq_df$Frequency >= 5)))
    log_message(sprintf("  Present in >=8 cancers: %d", sum(term_freq_df$Frequency >= 8)))
    
    # Visualize high-frequency terms.
    if (sum(term_freq_df$Frequency >= 5) > 0) {
      high_freq <- term_freq_df[term_freq_df$Frequency >= 5, ]
      
      p <- ggplot(head(high_freq, 30), 
                  aes(x = reorder(Term, Frequency), y = Frequency)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        geom_text(aes(label = Frequency), hjust = -0.2, size = 3) +
        coord_flip() +
        labs(title = paste("High-Frequency Terms:", ontology, "-", group),
             x = "GO Term", y = "Number of Cancers") +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 9),
              plot.title = element_text(hjust = 0.5, face = "bold"))
      
      out_plot <- file.path(cross_cancer_dir, 
                            paste0("Common_Terms_", ontology, "_", 
                                   group, "_barplot.pdf"))
      ggsave(out_plot, p, width = 7, height = max(10, nrow(head(high_freq, 30)) * 0.3))
    }
  }
}

# ------------------------------------------------
# 11. Generate overall summary report
# ------------------------------------------------

log_message("\nGenerating overall summary report...")

all_summaries <- list()

for (cancer in cancer_types) {
  summary_file <- file.path(output_base_dir, cancer, 
                            paste0(cancer, "_enrichment_summary.csv"))
  if (file.exists(summary_file)) {
    df <- fread(summary_file)
    all_summaries[[cancer]] <- df
  }
}

if (length(all_summaries) > 0) {
  final_summary <- rbindlist(all_summaries)
  final_summary_file <- file.path(output_base_dir, 
                                  "All_Cancers_Enrichment_Summary.csv")
  fwrite(final_summary, final_summary_file)
  
  log_message(sprintf("Overall summary report saved: %s", final_summary_file))
  
  # Print summary statistics.
  log_message("\nOverall statistics:")
  summary_stats <- final_summary[, .(
    Total_Terms = sum(Significant_Terms),
    Avg_Terms = round(mean(Significant_Terms), 1),
    N_Cancers = .N
  ), by = .(Ontology, Group)]
  
  print(summary_stats)
}

# Remove progress file after successful completion.
if (file.exists(progress_file)) {
  file.remove(progress_file)
}

log_message("\n")
log_message("=" %>% rep(60) %>% paste(collapse = ""))
log_message("All analyses completed.")
log_message("=" %>% rep(60) %>% paste(collapse = ""))
log_message(sprintf("Output directory: %s", output_base_dir))
log_message(sprintf("Log file: %s", log_file))

log_message("\nOutput file notes:")
log_message("1. {Cancer}_all_enrichment_results.csv - integrated data frame for each cancer type in CSV format")
log_message("2. {Cancer}_all_enrichResults.rds - integrated object list for each cancer type in RDS format")
log_message("3. {Cancer}_{Ontology}_{Group}_enrichResult.rds - object from a single analysis")
log_message("4. {Cancer}_{Group}_integrated_dotplot.pdf - dotplot integrating the three ontologies")
log_message("5. {Cancer}_{Group}_integrated_barplot.pdf - barplot integrating the three ontologies")
log_message("6. {Cancer}_enrichment_summary.csv - summary statistics for each cancer type")
log_message("7. Cross_Cancer_Analysis/ - cross-cancer commonality analysis results")

log_message("\n" %>% rep(60) %>% paste(collapse = "="))
log_message("Example code for reading the output data")
log_message("=" %>% rep(60) %>% paste(collapse = ""))

log_message("\n# Method 1: read integrated CSV, useful for quick inspection and filtering")
log_message("library(data.table)")
log_message("blca_results <- fread('BLCA/BLCA_all_enrichment_results.csv')")
log_message("tumor_bp <- blca_results[Ontology == 'BP' & Group == 'Tumor_genes']")

log_message("\n# Method 2: read integrated RDS, useful for advanced analysis and visualization")
log_message("library(clusterProfiler)")
log_message("blca_enrich <- readRDS('BLCA/BLCA_all_enrichResults.rds')")
log_message("# List structure: blca_enrich[['BP_Tumor_genes']] returns an enrichResult object")
log_message("# Use clusterProfiler functions directly:")
log_message("dotplot(blca_enrich[['BP_Tumor_genes']], showCategory=20)")

log_message("\n# Method 3: read a single RDS, most flexible for one analysis")
log_message("tumor_bp_obj <- readRDS('BLCA/BLCA_BP_Tumor_genes_enrichResult.rds')")
log_message("dotplot(tumor_bp_obj)")
log_message("enrichmap(tumor_bp_obj)")
log_message("cnetplot(tumor_bp_obj)")

log_message("\n# Method 4: batch-read one analysis across all cancer types")
log_message("all_tumor_bp <- list()")
log_message("for (cancer in cancer_types) {")
log_message("  file <- paste0(cancer, '/', cancer, '_BP_Tumor_genes_enrichResult.rds')")
log_message("  if (file.exists(file)) {")
log_message("    all_tumor_bp[[cancer]] <- readRDS(file)")
log_message("  }")
log_message("}")
