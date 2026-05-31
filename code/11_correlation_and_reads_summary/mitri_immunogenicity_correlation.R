# =========================
# MiTRI vs Immunogenicity Top-k Sensitivity
# Outputs:
#  - CosineSimilarity_MiTRI_vs_Immuno_Topk_facet.(pdf/svg)
#  - Spearman_MiTRI_vs_Immuno_Topk_facet.(pdf/svg)
#  - Spearman_overall_summary_Topk_byTool.csv
#  - QC_sparse_groups_npair_lt2.csv
# =========================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(proxy)
  library(ggplot2)
})

# ====== Input files (5: Top1/3/5/10/All immunogenicity aggregations) ======
in_files <- c(
  "/path/to/project/results_V3/summary/tumor_specific_MTR/tumor_specific_MiTRI_with_normalized_gene_diversity_peptides_immunoTop1.csv",
  "/path/to/project/results_V3/summary/tumor_specific_MTR/tumor_specific_MiTRI_with_normalized_gene_diversity_peptides_immunoTop3.csv",
  "/path/to/project/results_V3/summary/tumor_specific_MTR/tumor_specific_MiTRI_with_normalized_gene_diversity_peptides_immunoTop5.csv",
  "/path/to/project/results_V3/summary/tumor_specific_MTR/tumor_specific_MiTRI_with_normalized_gene_diversity_peptides_immunoTop10.csv",
  "/path/to/project/results_V3/summary/tumor_specific_MTR/tumor_specific_MiTRI_with_normalized_gene_diversity_peptides_immunoAll.csv"
)

# ====== Output directory ======
save_path <- "/path/to/project/results_V3/summary/figure3_cor/imm"
dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

# ====== Safe Spearman: return NA when fewer than 2 pairs to avoid cor() errors ======
safe_spearman <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 2) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

# ====== Safe p-value (cor.test is unstable for very small samples; require >=3) ======
safe_spearman_p <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) return(NA_real_)
  tryCatch(cor.test(x[ok], y[ok], method = "spearman")$p.value,
           error = function(e) NA_real_)
}

# ====== Cosine similarity (guard against NA, which would break dist) ======
safe_cosine <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 1) return(NA_real_)
  m <- rbind(x[ok], y[ok])
  d <- proxy::dist(m, method = "cosine")
  sim <- 1 - as.numeric(d)
  return(sim)
}

# ====== Read one file and pivot to long format: K x Tool x (Sample, Ref, ...) ======
read_one <- function(fp) {
  dat <- read.csv(fp) %>%
    mutate(
      Reference = gsub("_", " ", Reference),
      K = dplyr::case_when(
        grepl("Top1\\.csv$", fp) ~ "Top1",
        grepl("Top3\\.csv$", fp) ~ "Top3",
        grepl("Top5\\.csv$", fp) ~ "Top5",
        grepl("Top10\\.csv$", fp) ~ "Top10",
        grepl("All\\.csv$", fp) ~ "All",
        TRUE ~ "Unknown"
      )
    )

  # Match immunogenicity tool score columns: IEDB_top5mean / IEDB_allmean, etc.
  score_cols <- grep("^(IEDB|DeepImmuno|BigMHC)_(top\\d+mean|allmean)$", names(dat), value = TRUE)
  if (length(score_cols) == 0) {
    stop(paste0("No immunogenicity score columns found in file: ", fp))
  }

  dat %>%
    pivot_longer(
      cols = all_of(score_cols),
      names_to = "ToolRaw",
      values_to = "Score"
    ) %>%
    mutate(
      Tool = dplyr::case_when(
        grepl("^IEDB_", ToolRaw) ~ "IEDB",
        grepl("^DeepImmuno_", ToolRaw) ~ "DeepImmuno",
        grepl("^BigMHC_", ToolRaw) ~ "BigMHC",
        TRUE ~ ToolRaw
      )
    ) %>%
    select(-ToolRaw)
}

all_data <- bind_rows(lapply(in_files, read_one))

# ====== Harmonise factor levels so facet ordering is consistent ======
all_data$Cancer_Type <- factor(all_data$Cancer_Type, levels = sort(unique(all_data$Cancer_Type)))
all_data$Reference <- factor(all_data$Reference, levels = sort(unique(all_data$Reference), decreasing = TRUE))
all_data$K <- factor(all_data$K, levels = c("Top1", "Top3", "Top5", "Top10", "All"))
all_data$Tool <- factor(all_data$Tool, levels = c("IEDB", "DeepImmuno", "BigMHC"))

# =========================
# 1) Grouped cosine-similarity heatmap (facet)
# =========================
cosine_df <- all_data %>%
  group_by(K, Tool, Cancer_Type, Reference) %>%
  summarise(
    n_pair = sum(complete.cases(Activity, Score)),
    cosine_similarity = safe_cosine(Activity, Score),
    .groups = "drop"
  )

plot_cosine_facet <- function(df, title) {
  ggplot(df, aes(x = Cancer_Type, y = Reference, fill = cosine_similarity)) +
    geom_tile(width = 0.8, height = 0.8) +
    geom_text(
      aes(label = ifelse(is.na(cosine_similarity), "", sprintf("%.2f", cosine_similarity))),
      size = 2.2, color = "black", family = "Helvetica", fontface = "bold"
    ) +
    scale_fill_gradient2(
      low = "#2a749a", mid = "#d9d3c8", high = "#9d001f",
      midpoint = 0.5, limits = c(0.0, 1.0), na.value = "white"
    ) +
    facet_grid(Tool ~ K, scales = "fixed") +
    theme_minimal() +
    labs(title = title, x = "Cancer Type", y = "Reference") +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 9, family = "Helvetica", face = "bold", color = "black"),
      legend.key.size = unit(0.35, "cm"),
      legend.spacing.x = unit(0.35, "cm"),
      legend.margin = margin(t = -4, b = -4),
      plot.title = element_text(size = 10, family = "Helvetica", face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 6.8, family = "Helvetica", face = "bold", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6.8, family = "Helvetica", face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_line(color = "black", linewidth = 0.4),
      axis.ticks.y = element_line(color = "black", linewidth = 0.4),
      axis.line = element_line(linewidth = 0.4, color = "black"),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      strip.text = element_text(size = 9, face = "bold", family = "Helvetica")
    )
}

p_cos <- plot_cosine_facet(
  cosine_df,
  "Cosine Similarity: MiTRI Activity vs Immunogenicity (Top-k Sensitivity)"
)

ggsave(file.path(save_path, "CosineSimilarity_MiTRI_vs_Immuno_Topk_facet.pdf"),
       plot = p_cos, width = 14, height = 9, device = "pdf")
ggsave(file.path(save_path, "CosineSimilarity_MiTRI_vs_Immuno_Topk_facet.svg"),
       plot = p_cos, width = 14, height = 9, device = "svg")

# =========================
# 2) Grouped Spearman-correlation heatmap (facet)
# =========================
spearman_group_df <- all_data %>%
  group_by(K, Tool, Cancer_Type, Reference) %>%
  summarise(
    n_pair = sum(complete.cases(Activity, Score)),
    spearman = safe_spearman(Activity, Score),
    .groups = "drop"
  )

plot_spearman_facet <- function(df, title) {
  ggplot(df, aes(x = Cancer_Type, y = Reference, fill = spearman)) +
    geom_tile(width = 0.8, height = 0.8) +
    geom_text(
      aes(label = ifelse(is.na(spearman), "", sprintf("%.2f", spearman))),
      size = 2.2, color = "black", family = "Helvetica", fontface = "bold"
    ) +
    scale_fill_gradient2(
      low = "#2a749a", mid = "#d9d3c8", high = "#9d001f",
      midpoint = 0, limits = c(-1, 1), na.value = "white"
    ) +
    facet_grid(Tool ~ K, scales = "fixed") +
    theme_minimal() +
    labs(title = title, x = "Cancer Type", y = "Reference") +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 9, family = "Helvetica", face = "bold", color = "black"),
      legend.key.size = unit(0.35, "cm"),
      legend.spacing.x = unit(0.35, "cm"),
      legend.margin = margin(t = -4, b = -4),
      plot.title = element_text(size = 10, family = "Helvetica", face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 6.8, family = "Helvetica", face = "bold", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6.8, family = "Helvetica", face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_line(color = "black", linewidth = 0.4),
      axis.ticks.y = element_line(color = "black", linewidth = 0.4),
      axis.line = element_line(linewidth = 0.4, color = "black"),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      strip.text = element_text(size = 9, face = "bold", family = "Helvetica")
    )
}

p_spear <- plot_spearman_facet(
  spearman_group_df,
  "Spearman Correlation: MiTRI Activity vs Immunogenicity (Top-k Sensitivity)"
)

ggsave(file.path(save_path, "Spearman_MiTRI_vs_Immuno_Topk_facet.pdf"),
       plot = p_spear, width = 16, height = 12, device = "pdf")
ggsave(file.path(save_path, "Spearman_MiTRI_vs_Immuno_Topk_facet.svg"),
       plot = p_spear, width = 16, height = 12, device = "svg")

# =========================
# 3) Overall Spearman summary (one row per K x Tool)
# =========================
spearman_overall_df <- all_data %>%
  group_by(K, Tool) %>%
  summarise(
    n = sum(complete.cases(Activity, Score)),
    rho = safe_spearman(Activity, Score),
    p_value = safe_spearman_p(Activity, Score),
    .groups = "drop"
  )

write.csv(spearman_overall_df,
          file.path(save_path, "Spearman_overall_summary_Topk_byTool.csv"),
          row.names = FALSE)

# =========================
# 4) QC output: which groups have n_pair < 2 (causing NA Spearman)
# =========================
qc_sparse <- spearman_group_df %>%
  filter(is.na(spearman) | n_pair < 2)

write.csv(qc_sparse,
          file.path(save_path, "QC_sparse_groups_npair_lt2.csv"),
          row.names = FALSE)

cat("Done. Output directory:\n")
cat(save_path, "\n\n")
cat("Main outputs:\n")
cat("- CosineSimilarity_MiTRI_vs_Immuno_Topk_facet.(pdf/svg)\n")
cat("- Spearman_MiTRI_vs_Immuno_Topk_facet.(pdf/svg)\n")
cat("- Spearman_overall_summary_Topk_byTool.csv\n")
cat("- QC_sparse_groups_npair_lt2.csv\n")
