library(ggplot2)

# ============ Data paths ============
paths <- list(
  species = list(
    tumor = "/path/to/project/results_V3/summary/species_similarity_analysis/with_speciesINFO/tumor_cosine_similarities.csv",
    normal = "/path/to/project/results_V3/summary/species_similarity_analysis/with_speciesINFO/normal_cosine_similarities.csv"
  ),
  gene = list(
    tumor = "/path/to/project/results_V3/summary/gene_similarity_analysis/with_speciesINFO/tumor_cosine_similarities.csv",
    normal = "/path/to/project/results_V3/summary/gene_similarity_analysis/with_speciesINFO/normal_cosine_similarities.csv"
  ),
  peptide = list(
    tumor = "/path/to/project/results_V3/summary/peptide_similarity_analysis/with_speciesINFO/tumor_cosine_similarities.csv",
    normal = "/path/to/project/results_V3/summary/peptide_similarity_analysis/with_speciesINFO/normal_cosine_similarities.csv"
  )
)

output_dir <- "/path/to/project/results_V3/summary/heatmap/finnal"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============ Data processing ============
replace_names <- function(df) {
  df$Cancer_Type1 <- gsub("KIDNEY", "RCC", df$Cancer_Type1)
  df$Cancer_Type2 <- gsub("KIDNEY", "RCC", df$Cancer_Type2)
  df$Cancer_Type1 <- gsub("LUNG", "LC", df$Cancer_Type1)
  df$Cancer_Type2 <- gsub("LUNG", "LC", df$Cancer_Type2)
  df$Cancer_Type1 <- gsub("STAD", "GC", df$Cancer_Type1)
  df$Cancer_Type2 <- gsub("STAD", "GC", df$Cancer_Type2)
  df
}

stats_data <- data.frame()
all_diffs <- list()
all_keys <- list()

for (level in names(paths)) {
  df_tumor  <- replace_names(read.csv(paths[[level]]$tumor))
  df_normal <- replace_names(read.csv(paths[[level]]$normal))

  df_tumor  <- df_tumor[df_tumor$Cancer_Type1 != df_tumor$Cancer_Type2, ]
  df_normal <- df_normal[df_normal$Cancer_Type1 != df_normal$Cancer_Type2, ]


  # ---- Alignment check: tumor vs normal (within each level) ----
  key_t <- paste(df_tumor$Cancer_Type1, df_tumor$Cancer_Type2, sep = "||")
  key_n <- paste(df_normal$Cancer_Type1, df_normal$Cancer_Type2, sep = "||")

  all_keys[[level]] <- key_t

  # 1) same number of pairs|
  if (length(key_t) != length(key_n)) {
    message("[", level, "] tumor vs normal: different number of pairs! tumor=",
            length(key_t), " normal=", length(key_n))
  }

  # 2) same pair set, regardless of order|
  only_in_t <- setdiff(key_t, key_n)
  only_in_n <- setdiff(key_n, key_t)
  if (length(only_in_t) > 0 || length(only_in_n) > 0) {
    message("[", level, "] tumor vs normal: pair set mismatch!")
    message("  only in tumor: ", min(5, length(only_in_t)))
    if (length(only_in_t) > 0) print(head(only_in_t, 5))
    message("  only in normal: ", min(5, length(only_in_n)))
    if (length(only_in_n) > 0) print(head(only_in_n, 5))
  }

  # 3) same row order| (decides whether paired_diff is valid)
  order_match <- identical(key_t, key_n)
  if (!order_match) {
    message("[", level, "] tumor vs normal: NOT aligned by row order (paired subtraction is unsafe).")
    # show the first few mismatched positions
    idx <- which(key_t != key_n)
    message("  first mismatched indices:"); print(head(idx, 10))
    message("  tumor keys at mismatches:");  print(head(key_t[idx], 5))
    message("  normal keys at mismatches:"); print(head(key_n[idx], 5))
  } else {
    message("[", level, "] tumor vs normal: aligned by row order OK.")
  }


  paired_diff <- df_tumor$Cosine_Similarity - df_normal$Cosine_Similarity
  all_diffs[[level]] <- paired_diff

  # stats_data <- rbind(stats_data, data.frame(
  #   Level = level,
  #   Difference = mean(paired_diff),
  #   Diff_SE = sd(paired_diff) / sqrt(length(paired_diff))
  # ))

  n <- length(paired_diff)
  se <- sd(paired_diff) / sqrt(n)
  t_crit <- qt(0.975, df = n - 1)

  stats_data <- rbind(stats_data, data.frame(
    Level = level,
    Difference = mean(paired_diff),
    CI_lower = mean(paired_diff) - t_crit * se,
    CI_upper = mean(paired_diff) + t_crit * se
  ))

}

message("[Across layers] Species vs Gene order identical: ", identical(all_keys$species, all_keys$gene))
message("[Across layers] Gene vs Peptide order identical: ", identical(all_keys$gene, all_keys$peptide))
message("[Across layers] Species vs Peptide order identical: ", identical(all_keys$species, all_keys$peptide))

message("[Across layers] Species vs Gene set identical: ", setequal(all_keys$species, all_keys$gene))
message("[Across layers] Gene vs Peptide set identical: ", setequal(all_keys$gene, all_keys$peptide))
message("[Across layers] Species vs Peptide set identical: ", setequal(all_keys$species, all_keys$peptide))


stats_data$Level <- factor(stats_data$Level, levels = c("species", "gene", "peptide"))
stats_data$Level_Label <- c("Species", "Gene", "Peptide")
stats_data$x_pos <- c(0.8, 1.3, 1.8)

# ============ Statistical tests ============
diff_matrix <- cbind(
  Species  = all_diffs$species,
  Gene     = all_diffs$gene,
  Peptide  = all_diffs$peptide
)

friedman_result <- friedman.test(diff_matrix)

chi2_val <- round(friedman_result$statistic, 2)
p_val <- friedman_result$p.value
p_txt <- if (p_val < 0.001) {
  formatC(p_val, format = "e", digits = 1)
} else {
  print(p_val)
  round(p_val, 3)
}

pairwise_p <- data.frame(
  Comparison = c("Species vs Gene", "Gene vs Peptide", "Species vs Peptide"),
  P_value = c(
    wilcox.test(all_diffs$species, all_diffs$gene, paired = TRUE)$p.value,
    wilcox.test(all_diffs$gene, all_diffs$peptide, paired = TRUE)$p.value,
    wilcox.test(all_diffs$species, all_diffs$peptide, paired = TRUE)$p.value
  )
)

pairwise_p$Significance <- ifelse(pairwise_p$P_value < 1e-4, "****",
                                  ifelse(pairwise_p$P_value < 1e-3, "***",
                                         ifelse(pairwise_p$P_value < 1e-2, "**",
                                                ifelse(pairwise_p$P_value < 0.05, "*", "ns"))))

# ============ Plotting ============
p3 <- ggplot(stats_data, aes(x = x_pos, y = Difference)) +

  geom_hline(yintercept = 0, color = "gray30", linewidth = 0.8) +

  # trend line
  geom_line(linewidth = 0.8, color = "gray40", linetype = "solid") +

  geom_segment(aes(xend = x_pos, yend = Difference,
                   color = Difference > 0),
               y = 0, linewidth = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("#4e92b9", "#E74C3C"), guide = "none") +

  geom_errorbar(aes(min = CI_lower,
                    ymax = CI_upper),
                width = 0.10, linewidth = 0.8) +
  geom_point(size = 2.8, shape = 21, fill = "white", stroke = 1.2) +

  # Species vs Gene
  annotate("segment", x = 0.8, xend = 1.3, y = 0.085, yend = 0.085) +
  annotate("segment", x = 0.8, xend = 0.8, y = 0.085, yend = 0.082) +
  annotate("segment", x = 1.3, xend = 1.3, y = 0.085, yend = 0.082) +
  annotate("text", x = 1.05, y = 0.086,
           label = pairwise_p$Significance[1],
           family = "Helvetica", fontface = "bold", size = 3.5) +

  # Gene vs Peptide
  annotate("segment", x = 1.3, xend = 1.8, y = 0.101, yend = 0.101) +
  annotate("segment", x = 1.3, xend = 1.3, y = 0.101, yend = 0.098) +
  annotate("segment", x = 1.8, xend = 1.8, y = 0.101, yend = 0.098) +
  annotate("text", x = 1.55, y = 0.103,
           label = pairwise_p$Significance[2],
           family = "Helvetica", fontface = "bold", size = 3.5) +

  # Species vs Peptide
  annotate("segment", x = 0.8, xend = 1.8, y = 0.112, yend = 0.112) +
  annotate("segment", x = 0.8, xend = 0.8, y = 0.112, yend = 0.109) +
  annotate("segment", x = 1.8, xend = 1.8, y = 0.112, yend = 0.109) +
  annotate("text", x = 1.3, y = 0.114,
           label = pairwise_p$Significance[3],
           family = "Helvetica", fontface = "bold", size = 3.5) +

  scale_x_continuous(breaks = stats_data$x_pos,
                     labels = stats_data$Level_Label,
                     limits = c(0.3, 2.3)) +

  scale_y_continuous(breaks = seq(-0.10, 0.15, by = 0.05)) +

  coord_cartesian(ylim = c(-0.10, 0.12), clip = "off") +

  # Friedman test result annotation
  annotate("text", x = 0.9, y = -0.06, hjust = 0,
           label = paste0("p = ", p_txt),
           family = "Helvetica", size = 2.5, fontface = "bold") +

  theme_classic() +
  labs(
    x = "",
    y = expression(Delta ~ "Cosine similarity (Tumor - Normal)")
  ) +
  theme(
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    axis.line.x = element_blank(),
    axis.ticks.y = element_line(color = "black", linewidth = 0.8),
    axis.ticks.x = element_blank(),
    axis.text = element_text(face = "bold", family = "Helvetica", size = 8, color = "black"),
    axis.title.y = element_text(face = "bold", family = "Helvetica", size = 8),
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave(file.path(output_dir, "cross_level_comparison.pdf"),
       plot = p3, width = 2.5, height = 2.2, device = cairo_pdf)
ggsave(file.path(output_dir, "cross_level_comparison.svg"),
       plot = p3, width = 2.5, height = 2.2)
