# ============================================================================
# Bar plot of the number of unique peptides per cancer type
# ============================================================================

library(ggplot2)
library(dplyr)

# ============================================================================
# 1. Read data
# ============================================================================

peptide_type <- "8cancers_Tumor_peptides"
rds_file <- sprintf("/path/to/project/results_V3/summary/unpset_peptides/epitope_presence_matrix_%s.rds",
                    peptide_type)
epitope_matrix <- readRDS(rds_file)
output_dir <- "/path/to/project/results_V3/summary/unpset_peptides/figures"
# ============================================================================
# 2. Count peptides per cancer type
# ============================================================================

cancer_cols <- colnames(epitope_matrix)[colnames(epitope_matrix) != "Epitope"]

peptide_counts <- data.frame(
  Cancer = cancer_cols,
  Count = colSums(epitope_matrix[, cancer_cols]),
  stringsAsFactors = FALSE  # important: prevent automatic conversion to factor
)

# # sort by count
# peptide_counts <- peptide_counts %>%
#   arrange(desc(Count))

print(peptide_counts)


# ============================================================================
# 4. Bar plot (ggplot2)
# ============================================================================

# ensure Cancer is character, then set factor order explicitly
peptide_counts$Cancer <- as.character(peptide_counts$Cancer)
peptide_counts$Cancer <- factor(peptide_counts$Cancer,
                                levels = peptide_counts$Cancer)

# inspect data
print(str(peptide_counts))

# add a numeric position to control spacing
peptide_counts$Position <- seq(1, by = 0.5, length.out = nrow(peptide_counts))  # spacing factor 0.5; smaller = tighter

# plot
p <- ggplot(peptide_counts, aes(x = Position, y = Count / 1e6)) +  # divide by one million
  geom_bar(stat = "identity", fill = "#595959", width = 0.2) +  # keep original bar width
  scale_x_continuous(
    breaks = peptide_counts$Position,
    labels = peptide_counts$Cancer,
    expand = c(0.02, 0.02)
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "M")  # add "M" suffix
  ) +
  labs(
    title = paste(""),
    x = "",
    y = "Number of unique peptides (million)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 8, family="Helvetica", face = "bold", angle = 270, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(size = 8, family="Helvetica", face = "bold", color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 9, family="Helvetica", face = "bold"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.line.x = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    plot.margin = margin(t = 5, r = 10, b = 10, l = 5, unit = "pt")
  )
# display plot (in RStudio)
tryCatch({
  print(p)
}, error = function(e) {
  cat("Error displaying plot, but saving should be fine: ", e$message, "\n")
})

# ============================================================================
# 5. Save figures
# ============================================================================

# Save as PDF
ggsave(
  filename = file.path(output_dir, sprintf("peptide_counts_%s.pdf", peptide_type)),
  plot = p,
  width = 2,
  height = 1.5
)

# Save as SVG
ggsave(
  filename = file.path(output_dir, sprintf("peptide_counts_%s.svg", peptide_type)),
  plot = p,
  width = 2,
  height = 1.5
)

cat("ggplot2 figure saved as PDF and SVG\n")


tryCatch({
  print(p_horizontal)
}, error = function(e) {
  cat("Error displaying horizontal plot, but saving should be fine: ", e$message, "\n")
})

cat("Horizontal bar plot saved\n")
cat("All figures saved.\n")
