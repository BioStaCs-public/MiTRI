library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# ===== Paths =====
input_dir <- "/path/to/project/results_V3/special/CRC/01.group_matrix"
output_dir <- "/path/to/project/results_V3/special/CRC/01.group_matrix/figures"

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ===== Target microbes =====
target_microbes <- c(
  'Fusobacterium_animalis',
  'Fusobacterium_canifelinum',
  'Fusobacterium_hwasookii',
  'Fusobacterium_massiliense',
  'Fusobacterium_nucleatum',
  'Fusobacterium_periodonticum',
  'Fusobacterium_polymorphum',
  'Fusobacterium_pseudoperiodonticum',
  'Fusobacterium_vincentii'
)

# ===== Sample groups =====
sample_groups <- list(
  R = c('B1', 'B4', 'B7', 'E5'),
  NR = c('B2', 'B3', 'B6', 'C6')
)

# ===== Colour scheme =====
colors <- list(
  R_fill = "#d64863",
  NR_fill = "#4E92B9",
  R_edge = "#7e3d4e",
  NR_edge = "#20364F"
)

# ===== Read data =====
cat("============================================================\n")
cat("Boxplots of Fusobacterium species\n")
cat("============================================================\n\n")

r_file <- file.path(input_dir, "CRC_R_species_pathseq_Normalized.csv")
nr_file <- file.path(input_dir, "CRC_NR_species_pathseq_Normalized.csv")

cat("Reading data files...\n")
cat("R group:", r_file, "\n")
cat("NR group:", nr_file, "\n\n")

r_df <- read.csv(r_file, check.names = FALSE)
nr_df <- read.csv(nr_file, check.names = FALSE)

cat("R group samples:", paste(colnames(r_df)[-1], collapse = ", "), "\n")
cat("NR group samples:", paste(colnames(nr_df)[-1], collapse = ", "), "\n")
cat("\nOutput directory:", output_dir, "\n")
cat("\nStarting to plot...\n")
cat("------------------------------------------------------------\n")

# ===== Plot for each microbe =====
for (microbe in target_microbes) {

  # Look up this microbe's data
  r_microbe_data <- r_df[r_df$name == microbe, ]
  nr_microbe_data <- nr_df[nr_df$name == microbe, ]

  # Check whether data exist
  if (nrow(r_microbe_data) == 0 && nrow(nr_microbe_data) == 0) {
    cat("Warning:", microbe, "not found in data\n")
    next
  }

  # Extract numeric values (no log transform)
  if (nrow(r_microbe_data) > 0) {
    r_values <- as.numeric(r_microbe_data[1, -1])
  } else {
    r_values <- rep(0, length(sample_groups$R))
  }

  if (nrow(nr_microbe_data) > 0) {
    nr_values <- as.numeric(nr_microbe_data[1, -1])
  } else {
    nr_values <- rep(0, length(sample_groups$NR))
  }

  # Skip if all values are 0
  if (sum(c(r_values, nr_values)) == 0) {
    cat("Skipping", microbe, ": All values are 0\n")
    next
  }

  # Prepare plotting data (NR on the left, R on the right)
  plot_data <- data.frame(
    Value = c(nr_values, r_values),
    Group = factor(c(rep("NR", length(nr_values)), rep("R", length(r_values))),
                   levels = c("NR", "R"))
  )

  # Wilcoxon rank-sum test
  if (length(nr_values) > 0 && length(r_values) > 0) {
    wilcox_result <- wilcox.test(nr_values, r_values, alternative = "two.sided")
    p_value <- wilcox_result$p.value
  } else {
    p_value <- 1.0
  }

  # Format the P-value label (3 decimals)
  if (p_value < 0.001) {
    p_label <- "p < 0.001"
  } else {
    p_label <- sprintf("P = %.3f", p_value)
  }

  # Compute y position for the significance marker
  y_max <- max(plot_data$Value)
  y_min <- min(plot_data$Value)
  y_range <- y_max - y_min
  if (y_range == 0) y_range <- y_max * 0.1
  y_pos <- y_max + y_range * 0.15

  # Microbe name (displayed in italics)
  microbe_cleaned <- gsub("_", " ", microbe)
  title_expression <- bquote(italic(.(microbe_cleaned)) ~ "(CRC)")

  # Sample sizes
  n_nr <- length(nr_values)
  n_r <- length(r_values)

  # Draw the boxplot
  p <- ggplot(plot_data, aes(x = Group, y = Value, fill = Group, color = Group)) +
    # whiskers
    stat_boxplot(geom = "errorbar", width = 0.3, linewidth = 0.5) +
    # box body (reduced width)
    geom_boxplot(
      width = 0.45,
      outlier.shape = NA,  # hide outliers
      linewidth = 0.5,
      alpha = 0.8
    ) +
    # jittered points
    geom_jitter(
      width = 0.08,
      alpha = 0.9,
      size = 1.5,
      show.legend = FALSE
    ) +
    # significance bracket
    geom_segment(
      aes(x = 1, xend = 2, y = y_pos, yend = y_pos),
      color = "black",
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    geom_segment(
      aes(x = 1, xend = 1, y = y_pos - y_range * 0.04, yend = y_pos),
      color = "black",
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    geom_segment(
      aes(x = 2, xend = 2, y = y_pos - y_range * 0.04, yend = y_pos),
      color = "black",
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    # P-value text
    annotate(
      "text",
      x = 1.5,
      y = y_pos + y_range * 0.07,
      label = p_label,
      size = 3,
      fontface = "bold"
    ) +
    # colours
    scale_fill_manual(
      values = c("NR" = colors$NR_fill, "R" = colors$R_fill),
      labels = c("NR" = paste0("NR (n=", n_nr, ")"),
                 "R" = paste0("R (n=", n_r, ")"))
    ) +
    scale_color_manual(
      values = c("NR" = colors$NR_edge, "R" = colors$R_edge),
      labels = c("NR" = paste0("NR (n=", n_nr, ")"),
                 "R" = paste0("R (n=", n_r, ")"))
    ) +
    # x-axis labels (with sample sizes)
    scale_x_discrete(
      labels = c("NR" = paste0("NR\n(n=", n_nr, ")"),
                 "R" = paste0("R\n(n=", n_r, ")"))
    ) +
    # labels
    labs(
      title = title_expression,
      y = "Normalized Score (PathSeq)",
      x = NULL
    ) +
    # y-axis range
    ylim(y_min - y_range * 0.05,
         y_pos + y_range * 0.1) +
    # theme
    theme_minimal() +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 9, family = "Helvetica", face = "bold"),
      legend.key.size = unit(0.35, "cm"),
      legend.spacing.x = unit(0.3, "cm"),
      legend.margin = margin(t = -3, b = -3),
      plot.title = element_text(size = 9, family = "Helvetica", face = "bold",
                                hjust = 0.5, vjust = 1),
      axis.text.x = element_text(size = 8, family = "Helvetica", face = "bold",
                                 hjust = 0.5, color = "black"),
      axis.text.y = element_text(size = 8, family = "Helvetica", face = "bold",
                                 color = "black"),
      axis.title.y = element_text(size = 9, family = "Helvetica", face = "bold"),
      axis.ticks.x = element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = element_line(color = "black", linewidth = 0.5),
      axis.line = element_line(linewidth = 0.5, color = "black"),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      plot.margin = margin(t = 5, r = 10, b = 10, l = 5, unit = "pt")
    )

  # Save figures
  microbe_filename <- gsub(" ", "_", microbe)

  # Save SVG
  ggsave(
    filename = file.path(output_dir, paste0("CRC_", microbe_filename, ".svg")),
    plot = p,
    width = 1.8,
    height = 2.2,
    device = "svg"
  )

  # Save PDF
  ggsave(
    filename = file.path(output_dir, paste0("CRC_", microbe_filename, ".pdf")),
    plot = p,
    width = 1.8,
    height = 2.2,
    device = "pdf"
  )

  cat("OK: Saved:", microbe_filename, sprintf("(p-value: %.4f)\n", p_value))
}

cat("------------------------------------------------------------\n")
cat("\nAll figures complete.\n")
cat("Saved to:", output_dir, "\n")
cat("============================================================\n")
