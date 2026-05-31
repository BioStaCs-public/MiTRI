# Plot the per-cancer microbial reads ratio (PathSeq reads / total reads)
# as boxplots with mean +/- standard error and per-group sample sizes.

library(ggplot2)
library(dplyr)
library(ggrepel)  # smart label placement

# Input data and output directory
input_file <- "/path/to/project/results/reads_ratio/pathseq_reads_ratio.tsv"
save_path  <- "/path/to/project/results/reads_ratio/figures/"

if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

# Read the PathSeq reads-ratio table
data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Harmonise cancer-type labels
data$cancer_type <- gsub("LUNG", "LC", data$cancer_type)
data$cancer_type <- gsub("KIDNEY", "RCC", data$cancer_type)

# Convert the reads ratio to a percentage
data$pathseq_reads_ratio_pct <- data$pathseq_reads_ratio * 100

str(data)

# Cancer-type order (roughly head-to-toe anatomical position)
cancer_order <- c(
  "OSCC",   # oral squamous cell carcinoma
  "THCA",   # thyroid cancer
  "LC",     # lung cancer
  "BRCA",   # breast cancer
  "ESCA",   # esophageal cancer (if present)
  "STAD",   # stomach adenocarcinoma
  "LIHC",   # liver hepatocellular carcinoma
  "PAAD",   # pancreatic adenocarcinoma
  "CRC",    # colorectal cancer
  "RCC",    # renal cell carcinoma
  "BLCA",   # bladder cancer
  "CESC"    # cervical squamous cell carcinoma
)

# Order by cancer type and sample status
data <- data %>%
  mutate(cancer_type = factor(cancer_type, levels = cancer_order)) %>%
  mutate(status = factor(status, levels = c("Normal", "Tumor")))

# Colours: fill = light, line = dark
status_fill_colors <- c("Normal" = "#4E92B9", "Tumor" = "#d64863")
status_line_colors <- c("Normal" = "#20364F", "Tumor" = "#7e3d4e")

# Per-group sample sizes
sample_counts <- data %>%
  group_by(cancer_type, status) %>%
  summarise(n = n(), .groups = 'drop')

# Summary statistics (including sample size)
summary_data <- data %>%
  group_by(cancer_type, status) %>%
  summarise(
    mean = mean(pathseq_reads_ratio_pct, na.rm = TRUE),
    sd = sd(pathseq_reads_ratio_pct, na.rm = TRUE),
    se = sd / sqrt(n()),
    n = n(),
    .groups = 'drop'
  )

# Boxplot + error bars + sample-size annotation
p_boxplot_error <- ggplot(data, aes(x = cancer_type, y = pathseq_reads_ratio_pct, fill = status)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.6,
    outlier.shape = NA,  # hide outliers; raw points are drawn separately
    alpha = 0.7,
    color = "black",
    linewidth = 0.5
  ) +
  # Overlay raw data points
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.75
    ),
    size = 0.8,
    alpha = 0.3,
    shape = 16
  ) +
  # Error bars (mean +/- standard error)
  geom_errorbar(
    data = summary_data,
    aes(x = cancer_type, y = mean, ymin = mean - se, ymax = mean + se),
    position = position_dodge(width = 0.75),
    width = 0.25,
    linewidth = 0.8,
    color = "black"
  ) +
  # Mean point at the top of each error bar
  geom_point(
    data = summary_data,
    aes(x = cancer_type, y = mean),
    position = position_dodge(width = 0.75),
    size = 2,
    shape = 23,  # diamond
    fill = "white",
    color = "black",
    stroke = 1
  ) +
  # Sample-size labels
  geom_text(
    data = summary_data,
    aes(x = cancer_type, y = mean + se, label = paste0("n=", n)),
    position = position_dodge(width = 0.75),
    vjust = -0.8,
    size = 2.5,
    fontface = "bold",
    family = "Helvetica",
    color = "black"
  ) +
  scale_fill_manual(values = status_fill_colors) +
  theme_minimal() +
  labs(
    x = "Cancer Type",
    y = "PathSeq Reads Ratio (%)",
    fill = "Status"
  ) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 9, family = "Helvetica", face = "bold"),
    legend.text = element_text(size = 9, family = "Helvetica", face = "bold"),
    axis.text.x = element_text(size = 8, family = "Helvetica", face = "bold",
                               angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 8, family = "Helvetica", face = "bold", color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 9, family = "Helvetica", face = "bold"),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
  ) +
  # Expand the y-axis to fit the sample-size labels
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# Save the figure
ggsave(
  file.path(save_path, "pathseq_reads_ratio_boxplot_errorbar.svg"),
  plot = p_boxplot_error,
  width = 6,
  height = 4,
  device = "svg"
)

ggsave(
  file.path(save_path, "pathseq_reads_ratio_boxplot_errorbar.pdf"),
  plot = p_boxplot_error,
  width = 6,
  height = 4,
  device = "pdf"
)

print("Boxplot with error bars saved.")
print("File name: pathseq_reads_ratio_boxplot_errorbar")
