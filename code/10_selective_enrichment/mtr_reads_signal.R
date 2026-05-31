# Load required libraries
library(readxl)   # read Excel files
library(dplyr)    # data manipulation
library(ggplot2)  # plotting
library(tidyr)    # data reshaping
library(ggforce)

# Output directory
save_path <- "/path/to/project/results_V3/summary/MTR_enrichment_analysis/figures/"
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

# Load the combined enrichment-score table
data_file <- "/path/to/project/results_V3/summary/MTR_enrichment_analysis/mtr_enrichment_scores_all.csv"

if (file.exists(data_file)) {
  all_data <- read.csv(data_file)
}
# Recode abbreviations immediately after loading (before filtering)
all_data <- all_data %>%
  mutate(Cancer_Type = recode(Cancer_Type,
                              "LUNG"   = "LC",
                              "KIDNEY" = "RCC",
                              "STAD"   = "GC"))

# Cancer-type display order
cancer_order <- c(
  "BLCA", "BRCA", "CESC","CRC", "GC", "LC",
  "LIHC", "OSCC", "PAAD"
)

# Fix factor order
all_data$Cancer_Type <- factor(all_data$Cancer_Type, levels = cancer_order)

all_data <- all_data[all_data$P_Value < 0.05, ]


# Keep rows with P_Value < 0.05
p_value_filtered <- all_data[all_data$P_Value < 0.05, ]

# Of those, keep Log2_Enrichment > 2
log2_enrichment_filtered <- p_value_filtered[p_value_filtered$Log2_Enrichment > 2, ]

# Count qualifying rows vs. all P_Value < 0.05 rows
num_log2_enrichment_positive <- nrow(log2_enrichment_filtered)  # count with high enrichment
num_p_value_filtered <- nrow(p_value_filtered)  # count with P_Value < 0.05

# Proportion of qualifying rows
proportion <- num_log2_enrichment_positive / num_p_value_filtered

# As a percentage
proportion_percent <- sprintf("%.1f%%", proportion * 100)

# Print results
cat("Total with P_Value < 0.05: ", num_p_value_filtered, "\n")
cat("Count with high Log2_Enrichment: ", num_log2_enrichment_positive, "\n")
cat("Proportion qualifying: ", proportion_percent, "\n")


# Order by cancer-type levels
all_data <- all_data %>%
  mutate(Cancer_Type = factor(Cancer_Type, levels = cancer_order))

# Replace underscores with spaces
all_data$Reference <- gsub("_", " ", all_data$Microbe)

# # Optionally drop specific references
# all_data <- all_data %>%
#   filter(Reference != "Sphingobacterium shayense") %>%
#   filter(Reference != "Pasteurella multocida") %>%
#   # drop the CRC Fusobacterium periodonticum record
#   filter(!(Cancer_Type == "CRC" & Reference == "Fusobacterium periodonticum"))

# Unique microbes
microbes <- unique(all_data$Reference)
microbes <- sort(microbes)

microbe_colors <- setNames(
  c(
    # Achromobacter
    "#4e92b9",  # Achromobacter xylosoxidans

    # Acinetobacter
    "#7da77c",  # Acinetobacter johnsonii
    "#5b8c5a",  # Acinetobacter nosocomialis

    # Campylobacter
    "#d87c30",  # Campylobacter ureolyticus

    # Cutibacterium
    "#f3eed9",  # Cutibacterium acnes

    # Dialister
    "#9fb7c1",  # Dialister invisus

    # Enterobacter
    "#c9b27c",  # Enterobacter asburiae

    # Enterococcus
    "#96dfd2",  # Enterococcus casseliflavus
    "#6fbad0",  # Enterococcus innesii

    # Escherichia
    "#dcd9cf",  # Escherichia coli

    # Fusobacterium (genus shares a red palette)
    "#9d001f",  # Fusobacterium animalis
    "#7e3d4e",  # Fusobacterium canifelinum
    "#e19992",  # Fusobacterium hwasookii
    "#b3002d",  # Fusobacterium massiliense
    "#d14d65",  # Fusobacterium nucleatum
    "#f2a2a0",  # Fusobacterium periodonticum
    "#ffcccc",  # Fusobacterium polymorphum
    "#d87c30",  # Fusobacterium pseudoperiodonticum
    "#e4a26a",  # Fusobacterium varium
    "#c18c5d",  # Fusobacterium vincentii

    # Helicobacter
    "#4a7c4a",  # Helicobacter hepaticus
    "#6a9c6a",  # Helicobacter pullorum
    "#88bfae",  # Helicobacter pylori

    # Neisseria
    "#8c8c8c",  # Neisseria gonorrhoeae
    "#a5a5a5",  # Neisseria mucosa

    # Parvimonas
    "#8b5a2b",  # Parvimonas micra

    # Pasteurella
    "#a8cfa3",  # Pasteurella multocida

    # Peptostreptococcus
    "#7b6d5f",  # Peptostreptococcus stomatis

    # Porphyromonas
    "#3f6e3f",  # Porphyromonas gingivalis

    # Prevotella
    "#8c678a",  # Prevotella amnii
    "#a27b9e",  # Prevotella intermedia
    "#b497b7",  # Prevotella melaninogenica

    # Proteus
    "#5f4359",  # Proteus mirabilis

    # Pseudomonas
    "#2c4d6b",  # Pseudomonas putida

    # Shigella
    "#316480",  # Shigella sonnei

    # Sphingomonas
    "#20364f",  # Sphingomonas paucimobilis

    # Stenotrophomonas
    "#5a9bb1",  # Stenotrophomonas maltophilia

    # Streptococcus
    "#e6e2d3",  # Streptococcus anginosus

    # Veillonella
    "#6d6d6d"   # Veillonella parvula
  ),
  c(
    "Achromobacter xylosoxidans",

    "Acinetobacter johnsonii",
    "Acinetobacter nosocomialis",

    "Campylobacter ureolyticus",

    "Cutibacterium acnes",

    "Dialister invisus",

    "Enterobacter asburiae",

    "Enterococcus casseliflavus",
    "Enterococcus innesii",

    "Escherichia coli",

    "Fusobacterium animalis",
    "Fusobacterium canifelinum",
    "Fusobacterium hwasookii",
    "Fusobacterium massiliense",
    "Fusobacterium nucleatum",
    "Fusobacterium periodonticum",
    "Fusobacterium polymorphum",
    "Fusobacterium pseudoperiodonticum",
    "Fusobacterium varium",
    "Fusobacterium vincentii",

    "Helicobacter hepaticus",
    "Helicobacter pullorum",
    "Helicobacter pylori",

    "Neisseria gonorrhoeae",
    "Neisseria mucosa",

    "Parvimonas micra",

    "Pasteurella multocida",

    "Peptostreptococcus stomatis",

    "Porphyromonas gingivalis",

    "Prevotella amnii",
    "Prevotella intermedia",
    "Prevotella melaninogenica",

    "Proteus mirabilis",

    "Pseudomonas putida",

    "Shigella sonnei",

    "Sphingomonas paucimobilis",

    "Stenotrophomonas maltophilia",

    "Streptococcus anginosus",

    "Veillonella parvula"
  )
)


# Violin plot with a sina scatter overlay
p <- ggplot(all_data, aes(x = Cancer_Type, y = Log2_Enrichment)) +
  geom_violin(scale = "width") +  # violin
  geom_sina(size = 0.025, aes(color = Reference, group = Cancer_Type), alpha = 0.9, scale = FALSE) +
  scale_color_manual(values = microbe_colors) +
  theme_minimal() +  # minimal theme
  labs(x = "Cancer type", y = "log2(Enrichment)")  +
  scale_y_continuous(breaks = c(-5, 0, 5, 10,15)) +  # y-axis breaks
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 6, family = "Helvetica", face = "bold"), # Legend labels font
    legend.key.size = unit(0.25, "cm"),  # legend key size
    legend.spacing.x = unit(0.3, "cm"),  # spacing between legend items
    legend.margin = margin(t = -4, b = -4),  # tighten legend margins
    plot.title = element_text(size = 9, family="Helvetica", face = "bold", hjust = 0.5, vjust = 1),
    axis.text.x = element_text(size = 8, family="Helvetica", face = "bold", angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(size = 8, family="Helvetica", face = "bold", color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 9, family="Helvetica", face = "bold"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    # axis.line.x = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    plot.margin = margin(t = 5, r = 10, b = 10, l = 5, unit = "pt")
  )

p <- p +
  guides(
    color = guide_legend(
      override.aes = list(size = 2)
    )
  ) +
  theme(
    legend.key.size = unit(0.3, "cm")  # slightly larger row height
  )


# Save files
ggsave(
  file.path(save_path, paste0("Active_Microbe_MTR_reads_legend.svg")),
  plot = p,
  width = 7,
  height = 2.0,
  device = "svg"
)

ggsave(
  file.path(save_path, paste0("Active_Microbe_MTR_reads_legend.pdf")),
  plot = p,
  width = 7,
  height = 2.0,
  device = "pdf"
)
