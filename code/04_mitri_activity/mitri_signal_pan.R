library(ggplot2)
library(ggpubr)
library(grid)
library(dplyr)
library(data.table)
library(parallel)
library(openxlsx)
library(ggpubr)
library(extrafont)
library(readxl)

# # Import system fonts
# font_import()

# # Load fonts
# loadfonts()

# ===== Sample-count dictionary =====
sample_counts <- list(
  BLCA = list(PT = 44, T = 167),
  BRCA = list(PT = 25, T = 146),
  CESC = list(PT = 15, T = 45),
  CRC = list(PT = 89, T = 89),
  KIDNEY = list(PT = 10, T = 38),
  LIHC = list(PT = 92, T = 267),
  LUNG = list(PT = 137, T = 139),
  OSCC = list(PT = 35, T = 35),
  PAAD = list(PT = 22, T = 22),
  STAD = list(PT = 187, T = 211),
  THCA = list(PT = 56, T = 26)
)

# ===== Cancer-type display-name mapping =====
cancer_display_names <- list(
  KIDNEY = "RCC",
  LUNG = "LC",
  STAD = "GC"
)

# Helper: resolve a display name
get_display_name <- function(cancer_type) {
  if (cancer_type %in% names(cancer_display_names)) {
    return(cancer_display_names[[cancer_type]])
  } else {
    return(cancer_type)
  }
}
# ==========================================

# Load data
# cancer_type <- c("CRC")
# cancer_type <- c("BRCA", "KIDNEY", "OSCC", "BLCA", "CESC", "LIHC", "PAAD", "LUNG", "STAD", "THCA")
# , "KIDNEY", "PAAD", "THCA", "LUNG", "BRCA", "CESC"
cancer_type <- c("PAAD")


# ===== Look up the sample counts for this cancer type =====
if (!cancer_type %in% names(sample_counts)) {
  stop(paste0("Error: no sample-count information found for cancer type ", cancer_type, "!"))
}

n_PT <- sample_counts[[cancer_type]]$PT
n_T <- sample_counts[[cancer_type]]$T

# ===== Resolve the display name =====
cancer_display <- get_display_name(cancer_type)

# Build labels
label_T <- paste0("T (n=", n_T, ")")
label_PT <- paste0("PT (n=", n_PT, ")")

cat("Cancer type:", cancer_type, "\n")
cat("Display name:", cancer_display, "\n")
cat("T sample count:", n_T, "\n")
cat("PT sample count:", n_PT, "\n")
# ========================================

base_path <- "/path/to/project/results_V3/cancers_V3.1"

save_path <- file.path(base_path, cancer_type, "04.activity/figures")
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

file_name <- paste0(cancer_type, "_Tumor_Normal_MiTRI.xlsx")
file_path <- file.path(base_path, cancer_type, "04.activity", file_name)

tumor_data <- read_excel(file_path, sheet = "Tumor")
normal_data <- read_excel(file_path, sheet = "Normal")

tumor_data <- tumor_data %>% select(-'Sample_Name')
normal_data <- normal_data %>% select(-'Sample_Name')

# ===== Use the dynamic labels =====
tumor_data$Type <- label_T
normal_data$Type <- label_PT
# ================================

combined_data <- rbind(tumor_data, normal_data)

# Reshape from wide to long format
long_data <- reshape2::melt(combined_data,
                            id.vars = c("Reference", "Scenario", "Type"),
                            variable.name = "Group",
                            value.name = "Value")
# log-transform the Value column
long_data$Value <- log(long_data$Value + 1)  # note: +1 avoids log(0)


# For each microbe, draw a boxplot and add significance testing
unique_microbes <- unique(long_data$Reference)
# Build the file path
file_path <- paste0("/path/to/project/data/CSVs_20251114/", cancer_type, "_paired.csv")

# Read the taxon column from the CSV
crc_data <- read.csv(file_path, stringsAsFactors = FALSE)
valid_microbes <- intersect(unique_microbes, crc_data$Taxon)


colors <- c("#9d001f","#7e3d4e","#e19992", "#f3eed9",
            "#DCD9CF", "#96DFD2", "#4E92B9", "#316480", "#20364F")

colors <- c("#9d001f","#d64863","#ff82ab", "#540051","#d06bc6",
            "#ff70ab", "#0083cd", "#acecff", "#00beff","#002f50", "#2a749a")


for (microbe in valid_microbes) {

  microbe_data <- subset(long_data, Reference == microbe)
  # print(microbe)
  # Per-group max mean plus an offset, used for placing p-value labels
  y_positions <- aggregate(Value ~ Group, data = microbe_data, max)$Value + 0.01
  microbe_cleaned <- gsub("_", " ", microbe) # replace underscores with spaces

  # Test Type differences within each Group
  stat_test <- compare_means(
    Value ~ Type,                # compare across Type
    group.by = "Group",          # within each Group
    data = microbe_data,
    method = "wilcox.test"       # Wilcoxon test
  )

  # Build bracket data dynamically
  # Set bracket height based on each group's data maximum
  max_y <- microbe_data %>%
    group_by(Group) %>%
    summarise(max_value = max(Value, na.rm = TRUE)) %>%  # na.rm = TRUE
    filter(is.finite(max_value))  # drop non-finite values

  # Combine maxima with significance results to compute dynamic y positions
  stat_test <- stat_test %>%
    left_join(max_y, by = c("Group")) %>%
    mutate(
      y = max_value + 0.3 + (row_number() - 1) * 0.05, # dynamic y position
      xstart = as.numeric(Group) - 0.25,              # left bar offset
      xend = as.numeric(Group) + 0.25                 # right bar offset
    )

  p <- ggplot(microbe_data, aes(x = Group, y = Value, fill = Type)) +
    # whiskers
    stat_boxplot(position = position_dodge(width = 0.8),
                 geom = "errorbar", width = 0.4, aes(color = Type),) +  # whisker width
    # box body
    geom_boxplot(position = position_dodge(width = 0.8),
                 outlier.shape = NA, # hide outliers (avoid overlap)
                 width = 0.6, # box width
                 aes(color = Type),
                 lwd = 0.5,  # box line width
    ) +
    # jittered points
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2),
                alpha = 0.9, size = 0.5, aes(color = Type)) +
    # significance labels
    stat_compare_means(aes(group = Type),
                       method = "wilcox.test",
                       label = "p.format",
                       label.y = stat_test$y-0.3, # slightly above the bracket
                       size = 4) +
    # bracket horizontal line
    geom_segment(data = stat_test,
                 aes(x = xstart, xend = xend, y = y, yend = y),
                 inherit.aes = FALSE, # do not inherit global aes
                 linewidth = 0.5) +  # line width
    # bracket vertical ticks
    geom_segment(data = stat_test,
                 aes(x = xstart, xend = xstart, y = y - 0.15, yend = y+0.01),
                 inherit.aes = FALSE, linewidth = 0.5) +
    geom_segment(data = stat_test,
                 aes(x = xend, xend = xend, y = y - 0.15, yend = y+0.01),
                 inherit.aes = FALSE, linewidth = 0.5) +
    labs(
      # ===== Use the display name rather than the raw cancer code =====
      title = paste0(microbe_cleaned, " (", cancer_display, ")"),
      # ===============================================
      y = "log (MiTRI)",
      x = "Baseline") +
    scale_x_discrete(
      labels = c(
        "BG-1",
        "BG-2",
        "BG-3",
        "BG-4"
      ),
      expand = c(0, 0) # group spacing; increase for more gap
    )+
    # ===== Use dynamic labels and named vectors =====
    theme_minimal() +
    scale_fill_manual(values = setNames(c("#d64863", "#4E92B9"), c(label_T, label_PT))) +
    scale_color_manual(values = setNames(c("#7e3d4e", "#20364F"), c(label_T, label_PT))) +
    # ===========================================
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 9, family = "Helvetica", face = "bold"), # Legend labels font
      legend.key.size = unit(0.35, "cm"),  # legend key size
      legend.spacing.x = unit(0.3, "cm"),  # spacing between legend items
      legend.margin = margin(t = -3, b = -3),  # tighten legend margins
      plot.title = element_text(size = 9, family="Helvetica", face = "bold", hjust = 0.5, vjust = 1),
      axis.text.x = element_text(size = 8, family="Helvetica", face = "bold", hjust = 0.5, color = "black"),
      axis.text.y = element_text(size = 8, family="Helvetica", face = "bold", color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 9, family="Helvetica", face = "bold"),
      axis.ticks.x = element_line(color = "black", linewidth = 0.5), # x-axis ticks
      axis.ticks.y = element_line(color = "black", linewidth = 0.5), # y-axis ticks
      axis.line = element_line(linewidth = 0.5, color = "black"),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      plot.margin = margin(t = 5, r = 10, b = 10, l = 5, unit = "pt")
    )

  # p
  # print(p)

  # ===== Use the display name in the file name too =====
  # Save files
  ggsave(
    file.path(save_path, paste0(cancer_display,"_", microbe, "_format.svg")),
    plot = p,
    width = 2.4,
    height = 2.4,
    device = "svg"
  )

  ggsave(
    file.path(save_path, paste0(cancer_display,"_", microbe, "_format.pdf")),
    plot = p,
    width = 2.4,
    height = 2.4,
    device = "pdf"
  )
  # ========================================
}
