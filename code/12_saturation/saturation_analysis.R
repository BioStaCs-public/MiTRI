library(iNEXT)
library(stringr)
library(dplyr)
library(readr)
library(ggplot2)
#
# library(dplyr)
#
# # Cancer types
# cancer_types <- c("OSCC", "BRCA", "PAAD", "CESC", "CRC", "LUNG", "KIDNEY", "BLCA", "STAD", "LIHC", "THCA")
#
# # Path template for microbial data files
# microbe_base_path <- "/path/to/project/results_V3/cancers_V3.1/"
#
# # Read and process data
# data_list <- list()  # store microbial data for all cancer types
#
# for(cancer_type in cancer_types) {
#   message("Processing: ", cancer_type)
#
#   # 1. Read Normal
#   normal_file <- paste0(
#     microbe_base_path, cancer_type,
#     "/01.group_matrix/", cancer_type, "_Normal_species_pathseq.csv"
#   )
#   normal_data <- read.csv(normal_file, check.names = FALSE)
#
#   # 2. Read Tumor
#   tumor_file <- paste0(
#     microbe_base_path, cancer_type,
#     "/01.group_matrix/", cancer_type, "_Tumor_species_pathseq.csv"
#   )
#   tumor_data <- read.csv(tumor_file, check.names = FALSE)
#
#   # 3. Merge Normal + Tumor by microbe name
#   # assumes both have a 'name' column
#   merged_data <- full_join(normal_data, tumor_data, by = "name")
#
#   # 4. Extract microbe names
#   microbe_names <- merged_data$name
#
#   # 5. Convert any value > 0 to 1, otherwise 0
#   merged_clean <- merged_data %>%
#     select(-name) %>%
#     as.matrix() %>%
#     apply(2, function(x) ifelse(x > 0, 1, 0))
#
#   # 6. Set row names to species names and clean formatting
#   rownames(merged_clean) <- microbe_names
#   rownames(merged_clean) <- gsub("[\\[\\]]", "", rownames(merged_clean))
#   rownames(merged_clean) <- gsub(" ", "_", rownames(merged_clean))
#
#   # 7. Store in the list
#   data_list[[cancer_type]] <- merged_clean
# }
# 
# save(
#   data_list,
#   file = "/path/to/project/results_V3/summary/saturation/saturation_11_normal_tumor_species_pathseq.RData"
# )

# data_list <- lapply(data_list, function(m) {
#   m[is.na(m)] <- 0   # replace NA with 0
#   mode(m) <- "numeric"
#   return(m)
# })
# # 
# out <- iNEXT(data_list, q = 0, datatype = "incidence_raw")
# save(out, file = "/path/to/project/results/summary_V2/saturation/out.RData")

# Load the saved out object
load("/path/to/project/results_V3/summary/saturation/out.RData")
# Cancer-type to colour mapping
cancer_colors <- c(
  "OSCC" = "#9d001f",      # deep red
  "BRCA" = "#316480",       # deep sea blue
  "PAAD" = "#d87c30",       # ochre orange
  "CESC" = "#5f4359",       # deep purple
  "CRC" = "#4a7c4a",        # deep green
  "LUNG" = "#20364f",       # deep navy
  "KIDNEY" = "#96dfd2",     # ice blue
  "BLCA" = "#6fbad0",       # bright lake blue
  "STAD" = "#b57a4a",       # deep terracotta
  "LIHC" = "#3f6e3f",       # deep leaf green
  "THCA" = "#f2c27b"        # honey yellow
)

p <- ggiNEXT(x = out, type = 1, se = TRUE,
                grey=FALSE) +
  theme_classic() +
  ylab("Number of species") +
  xlab("Sample size") +
  scale_color_manual(values = cancer_colors) +   # apply colour mapping to the curves
  scale_shape_manual(values = rep(17, 11)) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 9, family = "Helvetica", face = "bold"), # Legend labels font
    legend.key.size = unit(0.35, "cm"),  # legend key size
    legend.spacing.x = unit(0.3, "cm"),  # spacing between legend items
    legend.margin = margin(t = -4, b = -4),  # tighten legend margins
    plot.title = element_text(size = 9, family="Helvetica", face = "bold", hjust = 0.5, vjust = 1),
    axis.text.x = element_text(size = 8, family="Helvetica", face = "bold", hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = 8, family="Helvetica", face = "bold", color = "black"),
    axis.title.x = element_text(size = 9, family="Helvetica", face = "bold"),
    axis.title.y = element_text(size = 9, family="Helvetica", face = "bold"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5), # x-axis ticks
    axis.ticks.y = element_line(color = "black", linewidth = 0.5), # y-axis ticks
    axis.line = element_line(linewidth = 0.5, color = "black"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    plot.margin = margin(t = 5, r = 10, b = 10, l = 5, unit = "pt")
  )

p
save_path = "/path/to/project/results_V3/summary/saturation/saturation"
# Save figures
ggsave(file.path(save_path, paste0("Saturation_pretein", ".svg")), plot = p, width = 4, height = 3, device = "svg")
ggsave(file.path(save_path, paste0("Saturation_pretein", ".pdf")), plot = p, width = 4, height = 3, device = "pdf")

# # 
# library(ggplot2)
# library(dplyr)
# library(patchwork)
# 
# # Load the saved out object
# load("/path/to/project/results_V3/summary/saturation/out.RData")
#
# # Extract rarefaction-curve data
# size_based_data <- out$iNextEst$size_based
# # Convert Assemblage to a factor and sort alphabetically
# size_based_data$Assemblage <- factor(size_based_data$Assemblage, levels = sort(unique(size_based_data$Assemblage)))
#
# # Slope scaling factor (to display slope and qD on the same plot)
# scale_factor <- 20
#
# # Container for the list of plots
# plot_list <- list()
#
# # Generate a subplot for each cancer type
# for (cancer in unique(size_based_data$Assemblage)) {
#   cancer_data <- size_based_data %>%
#     filter(Assemblage == cancer, Method == "Rarefaction") %>%
#     arrange(t)
#
#   n <- nrow(cancer_data)
#   slope <- rep(NA, n)
#   for (i in 2:(n - 1)) {
#     slope[i] <- (cancer_data$qD[i + 1] - cancer_data$qD[i - 1]) /
#       (cancer_data$t[i + 1] - cancer_data$t[i - 1])
#   }
#   cancer_data$slope <- slope
#
#   # Percent change in qD
#   cancer_data$qD_pct_change <- c(NA, diff(cancer_data$qD) / head(cancer_data$qD, -1) * 100)
#
#   # Percent-change threshold (e.g. 5%)
#   pct_threshold <- 5
#
#   # Points below the threshold are treated as saturation points
#   saturation_point <- which(abs(cancer_data$qD_pct_change) < pct_threshold)
#
#   # Set the grey-line yintercept to the qD at the saturation point
#   if (length(saturation_point) > 0) {
#     saturation_value <- cancer_data$qD[saturation_point[1]]  # qD at the first saturation point
#     saturation_t <- cancer_data$t[saturation_point[1]]  # corresponding sample size t
#   } else {
#     saturation_value <- max(cancer_data$qD)  # if none found, use the maximum
#     saturation_t <- max(cancer_data$t)  # corresponding maximum sample size
#   }
#
#   # Build the plot
#   p <- ggplot(cancer_data, aes(x = t)) +
#     geom_line(aes(y = qD), color = "#2a749a", size = 1) +
#     geom_point(aes(y = qD), color = "#2a749a", size = 1.5) +
#     geom_line(aes(y = slope * scale_factor), color = "#9d001f", size = 0.8, linetype = "dashed") +
#     geom_hline(yintercept = saturation_value, linetype = "dotted", color = "gray40") +  # qD at saturation point
#     geom_vline(xintercept = saturation_t, linetype = "dotted", color = "gray40") +  # add vertical dashed line
#     annotate("text", x = saturation_t, y = saturation_value, label = paste(round(saturation_t, 2)),
#              color = "black", angle = 0, hjust = -0.1, vjust = 1.5, size = 3,family = "Helvetica", face = "bold") +  # add label
#     scale_y_continuous(
#       name = "Species Richness (qD)",
#       sec.axis = sec_axis(~./scale_factor, name = "Slope (ΔqD / Δt)")
#     ) +
#     labs(
#       title = cancer,
#       x = "Sample Size (t)"
#     ) +
#     theme_minimal(base_size = 9) +
#     theme(
#       plot.title = element_text(size = 8, family = "Helvetica", face = "bold", hjust = 0.5),
#       axis.title.y.left = element_text(color = "#2a749a", size = 8, family = "Helvetica", face = "bold"),
#       axis.text.y.left = element_text(color = "#2a749a", size = 8, family = "Helvetica", face = "bold"),
#       axis.title.y.right = element_text(color = "#9d001f", size = 8, family = "Helvetica", face = "bold"),
#       axis.text.y.right = element_text(color = "#9d001f", size = 8, family = "Helvetica", face = "bold"),
#       axis.text.x = element_text(size = 8, family = "Helvetica", face = "bold")
#     )
# 
#   plot_list[[cancer]] <- p
# }
# # Sort cancer types alphabetically
# plot_list <- plot_list[order(names(plot_list))]
# # Combine into a 3-column layout
# final_plot <- wrap_plots(plot_list, ncol = 3) +
#   plot_annotation(
#     title = "Rarefaction Curves and Slope Trends Across Cancer Types",
#     theme = theme(plot.title = element_text(size = 9, family = "Helvetica", face = "bold", hjust = 0.5))
#   )
# 
# # Display
# # print(final_plot)
# save_path <- "/path/to/project/results_V3/summary/saturation/saturation"
#
# # Save as PDF or PNG
# ggsave(file.path(save_path, "Combined_Saturation_Slope_Plots.pdf"), final_plot, width = 11, height = 9, device = "pdf")
# ggsave(file.path(save_path, "Combined_Saturation_Slope_Plots.svg"), final_plot, width = 11, height = 9, device = "svg")
