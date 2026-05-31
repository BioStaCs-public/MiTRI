library(dplyr)
library(ggplot2)
library(ggbreak)

# ==============================================================================
# 1. Core statistical function (DFR 2x)
# ==============================================================================
# High-precision DFR statistic
elsdfr2x_raw <- function(data, nExp, nCtl, nameCtl, ho.c=log10(2), alpha=0.05, B=1000, seed=9456845){
  data <- as.data.frame(data)
  data <- data[order(data[,1],data[,2]),]
  ncdat <- data[data[,3]==nameCtl,]
  dat <- data.frame(matrix(NA,nrow=NROW(data[data[,3]!=nameCtl,]),ncol=3+nExp+nCtl))
  dat[,1:(3+nExp)] <- data[data[,3]!=nameCtl,1:(3+nExp)]
  names(dat) <- c(names(data)[1:(3+nExp)],paste("c",1:nCtl,sep=""))
  for(id in unique(data[,1])){
    for(day in unique(data[,2])){
      dat[dat[,1]==id & dat[,2]==day,3+nExp+(1:nCtl)] <- ncdat[ncdat[,1]==id & ncdat[,2]==day,-(1:3)]
    }}
  idx <- which(rowSums(dat[,-(1:3)]==0, na.rm=TRUE) > 0)
  if (length(idx)>0){ dat[idx,-(1:3)] <- dat[idx,-(1:3)] + 1 }
  dat[,-(1:3)] <- log10(dat[,-(1:3)])

  ind <- unique(dat[,1:2])
  adjp <- teststat <- pos.call <- NULL

  for(i in 1:NROW(ind)){
    temp.dat <- as.matrix(dat[dat[,1]==ind[i,1] & dat[,2]==ind[i,2],-(1:3)])
    temp.id <- ind[i,1]
    temp.day <- ind[i,2]
    temp.pep <- dat[dat[,1]==ind[i,1] & dat[,2]==ind[i,2],3]
    if (length(temp.pep)==1) temp.dat <- t(temp.dat)
    temp.adjp <- temp.teststat <- rep(NA, NROW(temp.dat))

    nas <- is.na(temp.dat)
    nCtlc <- nCtl - sum(nas[1,nExp+(1:nCtl)])
    nExpC <- nExp - apply(nas[,1:nExp],1,sum,na.rm=TRUE)
    idx <- (nExpC>=3 & nCtlc>=3) | (nExpC>=2 & nCtlc>=4)
    if (sum(idx)>0){
      temp.dat <- as.matrix(temp.dat[idx,])
      if (sum(idx)==1 & length(idx)>1) temp.dat <- t(temp.dat)
      out <- bootfn(dat=temp.dat, nExp=nExp ,nCtl=nCtl, ho.c=ho.c, B=B, seed=seed)
      temp.teststat[idx] <- out$tstat
      temp.adjp[idx] <- out$tadjp
    }
    teststat <- c(teststat, temp.teststat)
    adjp <- c(adjp, temp.adjp)
  }
  pos.call <- ifelse(adjp <= alpha,1,0)
  output <- data.frame(cbind(dat[,1:3], teststat, adjp, pos.call))
  colnames(output) <- c("ptid","day","peptide","t-stat","adjp","pos")
  output
}

# ==============================================================================
# 2. Data entry and statistical computation
# ==============================================================================

# Raw data (used to compute P values)
# Note: only Neg and Pool are included, since the DFR method compares the two
my_data_stats <- data.frame(
  id = rep("Patient_01", 4),
  day = rep("Batch_1", 4),
  antigen = c("Peptide_A", "Peptide_B", "Neg_Ctrl", "Neg_Ctrl"),
  rep1 = c(795, 961, 270, 538),
  rep2 = c(867, 944, 271, 528),
  rep3 = c(725, 870, 242, 468),
  stringsAsFactors = FALSE
)

# Compute P values (1.5x threshold, 100,000 resamples)
result_stats <- elsdfr2x_raw(
  data = my_data_stats,
  nExp = 3, nCtl = 3, nameCtl = "Neg_Ctrl",
  ho.c = log10(1.5), B = 100000
)
# Correct P = 0
result_stats$adjp_corrected <- ifelse(result_stats$adjp == 0, 1/100000, result_stats$adjp)

# Extract P values
pval_A <- result_stats$adjp_corrected[result_stats$peptide == "Peptide_A"]
pval_B <- result_stats$adjp_corrected[result_stats$peptide == "Peptide_B"]

# ==============================================================================
# 3. Build plotting data (strict ordering)
# ==============================================================================

# Manually build the plotting data frame.
# Plot_Group controls the X-axis order (1-6):
# 1: Neg 1 (Peptide A Neg)
# 2: Pool 1 (Peptide A Stim)
# 3: Pos (Peptide A Pos)
# 4: Neg 2 (Peptide B Neg)
# 5: Pool 2 (Peptide B Stim)
# 6: Pos (Peptide B Pos)

plot_data <- rbind(
  # --- Group A ---
  data.frame(Plot_Group = "A_Neg",  Label = "Neg 1",  Spots = c(270, 271, 242)),
  data.frame(Plot_Group = "A_Pool", Label = "Pool 1", Spots = c(795, 867, 725)),
  data.frame(Plot_Group = "A_Pos",  Label = "Pos",    Spots = c(2880, 3311)),

  # --- Group B ---
  data.frame(Plot_Group = "B_Neg",  Label = "Neg 2",  Spots = c(538, 528, 468)),
  data.frame(Plot_Group = "B_Pool", Label = "Pool 2", Spots = c(961, 944, 870)),
  data.frame(Plot_Group = "B_Pos",  Label = "Pos",    Spots = c(2764, 2619))
)

# Key: set factor order so bars are drawn left-to-right as intended
plot_data$Plot_Group <- factor(plot_data$Plot_Group,
                               levels = c("A_Neg", "A_Pool", "A_Pos", "B_Neg", "B_Pool", "B_Pos"))

# Compute mean and SD
summary_df <- plot_data %>%
  group_by(Plot_Group, Label) %>%
  summarise(
    mean_spots = mean(Spots),
    sd_spots = sd(Spots),
    .groups = 'drop'
  )

# ==============================================================================
# 4. Prepare plotting parameters (shapes, positions, annotations)
# ==============================================================================

# Significance-star helper
get_stars <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

label_A <- get_stars(pval_A)
label_B <- get_stars(pval_B)

# Bracket line position (based on max of Neg and Pool, excluding Pos)
# Group A (X=1, X=2)
max_A <- max(plot_data$Spots[plot_data$Plot_Group %in% c("A_Neg", "A_Pool")])
y_line_A <- max_A + 50
y_text_A <- y_line_A + 50

# Group B (X=4, X=5)
max_B <- max(plot_data$Spots[plot_data$Plot_Group %in% c("B_Neg", "B_Pool")])
y_line_B <- max_B + 50
y_text_B <- y_line_B + 50

# ==============================================================================
# 5. Draw the figure
# ==============================================================================

p_final <- ggplot(summary_df, aes(x = Plot_Group, y = mean_spots)) +

  # 1. Bars (white fill, black border)
  geom_bar(stat = "identity", fill = "white", color = "black", width = 0.6) +

  # 2. Error bars
  geom_errorbar(aes(ymin = ifelse(mean_spots - sd_spots < 0, 0, mean_spots - sd_spots),
                    ymax = mean_spots + sd_spots),
                width = 0.25, linewidth = 0.8) +

  # 3. Jittered points (custom shapes)
  # aes(shape = Plot_Group) assigns shapes automatically
  geom_jitter(data = plot_data, aes(x = Plot_Group, y = Spots, shape = Plot_Group),
              size = 1, width = 0.1, stroke = 1) +

  # Assign shapes manually (circle, square, diamond, triangle, etc.)
  scale_shape_manual(values = c(
    "A_Neg"  = 16, # filled circle
    "A_Pool" = 15, # filled square
    "A_Pos"  = 18, # filled diamond
    "B_Neg"  = 1,  # open circle (distinguishes Group 2)
    "B_Pool" = 17, # filled triangle
    "B_Pos"  = 18  # filled diamond
  )) +

  # 4. X-axis labels (Neg 1, Pool 1, Pos, Neg 2, Pool 2, Pos)
  scale_x_discrete(labels = c("Neg 1", "Pool 1", "Pos", "Neg 2", "Pool 2", "Pos")) +

  # 5. Significance markers (Neg vs Pool)
  # Group A: X=1 vs X=2
  annotate("segment", x = 1, xend = 2, y = y_line_A, yend = y_line_A, linewidth = 0.6) +
  annotate("text", x = 1.5, y = y_text_A, label = label_A, size = 6, family = "Helvetica", fontface = "bold") +

  # Group B: X=4 vs X=5
  annotate("segment", x = 4, xend = 5, y = y_line_B, yend = y_line_B, linewidth = 0.6) +
  annotate("text", x = 4.5, y = y_text_B, label = label_B, size = 6, family = "Helvetica", fontface = "bold") +

  # 6. Theme
  labs(x = "", y = expression("IFN-"*gamma*" spots/2x10"^5*" cells")) +
  # Set Y-axis upper limit, leaving room for Pos
  scale_y_continuous(limits = c(0, 3600), expand = c(0, 0)) +
  theme_classic() +
  theme(
    legend.position = "none", # hide legend (labels are already on the X axis)
    axis.text.x = element_text(size = 10, family = "Helvetica", face = "bold", color = "black", angle = 45, hjust = 1), # angled labels
    axis.text.y = element_text(size = 10, family = "Helvetica", face = "bold", color = "black"),
    axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold", color = "black"),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    axis.line.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank()
  ) +

  # 7. Y-axis scale break
  # Break interval 1200 - 2400, compression 0.4
  scale_y_break(c(1200, 2400), scales = 0.4, ticklabels = c(2500, 3000, 3500))

# ==============================================================================
# 6. Save figures
# ==============================================================================
output_dir <- "/path/to/project/results_V3/summary/elispot/pep_21" # adjust path as needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(file.path(output_dir, "Custom_Elispot_Final.pdf"), p_final, width = 3, height = 2.5, device = "pdf")
ggsave(file.path(output_dir, "Custom_Elispot_Final.svg"), p_final, width = 3, height = 2.5, device = "svg")

print("Plotting complete.")
print(p_final)
