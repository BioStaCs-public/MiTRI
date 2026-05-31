# Test whether per-window coverage over candidate microbial genomes differs
# between responder (R) and non-responder (NR) groups, using the mapped
# bedgraph windows produced by bedgraph_significance.sh.

# Microbial species to test
microorganisms <- c("Fusobacterium_nucleatum", "Fusobacterium_animalis", "Fusobacterium_canifelinum",
                    "Fusobacterium_hwasookii", "Fusobacterium_massiliense", "Fusobacterium_polymorphum",
                    "Fusobacterium_pseudoperiodonticum", "Fusobacterium_vincentii", "Fusobacterium_periodonticum")

# Input / output directory
output_dir <- "/path/to/project/results/special/CRC/bedgraph_signif"
R_dir <- output_dir   # responder-group files
NR_dir <- output_dir  # non-responder-group files

# Initialise an empty data frame for the results
results <- data.frame(Microbe = character(),
                      TTest_p_value = numeric(),
                      Wilcoxon_p_value = numeric(),
                      TTest_one_sided_p_value = numeric(),
                      Wilcoxon_one_sided_p_value = numeric(),
                      stringsAsFactors = FALSE)

# Process each microbe
for (microbe in microorganisms) {
  # Read responder-group data
  R_file <- paste0(R_dir, "/", microbe, "_R_mapped.bed")
  R_data <- read.table(R_file, header = FALSE)
  R_cov <- R_data$V4  # column 4 holds coverage

  # Read non-responder-group data
  NR_file <- paste0(NR_dir, "/", microbe, "_NR_mapped.bed")
  NR_data <- read.table(NR_file, header = FALSE)
  NR_cov <- NR_data$V4  # column 4 holds coverage

  # Combine R and NR coverage
  cov_data <- data.frame(
    Coverage = c(R_cov, NR_cov),
    Group = rep(c("R", "NR"), times = c(length(R_cov), length(NR_cov)))
  )

  # Set Group as a factor so that "R" is the reference level
  cov_data$Group <- factor(cov_data$Group, levels = c("R", "NR"))

  # Two-sided t-test
  t_test_result <- t.test(Coverage ~ Group, data = cov_data)

  # One-sided t-test (is R significantly higher than NR|)
  t_test_one_sided_result <- t.test(Coverage ~ Group, data = cov_data, alternative = "greater")

  # Mann-Whitney U test (non-parametric)
  wilcoxon_result <- wilcox.test(Coverage ~ Group, data = cov_data)

  # One-sided Mann-Whitney U test (is R significantly higher than NR|)
  wilcoxon_one_sided_result <- wilcox.test(Coverage ~ Group, data = cov_data, alternative = "greater")

  # Store the result
  results <- rbind(results, data.frame(
    Microbe = microbe,
    # TTest_p_value = t_test_result$p.value,
    Wilcoxon_p_value = wilcoxon_result$p.value
    # TTest_one_sided_p_value = t_test_one_sided_result$p.value,
    # Wilcoxon_one_sided_p_value = wilcoxon_one_sided_result$p.value
  ))
}

# Inspect the results
print(results)

# Save the results as CSV
write.csv(results, file = "/path/to/project/results/special/CRC/bedgraph_signif/microbial_coverage_results.csv", row.names = FALSE)
