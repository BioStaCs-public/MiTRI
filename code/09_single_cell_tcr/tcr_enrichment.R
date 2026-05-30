# clonally expanded, high-affinity TCRs (events)
counts <- c(769, 494)

# number of candidate epitopes/peptides (exposure)
exposure <- c(7101, 446)

res <- poisson.test(counts, exposure, alternative = "two.sided")

res$p.value
res$estimate   # rate ratio
res$conf.int   # 95% CI

# Panel B (top)
# events = count of high-affinity TCRs
counts_B   <- c(769, 494)     # Genome-wide, MTR-derived
# exposure = count of epitopes
exposure_B <- c(3979, 238)    # Genome-wide, MTR-derived

res_B <- poisson.test(counts_B, exposure_B, alternative = "two.sided")

res_B$p.value
res_B$estimate      # rate ratio = (Genome-wide rate) / (MTR-derived rate)
res_B$conf.int

# Fold enrichment of MTR-derived relative to genome-wide (~10x)
fold_B <- 1 / as.numeric(res_B$estimate)
fold_B


# Panel D (bottom)
counts_D   <- c(413, 221)     # Genome-wide, MTR-derived
exposure_D <- c(3122, 208)    # Genome-wide, MTR-derived

res_D <- poisson.test(counts_D, exposure_D, alternative = "two.sided")

res_D$p.value
res_D$estimate
res_D$conf.int

fold_D <- 1 / as.numeric(res_D$estimate)
fold_D

