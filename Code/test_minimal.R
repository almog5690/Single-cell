#!/usr/bin/env Rscript
# Minimal test script - runs quickly on one tissue with mean expression only
# Tests the new do.interaction parameter

# Set environment before sourcing config
run.env <- "or.wsl"  # Change to "or.windows" if running in Windows R
data.type <- "TM.facs"

setwd("/mnt/c/Code/Github/Single-cell/Code")
source("scRNA_seq_config.R")

# Test 1: Basic regression on Liver only (fast)
cat("\n=== Test 1: Basic mean regression on Liver ===\n")
result <- expression_regression_analysis(
  data.type = "TM.facs",
  expression.stat.y = "mean",
  feature.types = "gene.len",
  reg.params = list(fc.method = "log_old_minus_young"),
  filter.params = filter.params,
  force.rerun = TRUE,
  run.organs = "Liver"
)
cat("Result dimensions:", dim(result), "\n")
cat("Columns:", paste(head(names(result), 10), collapse=", "), "...\n")
print(result[, c("Organs", "Cell_type", "gene.len_young_cor", "gene.len_old_cor")])

# Test 2: With interaction (tests new feature)
cat("\n=== Test 2: With do.interaction=TRUE ===\n")
result_inter <- expression_regression_analysis(
  data.type = "TM.facs",
  expression.stat.y = "mean",
  feature.types = "gene.len",
  reg.params = list(fc.method = "log_old_minus_young"),
  filter.params = filter.params,
  force.rerun = TRUE,
  run.organs = "Liver",
  do.interaction = TRUE
)
cat("Result dimensions:", dim(result_inter), "\n")
# Check interaction columns exist
inter_cols <- grep("_beta_Old|_beta_Inter|_beta_Young", names(result_inter), value=TRUE)
cat("Interaction columns found:", length(inter_cols), "\n")
if(length(inter_cols) > 0) {
  print(result_inter[, c("Organs", "Cell_type", inter_cols)])
  cat("\n=== SUCCESS: Interaction analysis working! ===\n")
} else {
  cat("\n=== WARNING: No interaction columns found ===\n")
}
