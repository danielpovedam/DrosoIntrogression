#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(poolfstat)
})

# ============================================================
# 05_compute_all_fstats.R
#
# Drosophila suzukii / D. subpulchrella introgression project
#
# Genome-wide F-statistics on the FINAL population-level dataset
#
# Frozen input:
#   50 population units
#   44 D. suzukii
#    4 D. subpulchrella
#    2 D. biarmipes outgroup representations
#   8,004,305 autosomal SNPs
#
# Analyses:
#   - all F3 configurations
#   - all F4 configurations
#   - block-jackknife SE and Z-scores
#
# Following:
#   Gautier et al. 2022
#   Gay et al. 2026
#
# Block size:
#   200,000 consecutive SNPs
#
# IMPORTANT:
#   - D-statistics are NOT calculated here
#   - Q covariance matrix is NOT calculated
#   - F2 block samples are NOT retained
#
# Expected:
#   F3 = 58,800 configurations
#   F4 = 676,900 configurations
# ============================================================


# ============================================================
# 1. COMMAND-LINE ARGUMENTS
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {

  stop(
    paste0(
      "\nUsage:\n",
      "Rscript 05_compute_all_fstats.R ",
      "<population_pooldata.rds> ",
      "<output_fstats.rds> ",
      "<output_f3.tsv.gz> ",
      "<output_f4.tsv.gz> ",
      "<summary.tsv>\n"
    )
  )
}

pop_rds     <- args[1]
out_rds     <- args[2]
out_f3      <- args[3]
out_f4      <- args[4]
summary_tsv <- args[5]


# ============================================================
# 2. FIXED ANALYSIS PARAMETERS
# ============================================================

NSNP_BLOCK <- 200000L

EXPECTED_NPOPS <- 50L
EXPECTED_NSNP  <- 8004305L

EXPECTED_DSUZ <- 44L
EXPECTED_DSUB <- 4L
EXPECTED_OUT  <- 2L

SPECIAL_POPS <- c(
  "Dsub_KO",
  "Dsub_CN",
  "DsubPool",
  "DsubPool2",
  "DbiarRef",
  "DbiarPool"
)


# ============================================================
# 3. START REPORT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("GENOME-WIDE F3 AND F4 STATISTICS\n")
cat("============================================================\n")

cat(
  "Date:              ",
  format(Sys.time()),
  "\n"
)

cat(
  "poolfstat version: ",
  as.character(packageVersion("poolfstat")),
  "\n"
)

cat(
  "Input RDS:         ",
  pop_rds,
  "\n"
)

cat(
  "Jackknife block:   ",
  format(NSNP_BLOCK, big.mark = ","),
  " consecutive SNPs\n\n"
)


# ============================================================
# 4. INPUT CHECK
# ============================================================

if (!file.exists(pop_rds)) {

  stop(
    "Population-level RDS not found: ",
    pop_rds
  )
}

pop <- readRDS(pop_rds)

if (!isS4(pop)) {

  stop(
    "Input object is not an S4 object."
  )
}

if (!inherits(pop, "pooldata")) {

  stop(
    "Input object is not of class pooldata. Class = ",
    paste(class(pop), collapse = ", ")
  )
}


# ============================================================
# 5. VALIDATE FROZEN DATASET
# ============================================================

cat("Input dataset:\n")
cat(
  "  Populations: ",
  pop@npools,
  "\n"
)

cat(
  "  SNPs:        ",
  format(pop@nsnp, big.mark = ","),
  "\n\n"
)


if (pop@npools != EXPECTED_NPOPS) {

  stop(
    "Expected ",
    EXPECTED_NPOPS,
    " populations but found ",
    pop@npools
  )
}


if (pop@nsnp != EXPECTED_NSNP) {

  stop(
    "Expected ",
    format(EXPECTED_NSNP, big.mark = ","),
    " SNPs but found ",
    format(pop@nsnp, big.mark = ",")
  )
}


if (anyDuplicated(pop@poolnames)) {

  dups <- unique(
    pop@poolnames[
      duplicated(pop@poolnames)
    ]
  )

  stop(
    "Duplicated population names: ",
    paste(dups, collapse = ", ")
  )
}


missing_special <- SPECIAL_POPS[
  !SPECIAL_POPS %in% pop@poolnames
]

if (length(missing_special) > 0) {

  stop(
    "Missing required populations: ",
    paste(
      missing_special,
      collapse = ", "
    )
  )
}


cat("Required focal populations present:\n")

for (p in SPECIAL_POPS) {

  cat(
    "  ",
    p,
    "\n",
    sep = ""
  )
}

cat("\n")


# ============================================================
# 6. EXPECTED NUMBER OF CONFIGURATIONS
# ============================================================

n <- pop@npools

expected_f3 <- (
  n *
  (n - 1) *
  (n - 2)
) / 2


expected_f4 <- (
  n *
  (n - 1) *
  (n - 2) *
  (n - 3)
) / 8


# Approximate number of consecutive-SNP blocks.
# compute.fstats() performs the actual internal block assignment.

approx_nblocks <- ceiling(
  pop@nsnp /
    NSNP_BLOCK
)


cat("Expected calculations:\n")

cat(
  "  F3 configurations: ",
  format(
    expected_f3,
    big.mark = ","
  ),
  "\n"
)

cat(
  "  F4 configurations: ",
  format(
    expected_f4,
    big.mark = ","
  ),
  "\n"
)

cat(
  "  Approx. blocks:     ",
  approx_nblocks,
  "\n\n"
)


if (expected_f3 != 58800) {

  stop(
    "Internal F3-count check failed."
  )
}


if (expected_f4 != 676900) {

  stop(
    "Internal F4-count check failed."
  )
}


# ============================================================
# 7. MEMORY CLEANUP BEFORE ANALYSIS
# ============================================================

gc(verbose = TRUE)


# ============================================================
# 8. COMPUTE COMPLETE F-STATISTICS
# ============================================================

cat("============================================================\n")
cat("RUNNING compute.fstats()\n")
cat("============================================================\n\n")

start_time <- Sys.time()


fs <- compute.fstats(

  x = pop,

  nsnp.per.bjack.block = NSNP_BLOCK,

  # Not required for the treeness analysis.
  computeDstat = FALSE,

  # Estimate all F3.
  computeF3 = TRUE,

  # Estimate all F4.
  computeF4 = TRUE,

  # Not required for this manuscript step.
  output.pairwise.fst = FALSE,

  # Not required for this manuscript step.
  output.pairwise.div = FALSE,

  # Critical for memory:
  # do NOT compute the huge F2/F3 covariance matrix.
  computeQmat = FALSE,

  # Not needed for our downstream analyses.
  return.F2.blockjackknife.samples = FALSE,

  verbose = TRUE
)


end_time <- Sys.time()

runtime_hours <- as.numeric(
  difftime(
    end_time,
    start_time,
    units = "hours"
  )
)


cat("\n")
cat("compute.fstats() finished.\n")

cat(
  "Runtime (hours): ",
  round(runtime_hours, 3),
  "\n\n"
)


# ============================================================
# 9. VALIDATE FSTATS OBJECT
# ============================================================

if (!isS4(fs)) {

  stop(
    "compute.fstats() did not return an S4 object."
  )
}


n_f3 <- nrow(
  fs@f3.values
)

n_f4 <- nrow(
  fs@f4.values
)


cat("Returned statistics:\n")

cat(
  "  F3 rows: ",
  format(
    n_f3,
    big.mark = ","
  ),
  "\n"
)

cat(
  "  F4 rows: ",
  format(
    n_f4,
    big.mark = ","
  ),
  "\n\n"
)


if (n_f3 != expected_f3) {

  stop(
    "Unexpected number of F3 configurations: ",
    n_f3,
    ". Expected ",
    expected_f3
  )
}


if (n_f4 != expected_f4) {

  stop(
    "Unexpected number of F4 configurations: ",
    n_f4,
    ". Expected ",
    expected_f4
  )
}


# ============================================================
# 10. CHECK COMPARISON TABLES
# ============================================================

if (
  is.null(
    fs@comparisons[["F3"]]
  )
) {

  stop(
    "F3 comparison table missing from fstats object."
  )
}


if (
  is.null(
    fs@comparisons[["F4"]]
  )
) {

  stop(
    "F4 comparison table missing from fstats object."
  )
}


f3_comp <- as.data.frame(
  fs@comparisons[["F3"]],
  stringsAsFactors = FALSE
)


f4_comp <- as.data.frame(
  fs@comparisons[["F4"]],
  stringsAsFactors = FALSE
)


if (
  !identical(
    colnames(f3_comp),
    c(
      "Px",
      "P1",
      "P2"
    )
  )
) {

  stop(
    "Unexpected F3 comparison columns: ",
    paste(
      colnames(f3_comp),
      collapse = ", "
    )
  )
}


if (
  !identical(
    colnames(f4_comp),
    c(
      "P1",
      "P2",
      "P3",
      "P4"
    )
  )
) {

  stop(
    "Unexpected F4 comparison columns: ",
    paste(
      colnames(f4_comp),
      collapse = ", "
    )
  )
}


if (
  nrow(f3_comp) != n_f3
) {

  stop(
    "F3 comparison/value row mismatch."
  )
}


if (
  nrow(f4_comp) != n_f4
) {

  stop(
    "F4 comparison/value row mismatch."
  )
}


# ============================================================
# 11. BUILD FULL F3 TABLE
# ============================================================

f3_values <- as.data.frame(
  fs@f3.values,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


expected_value_cols <- c(
  "Estimate",
  "bjack mean",
  "bjack s.e.",
  "Z-score"
)


if (
  !all(
    expected_value_cols %in%
      colnames(f3_values)
  )
) {

  stop(
    "Unexpected F3 value columns: ",
    paste(
      colnames(f3_values),
      collapse = ", "
    )
  )
}


f3_out <- cbind(
  f3_comp,
  f3_values
)


# Add convenience columns.
# These do not alter the poolfstat estimates.

f3_out$P_one_sided_negative <- pnorm(
  f3_out[["Z-score"]]
)


f3_out$Significant_negative_Zlt_minus1.65 <- (
  f3_out[["Z-score"]] < -1.65
)


# ============================================================
# 12. BUILD FULL F4 TABLE
# ============================================================

f4_values <- as.data.frame(
  fs@f4.values,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


if (
  !all(
    expected_value_cols %in%
      colnames(f4_values)
  )
) {

  stop(
    "Unexpected F4 value columns: ",
    paste(
      colnames(f4_values),
      collapse = ", "
    )
  )
}


f4_out <- cbind(
  f4_comp,
  f4_values
)


# Convenience variables for later treeness extraction.

f4_out$AbsZ <- abs(
  f4_out[["Z-score"]]
)


f4_out$Treeness_compatible_95 <- (
  f4_out$AbsZ < 1.96
)


# ============================================================
# 13. SAVE COMPLETE FSTATS OBJECT FIRST
#
# This is the most important output.
# All downstream extraction should use this frozen object.
# ============================================================

dir.create(
  dirname(out_rds),
  recursive = TRUE,
  showWarnings = FALSE
)


saveRDS(
  fs,
  out_rds
)


cat(
  "Saved complete fstats RDS:\n  ",
  out_rds,
  "\n\n"
)


# ============================================================
# 14. EXPORT COMPLETE F3 TABLE
# ============================================================

cat(
  "Writing complete F3 table...\n"
)


con_f3 <- gzfile(
  out_f3,
  open = "wt"
)


write.table(
  f3_out,
  file = con_f3,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


close(
  con_f3
)


cat(
  "Saved:\n  ",
  out_f3,
  "\n\n"
)


# ============================================================
# 15. EXPORT COMPLETE F4 TABLE
# ============================================================

cat(
  "Writing complete F4 table...\n"
)


con_f4 <- gzfile(
  out_f4,
  open = "wt"
)


write.table(
  f4_out,
  file = con_f4,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


close(
  con_f4
)


cat(
  "Saved:\n  ",
  out_f4,
  "\n\n"
)


# ============================================================
# 16. BASIC QC COUNTS
#
# These are descriptive only.
# Biological reciprocal subsets are extracted in script 06.
# ============================================================

n_f3_negative <- sum(
  f3_out[["Significant_negative_Zlt_minus1.65"]],
  na.rm = TRUE
)


n_f4_treeness <- sum(
  f4_out[["Treeness_compatible_95"]],
  na.rm = TRUE
)


cat("Global descriptive QC:\n")

cat(
  "  F3 with Z < -1.65: ",
  format(
    n_f3_negative,
    big.mark = ","
  ),
  "\n"
)

cat(
  "  F4 configurations with |Z| < 1.96: ",
  format(
    n_f4_treeness,
    big.mark = ","
  ),
  "\n\n"
)


# ============================================================
# 17. WRITE ANALYSIS SUMMARY
# ============================================================

summary_df <- data.frame(

  Metric = c(

    "poolfstat_version",

    "population_units",

    "autosomal_SNPs",

    "jackknife_block_size_SNPs",

    "approximate_number_blocks",

    "F3_expected",

    "F3_observed",

    "F4_expected",

    "F4_observed",

    "F3_Z_lt_minus1.65_all_configurations",

    "F4_absZ_lt_1.96_all_configurations",

    "runtime_hours"
  ),

  Value = c(

    as.character(
      packageVersion(
        "poolfstat"
      )
    ),

    pop@npools,

    pop@nsnp,

    NSNP_BLOCK,

    approx_nblocks,

    expected_f3,

    n_f3,

    expected_f4,

    n_f4,

    n_f3_negative,

    n_f4_treeness,

    runtime_hours
  ),

  stringsAsFactors = FALSE
)


write.table(

  summary_df,

  file = summary_tsv,

  sep = "\t",

  quote = FALSE,

  row.names = FALSE
)


# ============================================================
# 18. SAVE POPULATION ORDER
#
# Important because poolfstat uses population order when
# generating canonical F3/F4 configurations.
# ============================================================

population_order_file <- file.path(

  dirname(out_rds),

  "population_order_50.txt"
)


writeLines(

  pop@poolnames,

  con = population_order_file
)


# ============================================================
# 19. SAVE SESSION INFORMATION
# ============================================================

session_file <- file.path(

  dirname(out_rds),

  "05_compute_all_fstats_sessionInfo.txt"
)


sink(
  session_file
)

cat(
  "Analysis completed:\n",
  format(
    Sys.time()
  ),
  "\n\n"
)

print(
  sessionInfo()
)

sink()


# ============================================================
# 20. FINAL REPORT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("COMPLETE F-STATISTICS DATASET FINISHED\n")
cat("============================================================\n")

cat(
  "Population units: ",
  pop@npools,
  "\n"
)

cat(
  "SNPs:             ",
  format(
    pop@nsnp,
    big.mark = ","
  ),
  "\n"
)

cat(
  "F3 configurations:",
  " ",
  format(
    n_f3,
    big.mark = ","
  ),
  "\n"
)

cat(
  "F4 configurations:",
  " ",
  format(
    n_f4,
    big.mark = ","
  ),
  "\n"
)

cat(
  "Block size:       ",
  format(
    NSNP_BLOCK,
    big.mark = ","
  ),
  " SNPs\n"
)

cat(
  "Approx. blocks:   ",
  approx_nblocks,
  "\n"
)

cat(
  "Runtime:          ",
  round(
    runtime_hours,
    3
  ),
  " h\n"
)

cat("\nCore output:\n")
cat(
  "  ",
  out_rds,
  "\n"
)

cat("\nTables:\n")
cat(
  "  ",
  out_f3,
  "\n"
)
cat(
  "  ",
  out_f4,
  "\n"
)

cat("\nSummary:\n")
cat(
  "  ",
  summary_tsv,
  "\n"
)

cat("\nPopulation order:\n")
cat(
  "  ",
  population_order_file,
  "\n"
)

cat("\nSession information:\n")
cat(
  "  ",
  session_file,
  "\n"
)

cat("\nDONE\n")
