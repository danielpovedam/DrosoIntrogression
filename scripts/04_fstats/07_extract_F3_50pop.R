#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(poolfstat)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Usage: Rscript 07_extract_F3_50pop.R ",
    "<fstats_50pop.rds> <population_metadata.tsv> <output_directory>"
  )
}

fstats_file <- args[1]
meta_file   <- args[2]
outdir      <- args[3]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

EXPECTED_NPOP      <- 50L
EXPECTED_F3_TOTAL  <- 58800L
EXPECTED_F3_DSUZ   <- 7568L
EXPECTED_F3_DSUB   <- 528L
F3_Z_THRESHOLD     <- -1.65

DSUB_EXPECTED <- c(
  "Dsub_CN",
  "Dsub_KO",
  "DsubPool",
  "DsubPool2"
)

cat("============================================================\n")
cat("F3 EXTRACTION — FINAL 50-POPULATION DATASET\n")
cat("============================================================\n")

fs <- readRDS(fstats_file)

if (!inherits(fs, "fstats")) {
  stop("Input object is not of class fstats.")
}

meta <- read.table(
  meta_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)

required_meta <- c(
  "PopID",
  "Species",
  "Group",
  "Type",
  "Role"
)

if (!all(required_meta %in% names(meta))) {
  stop("Population metadata missing required columns.")
}

if (nrow(meta) != EXPECTED_NPOP) {
  stop("Expected 50 metadata rows; found ", nrow(meta))
}

if (sum(meta$Role == "Dsuz") != 44) {
  stop("Expected 44 D. suzukii populations.")
}

if (sum(meta$Role == "Dsub") != 4) {
  stop("Expected 4 D. subpulchrella populations.")
}

if (sum(meta$Role == "Outgroup") != 2) {
  stop("Expected 2 D. biarmipes outgroup representations.")
}

if (!setequal(
  meta$PopID[meta$Role == "Dsub"],
  DSUB_EXPECTED
)) {
  stop(
    "Unexpected Dsub populations: ",
    paste(meta$PopID[meta$Role == "Dsub"], collapse = ", ")
  )
}

role_map  <- setNames(meta$Role,  meta$PopID)
group_map <- setNames(meta$Group, meta$PopID)
type_map  <- setNames(meta$Type,  meta$PopID)

f3_comp <- as.data.frame(
  fs@comparisons[["F3"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

f3_val <- as.data.frame(
  fs@f3.values,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!identical(colnames(f3_comp), c("Px", "P1", "P2"))) {
  stop("Unexpected F3 comparison columns.")
}

if (nrow(f3_comp) != EXPECTED_F3_TOTAL) {
  stop(
    "Expected ", EXPECTED_F3_TOTAL,
    " complete F3 configurations; found ",
    nrow(f3_comp)
  )
}

f3 <- cbind(f3_comp, f3_val)

f3$TargetRole <- unname(role_map[f3$Px])
f3$P1Role     <- unname(role_map[f3$P1])
f3$P2Role     <- unname(role_map[f3$P2])

if (
  anyNA(f3$TargetRole) ||
  anyNA(f3$P1Role) ||
  anyNA(f3$P2Role)
) {
  stop("At least one F3 population is absent from metadata.")
}

# ============================================================
# PRIMARY F3
# F3(Dsuz_target ; Dsuz_source, Dsub_source)
# ============================================================

keep_dsuz <- (
  f3$TargetRole == "Dsuz" &
  (
    (f3$P1Role == "Dsuz" & f3$P2Role == "Dsub") |
    (f3$P1Role == "Dsub" & f3$P2Role == "Dsuz")
  )
)

f3_dsuz <- f3[keep_dsuz, , drop = FALSE]

if (nrow(f3_dsuz) != EXPECTED_F3_DSUZ) {
  stop(
    "Expected ", EXPECTED_F3_DSUZ,
    " primary Dsuz-target F3 configurations; found ",
    nrow(f3_dsuz)
  )
}

f3_dsuz$Target <- f3_dsuz$Px

f3_dsuz$SameSpeciesSource <- ifelse(
  f3_dsuz$P1Role == "Dsuz",
  f3_dsuz$P1,
  f3_dsuz$P2
)

f3_dsuz$OtherSpeciesSource <- ifelse(
  f3_dsuz$P1Role == "Dsub",
  f3_dsuz$P1,
  f3_dsuz$P2
)

f3_dsuz$TargetGroup <- unname(group_map[f3_dsuz$Target])
f3_dsuz$TargetType  <- unname(type_map[f3_dsuz$Target])

f3_dsuz$P_one_sided_negative <- pnorm(
  f3_dsuz[["Z-score"]]
)

f3_dsuz$FDR_BH <- p.adjust(
  f3_dsuz$P_one_sided_negative,
  method = "BH"
)

f3_dsuz$Significant_Zlt_minus1.65 <- (
  f3_dsuz[["Z-score"]] < F3_Z_THRESHOLD
)

f3_dsuz$Significant_FDR05 <- (
  f3_dsuz$Estimate < 0 &
  f3_dsuz$FDR_BH < 0.05
)

tests_per_dsuz <- table(f3_dsuz$Target)

if (
  length(tests_per_dsuz) != 44 ||
  any(tests_per_dsuz != 172)
) {
  stop("Expected 172 tests per D. suzukii focal population.")
}

# ============================================================
# RECIPROCAL F3
# F3(Dsub_target ; Dsub_source, Dsuz_source)
# ============================================================

keep_dsub <- (
  f3$TargetRole == "Dsub" &
  (
    (f3$P1Role == "Dsub" & f3$P2Role == "Dsuz") |
    (f3$P1Role == "Dsuz" & f3$P2Role == "Dsub")
  )
)

f3_dsub <- f3[keep_dsub, , drop = FALSE]

if (nrow(f3_dsub) != EXPECTED_F3_DSUB) {
  stop(
    "Expected ", EXPECTED_F3_DSUB,
    " reciprocal Dsub-target F3 configurations; found ",
    nrow(f3_dsub)
  )
}

f3_dsub$Target <- f3_dsub$Px

f3_dsub$SameSpeciesSource <- ifelse(
  f3_dsub$P1Role == "Dsub",
  f3_dsub$P1,
  f3_dsub$P2
)

f3_dsub$OtherSpeciesSource <- ifelse(
  f3_dsub$P1Role == "Dsuz",
  f3_dsub$P1,
  f3_dsub$P2
)

f3_dsub$TargetGroup <- unname(group_map[f3_dsub$Target])
f3_dsub$TargetType  <- unname(type_map[f3_dsub$Target])

f3_dsub$P_one_sided_negative <- pnorm(
  f3_dsub[["Z-score"]]
)

f3_dsub$FDR_BH <- p.adjust(
  f3_dsub$P_one_sided_negative,
  method = "BH"
)

f3_dsub$Significant_Zlt_minus1.65 <- (
  f3_dsub[["Z-score"]] < F3_Z_THRESHOLD
)

f3_dsub$Significant_FDR05 <- (
  f3_dsub$Estimate < 0 &
  f3_dsub$FDR_BH < 0.05
)

tests_per_dsub <- table(f3_dsub$Target)

if (
  length(tests_per_dsub) != 4 ||
  any(tests_per_dsub != 132)
) {
  stop("Expected 132 reciprocal tests per Dsub focal population.")
}

# ============================================================
# TARGET SUMMARIES
# ============================================================

summarise_targets <- function(dat) {

  spl <- split(dat, dat$Target)

  out <- lapply(spl, function(z) {

    best <- which.min(z[["Z-score"]])

    data.frame(
      Target = z$Target[1],
      TargetGroup = z$TargetGroup[1],
      TargetType = z$TargetType[1],
      N_tests = nrow(z),
      N_negative = sum(z$Estimate < 0, na.rm = TRUE),
      N_significant_Zlt_minus1.65 =
        sum(z$Significant_Zlt_minus1.65, na.rm = TRUE),
      N_significant_FDR05 =
        sum(z$Significant_FDR05, na.rm = TRUE),
      Minimum_Z = z[["Z-score"]][best],
      F3_at_minimum_Z = z$Estimate[best],
      SameSpeciesSource_at_minimum_Z =
        z$SameSpeciesSource[best],
      OtherSpeciesSource_at_minimum_Z =
        z$OtherSpeciesSource[best],
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

f3_dsuz_summary <- summarise_targets(f3_dsuz)
f3_dsub_summary <- summarise_targets(f3_dsub)

# ============================================================
# OUTPUTS
# ============================================================

write.table(
  f3_dsuz,
  file.path(outdir, "F3_Dsuz_targets_Dsub_affinity_50pop.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  f3_dsuz_summary,
  file.path(outdir, "F3_Dsuz_target_summary_50pop.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  f3_dsub,
  file.path(outdir, "F3_Dsub_targets_Dsuz_affinity_50pop.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  f3_dsub_summary,
  file.path(outdir, "F3_Dsub_target_summary_50pop.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_df <- data.frame(
  Analysis = c(
    "F3_Dsuz_targets",
    "F3_Dsub_targets"
  ),
  N_tests = c(
    nrow(f3_dsuz),
    nrow(f3_dsub)
  ),
  N_significant_Zlt_minus1.65 = c(
    sum(f3_dsuz$Significant_Zlt_minus1.65, na.rm = TRUE),
    sum(f3_dsub$Significant_Zlt_minus1.65, na.rm = TRUE)
  ),
  N_targets_with_signal = c(
    length(unique(
      f3_dsuz$Target[f3_dsuz$Significant_Zlt_minus1.65]
    )),
    length(unique(
      f3_dsub$Target[f3_dsub$Significant_Zlt_minus1.65]
    ))
  ),
  stringsAsFactors = FALSE
)

write.table(
  summary_df,
  file.path(outdir, "F3_50pop_analysis_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n============================================================\n")
cat("F3 EXTRACTION COMPLETE\n")
cat("============================================================\n\n")

print(summary_df, row.names = FALSE)

cat("\nPrimary tests per Dsuz target: 172\n")
cat("Reciprocal tests per Dsub target: 132\n")
cat("\nDONE\n")
