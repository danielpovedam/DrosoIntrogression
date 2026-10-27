#!/usr/bin/env Rscript

# ============================================================
# 10_build_F3_publication_tables.R
#
# Final publication tables for revised genome-wide F3 analyses
#
# Outputs:
#   Table S3  Complete primary Dsuz-target F3 tests
#   Table S4  Primary F3 sensitivity analyses
#   Table S5  Complete reciprocal population-level F3 tests
#   Table S6  Individual-level reciprocal Dsub sensitivity summary
#   Data S1   Complete individual-level reciprocal F3 tests
#   captions  Final manuscript-ready captions
# ============================================================

rm(list = ls())

BASE <- "/home/dpoveda/work/micetral/DSU/manuscript_introgression_2026"

F3DIR <- file.path(
  BASE,
  "final_results_2026",
  "Fstats"
)

INPUT_S3 <- file.path(
  F3DIR,
  "02_F3_population",
  "F3_Dsuz_targets_Dsub_affinity_50pop.tsv"
)

INPUT_S4 <- file.path(
  F3DIR,
  "03_F3_sensitivity",
  "F3_primary_sensitivity_summary.tsv"
)

INPUT_S5 <- file.path(
  F3DIR,
  "02_F3_population",
  "F3_Dsub_targets_Dsuz_affinity_50pop.tsv"
)

INPUT_S6 <- file.path(
  F3DIR,
  "04_F3_Dsub_individual",
  "F3_Dsub_individual_target_summary.tsv"
)

INPUT_DATA_S1 <- file.path(
  F3DIR,
  "04_F3_Dsub_individual",
  "F3_Dsub_individual_targets_Dsuz_affinity.tsv"
)

OUTDIR <- file.path(
  BASE,
  "final_results_2026",
  "publication_tables",
  "F3"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Helpers
# ============================================================

check_file <- function(x) {
  if (!file.exists(x)) {
    stop("Missing input file: ", x)
  }
}

for (f in c(
  INPUT_S3,
  INPUT_S4,
  INPUT_S5,
  INPUT_S6,
  INPUT_DATA_S1
)) {
  check_file(f)
}

read_tsv <- function(x) {
  read.delim(
    x,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
}

clean_seq <- function(x) {
  x[x == "PoolSeq"] <- "Pool-Seq"
  x[x == "IndSeq"]  <- "Ind-Seq"
  x
}

fmt_f3 <- function(x) {
  formatC(
    as.numeric(x),
    format = "f",
    digits = 7
  )
}

fmt_z <- function(x) {
  formatC(
    as.numeric(x),
    format = "f",
    digits = 2
  )
}

fmt_pct <- function(x) {
  formatC(
    as.numeric(x),
    format = "f",
    digits = 2
  )
}

sig_label <- function(x) {
  ifelse(
    as.numeric(x) < -1.65,
    "Yes",
    "No"
  )
}

# ============================================================
# TABLE S3
# Complete primary F3 tests
# ============================================================

s3 <- read_tsv(INPUT_S3)

required_s3 <- c(
  "Target",
  "SameSpeciesSource",
  "OtherSpeciesSource",
  "TargetType",
  "Estimate",
  "bjack s.e.",
  "Z-score"
)

if (!all(required_s3 %in% names(s3))) {
  stop("Table S3 input missing required columns.")
}

if (nrow(s3) != 7568) {
  stop(
    "Table S3 expected 7,568 rows; found ",
    nrow(s3)
  )
}

table_s3 <- data.frame(
  `Focal D. suzukii population` =
    s3$Target,

  Sequencing =
    clean_seq(s3$TargetType),

  `D. suzukii reference` =
    s3$SameSpeciesSource,

  `D. subpulchrella reference` =
    s3$OtherSpeciesSource,

  `Reference populations` =
    paste(
      s3$SameSpeciesSource,
      s3$OtherSpeciesSource,
      sep = ", "
    ),

  F3 =
    fmt_f3(s3$Estimate),

  `Jackknife SE` =
    fmt_f3(s3[["bjack s.e."]]),

  `Z-score` =
    fmt_z(s3[["Z-score"]]),

  `Significant (Z < -1.65)` =
    sig_label(s3[["Z-score"]]),

  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Publication ordering:
# focal population alphabetically, strongest negative Z first.
ord_s3 <- order(
  table_s3[["Focal D. suzukii population"]],
  as.numeric(table_s3[["Z-score"]])
)

table_s3 <- table_s3[
  ord_s3,
  ,
  drop = FALSE
]

write.table(
  table_s3,
  file.path(
    OUTDIR,
    "Table_S3_complete_primary_F3.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# TABLE S4
# Primary F3 sensitivity analyses
# ============================================================

s4 <- read_tsv(INPUT_S4)

if (nrow(s4) != 8) {
  stop(
    "Expected eight sensitivity scenarios; found ",
    nrow(s4)
  )
}

scenario_labels <- c(
  Full =
    "Full dataset",

  Exclude_CN_Nin =
    "Excluding CN-Nin",

  Wild_Dsub_only =
    "Wild D. subpulchrella only",

  PoolSeq_Dsub_only =
    "Pool-Seq D. subpulchrella only",

  Exclude_DsubPool2 =
    "Excluding DsubPool2",

  DsubPool_only =
    "DsubPool only",

  DsubPool2_only =
    "DsubPool2 only",

  Exclude_CN_Nin_Wild_Dsub =
    "Excluding CN-Nin + wild D. subpulchrella only"
)

if (!all(s4$Scenario %in% names(scenario_labels))) {
  stop("Unexpected sensitivity scenario names.")
}

table_s4 <- data.frame(
  `Sensitivity analysis` =
    unname(
      scenario_labels[s4$Scenario]
    ),

  `D. suzukii focal populations` =
    s4$N_Dsuz_focal_populations,

  `D. subpulchrella references` =
    s4$N_Dsub_references,

  `N tests` =
    s4$N_tests,

  `Significant F3 tests` =
    paste0(
      s4$N_significant,
      " (",
      fmt_pct(s4$Percent_significant),
      "%)"
    ),

  `Focal populations with signal` =
    paste0(
      s4$N_focal_populations_with_signal,
      "/",
      s4$N_Dsuz_focal_populations,
      " (",
      fmt_pct(
        s4$Percent_focal_populations_with_signal
      ),
      "%)"
    ),

  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Preserve analytical scenario order.
desired_order <- names(scenario_labels)

table_s4 <- table_s4[
  match(
    unname(scenario_labels[desired_order]),
    table_s4[["Sensitivity analysis"]]
  ),
  ,
  drop = FALSE
]

write.table(
  table_s4,
  file.path(
    OUTDIR,
    "Table_S4_primary_F3_sensitivity.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# TABLE S5
# Complete reciprocal population-level F3 tests
# ============================================================

s5 <- read_tsv(INPUT_S5)

required_s5 <- c(
  "Target",
  "SameSpeciesSource",
  "OtherSpeciesSource",
  "TargetType",
  "Estimate",
  "bjack s.e.",
  "Z-score"
)

if (!all(required_s5 %in% names(s5))) {
  stop("Table S5 input missing required columns.")
}

if (nrow(s5) != 528) {
  stop(
    "Table S5 expected 528 rows; found ",
    nrow(s5)
  )
}

table_s5 <- data.frame(
  `Focal D. subpulchrella population` =
    s5$Target,

  Sequencing =
    clean_seq(s5$TargetType),

  `D. subpulchrella reference` =
    s5$SameSpeciesSource,

  `D. suzukii reference` =
    s5$OtherSpeciesSource,

  `Reference populations` =
    paste(
      s5$SameSpeciesSource,
      s5$OtherSpeciesSource,
      sep = ", "
    ),

  F3 =
    fmt_f3(s5$Estimate),

  `Jackknife SE` =
    fmt_f3(s5[["bjack s.e."]]),

  `Z-score` =
    fmt_z(s5[["Z-score"]]),

  `Significant (Z < -1.65)` =
    sig_label(s5[["Z-score"]]),

  check.names = FALSE,
  stringsAsFactors = FALSE
)

ord_s5 <- order(
  table_s5[["Focal D. subpulchrella population"]],
  as.numeric(table_s5[["Z-score"]])
)

table_s5 <- table_s5[
  ord_s5,
  ,
  drop = FALSE
]

write.table(
  table_s5,
  file.path(
    OUTDIR,
    "Table_S5_complete_reciprocal_F3.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# TABLE S6
# Individual-level reciprocal sensitivity summary
# ============================================================

s6 <- read_tsv(INPUT_S6)

if (nrow(s6) != 11) {
  stop(
    "Table S6 expected 11 individuals; found ",
    nrow(s6)
  )
}

table_s6 <- data.frame(
  `Focal D. subpulchrella individual` =
    s6$TargetIndividual,

  Population =
    s6$TargetPopulation,

  `N tests` =
    s6$N_tests,

  `Negative F3 tests` =
    s6$N_negative,

  `Significant F3 tests` =
    paste0(
      s6$N_significant_Zlt_minus1.65,
      " (",
      fmt_pct(s6$Percent_significant),
      "%)"
    ),

  `FDR-significant tests` =
    s6$N_significant_FDR05,

  `Minimum Z-score` =
    fmt_z(s6$Minimum_Z),

  `F3 at minimum Z` =
    fmt_f3(s6$F3_at_minimum_Z),

  `D. subpulchrella reference at minimum Z` =
    s6$Dsub_source_at_minimum_Z,

  `D. suzukii reference at minimum Z` =
    s6$Dsuz_source_at_minimum_Z,

  check.names = FALSE,
  stringsAsFactors = FALSE
)

table_s6 <- table_s6[
  order(
    table_s6$Population,
    table_s6[["Focal D. subpulchrella individual"]]
  ),
  ,
  drop = FALSE
]

write.table(
  table_s6,
  file.path(
    OUTDIR,
    "Table_S6_individual_Dsub_reciprocal_F3_sensitivity.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# DATA S1
# Complete individual-level reciprocal F3 configurations
# ============================================================

ds1 <- read_tsv(INPUT_DATA_S1)

if (nrow(ds1) != 4840) {
  stop(
    "Data S1 expected 4,840 rows; found ",
    nrow(ds1)
  )
}

required_ds1 <- c(
  "TargetIndividual",
  "TargetPopulation",
  "DsubSourceIndividual",
  "DsubSourcePopulation",
  "DsuzSourcePopulation",
  "Estimate",
  "bjack s.e.",
  "Z-score"
)

if (!all(required_ds1 %in% names(ds1))) {
  stop("Data S1 input missing required columns.")
}

data_s1 <- data.frame(
  `Focal D. subpulchrella individual` =
    ds1$TargetIndividual,

  `Focal population` =
    ds1$TargetPopulation,

  `D. subpulchrella reference individual` =
    ds1$DsubSourceIndividual,

  `D. subpulchrella reference population` =
    ds1$DsubSourcePopulation,

  `D. suzukii reference population` =
    ds1$DsuzSourcePopulation,

  F3 =
    fmt_f3(ds1$Estimate),

  `Jackknife SE` =
    fmt_f3(ds1[["bjack s.e."]]),

  `Z-score` =
    fmt_z(ds1[["Z-score"]]),

  `Significant (Z < -1.65)` =
    sig_label(ds1[["Z-score"]]),

  check.names = FALSE,
  stringsAsFactors = FALSE
)

data_s1 <- data_s1[
  order(
    data_s1[["Focal population"]],
    data_s1[["Focal D. subpulchrella individual"]],
    as.numeric(data_s1[["Z-score"]])
  ),
  ,
  drop = FALSE
]

write.table(
  data_s1,
  file.path(
    OUTDIR,
    "Data_S1_complete_individual_Dsub_reciprocal_F3.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# CAPTIONS
# ============================================================

captions <- c(

"Table S3. Complete results of the primary F3-based admixture tests in Drosophila suzukii. Results are reported for all 7,568 configurations of the form F3(Dsuz_i; Dsuz_j, Dsub_k), where Dsuz_i is the focal D. suzukii population, Dsuz_j is another D. suzukii reference population, and Dsub_k is one of the four D. subpulchrella reference populations (Dsub_CN, Dsub_KO, DsubPool, or DsubPool2). Sequencing indicates whether the focal population was represented by Pool-Seq or aggregated individual whole-genome sequencing (Ind-Seq) data. F3 estimates, block-jackknife standard errors, and Z-scores are reported for each configuration. Significantly negative F3 statistics were defined using the one-sided 95% threshold Z < -1.65.",

"Table S4. Sensitivity analyses of the primary F3-based admixture tests. Results are shown for the complete dataset and for analyses excluding CN-Nin, restricting D. subpulchrella references to the two wild populations (Dsub_CN and Dsub_KO), restricting references to the two Pool-Seq samples (DsubPool and DsubPool2), excluding DsubPool2, using either Pool-Seq D. subpulchrella representation individually, or jointly excluding CN-Nin and restricting the analysis to wild D. subpulchrella. Significantly negative F3 statistics were defined as Z < -1.65.",

"Table S5. Complete results of the reciprocal F3-based admixture tests with Drosophila subpulchrella as the focal population. Results are reported for all 528 configurations of the form F3(Dsub_i; Dsub_j, Dsuz_k), where Dsub_i is the focal D. subpulchrella population, Dsub_j is another D. subpulchrella reference population, and Dsuz_k is a D. suzukii reference population. Sequencing indicates whether the focal D. subpulchrella representation was derived from Pool-Seq or aggregated individual whole-genome sequencing (Ind-Seq) data. F3 estimates, block-jackknife standard errors, and Z-scores are reported for each configuration. Significantly negative F3 statistics were defined using the one-sided 95% threshold Z < -1.65.",

"Table S6. Individual-level sensitivity analysis of reciprocal F3 signals in wild Drosophila subpulchrella. The 11 wild Ind-Seq D. subpulchrella genomes were analysed separately as focal units against the other wild D. subpulchrella individuals and 44 population-level D. suzukii reference units. The analysis was based on 7,708,587 SNPs retained across the harmonized dataset. For each focal D. subpulchrella individual, the table reports the number and proportion of significantly negative F3 configurations, the minimum Z-score, and the reference populations associated with the strongest signal. Significantly negative F3 statistics were defined as Z < -1.65.",

"Data S1. Complete results of the individual-level reciprocal F3 sensitivity analysis. Results are reported for all 4,840 configurations of the form F3(DsubIndividual_i; DsubIndividual_j, Dsuz_k), based on 11 wild D. subpulchrella individuals and 44 population-level D. suzukii reference units. F3 estimates, block-jackknife standard errors, and Z-scores are reported for each configuration. Significantly negative F3 statistics were defined using the one-sided 95% threshold Z < -1.65."
)

writeLines(
  captions,
  file.path(
    OUTDIR,
    "F3_supplementary_table_captions.txt"
  )
)

# ============================================================
# QC SUMMARY
# ============================================================

qc <- data.frame(
  Output = c(
    "Table_S3",
    "Table_S4",
    "Table_S5",
    "Table_S6",
    "Data_S1"
  ),
  Rows = c(
    nrow(table_s3),
    nrow(table_s4),
    nrow(table_s5),
    nrow(table_s6),
    nrow(data_s1)
  ),
  Expected = c(
    7568,
    8,
    528,
    11,
    4840
  ),
  stringsAsFactors = FALSE
)

if (!all(qc$Rows == qc$Expected)) {
  stop("Final publication-table row-count QC failed.")
}

write.table(
  qc,
  file.path(
    OUTDIR,
    "F3_publication_tables_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n============================================================\n")
cat("FINAL F3 PUBLICATION TABLES COMPLETE\n")
cat("============================================================\n\n")

print(qc, row.names = FALSE)

cat("\nFiles written to:\n")
cat(OUTDIR, "\n\n")

print(
  list.files(
    OUTDIR,
    full.names = FALSE
  )
)

cat("\nDONE\n")
