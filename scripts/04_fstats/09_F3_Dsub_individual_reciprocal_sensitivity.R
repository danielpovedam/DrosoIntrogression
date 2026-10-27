#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(poolfstat)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop(
    "Usage: Rscript 09_F3_Dsub_individual_reciprocal_sensitivity.R ",
    "<sample_pooldata.rds> <population_pooldata.rds> ",
    "<metadata_master.tsv> <population_plan.tsv> <output_directory>"
  )
}

sample_rds <- args[1]
pop_rds    <- args[2]
meta_file  <- args[3]
plan_file  <- args[4]
outdir     <- args[5]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

NSNP_BLOCK <- 200000L
Z_THRESHOLD <- -1.65

EXPECTED_NSNP <- 8004305L
EXPECTED_DSUZ <- 44L
EXPECTED_DSUB_IND <- 11L
EXPECTED_TOTAL_UNITS <- 55L
EXPECTED_TARGETED_F3 <- 4840L
EXPECTED_TESTS_PER_DSUB_IND <- 440L

cat("============================================================\n")
cat("INDIVIDUAL-LEVEL RECIPROCAL Dsub F3 SENSITIVITY\n")
cat("============================================================\n\n")

# ============================================================
# 1. INPUTS
# ============================================================

sample_pool <- readRDS(sample_rds)
pop_pool    <- readRDS(pop_rds)

if (!inherits(sample_pool, "pooldata")) {
  stop("sample_rds is not a pooldata object.")
}

if (!inherits(pop_pool, "pooldata")) {
  stop("pop_rds is not a pooldata object.")
}

if (sample_pool@nsnp != EXPECTED_NSNP) {
  stop(
    "Unexpected SNP count in sample object: ",
    sample_pool@nsnp
  )
}

if (pop_pool@nsnp != EXPECTED_NSNP) {
  stop(
    "Unexpected SNP count in population object: ",
    pop_pool@nsnp
  )
}

# The sample- and population-level objects must originate from
# exactly the same ordered SNP universe before any subsetting.
if (!identical(
  sample_pool@snp.info,
  pop_pool@snp.info
)) {
  stop(
    "Sample-level and population-level SNP information/order ",
    "are not identical."
  )
}

cat("Sample/population SNP universes are identical: PASSED\n\n")

meta <- read.table(
  meta_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)

plan <- read.table(
  plan_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)

required_meta <- c(
  "SampleVCF",
  "SampleLabel",
  "PopID",
  "Species",
  "Group",
  "Type",
  "Role"
)

required_plan <- c(
  "PopID",
  "Species",
  "Group",
  "Type",
  "Role"
)

if (!all(required_meta %in% names(meta))) {
  stop("Master metadata incomplete.")
}

if (!all(required_plan %in% names(plan))) {
  stop("Population plan incomplete.")
}

# ============================================================
# 2. 44 D. suzukii POPULATION UNITS
# ============================================================

dsuz_names <- plan$PopID[
  plan$Role == "Dsuz"
]

if (length(dsuz_names) != EXPECTED_DSUZ) {
  stop(
    "Expected 44 D. suzukii populations; found ",
    length(dsuz_names)
  )
}

idx_dsuz <- match(
  dsuz_names,
  pop_pool@poolnames
)

if (anyNA(idx_dsuz)) {
  stop(
    "Missing D. suzukii population units: ",
    paste(dsuz_names[is.na(idx_dsuz)], collapse = ", ")
  )
}

dsuz_pool <- pooldata.subset(
  pop_pool,
  pool.index = idx_dsuz
)

if (!identical(
  dsuz_pool@poolnames,
  dsuz_names
)) {
  stop("D. suzukii population order mismatch.")
}

# ============================================================
# 3. 11 WILD D. subpulchrella INDIVIDUALS
# ============================================================

dsub_meta <- meta[
  meta$Role == "Dsub" &
  meta$Type == "IndSeq",
  ,
  drop = FALSE
]

if (nrow(dsub_meta) != EXPECTED_DSUB_IND) {
  stop(
    "Expected 11 wild Dsub Ind-Seq individuals; found ",
    nrow(dsub_meta)
  )
}

# Explicitly confirm these are the two wild population groups.
if (!setequal(
  unique(dsub_meta$PopID),
  c("Dsub_CN", "Dsub_KO")
)) {
  stop(
    "Unexpected Dsub Ind-Seq populations: ",
    paste(unique(dsub_meta$PopID), collapse = ", ")
  )
}

idx_dsub <- match(
  dsub_meta$SampleLabel,
  sample_pool@poolnames
)

if (anyNA(idx_dsub)) {
  stop(
    "Missing Dsub individuals in sample-level pooldata: ",
    paste(
      dsub_meta$SampleLabel[is.na(idx_dsub)],
      collapse = ", "
    )
  )
}

dsub_pool <- pooldata.subset(
  sample_pool,
  pool.index = idx_dsub,
  return.snp.idx = TRUE
)

# Use individual sample IDs as analytical unit names.
dsub_pool@poolnames <- dsub_meta$SampleLabel

if (anyDuplicated(dsub_pool@poolnames)) {
  stop("Duplicated Dsub individual analytical names.")
}

# pooldata.subset() removes SNPs with zero total coverage across
# the selected pools. Therefore the 11 Dsub individuals may retain
# fewer SNPs than the complete 8,004,305-SNP analytical universe.
#
# Because return.snp.idx=TRUE was used above, row names identify
# each retained SNP by its index in the original sample-level
# pooldata object. The population-level object has already been
# verified to have exactly the same ordered SNP universe.
dsub_snp_idx <- as.integer(
  sub(
    "^rs",
    "",
    rownames(dsub_pool@snp.info)
  )
)

if (
  length(dsub_snp_idx) != dsub_pool@nsnp ||
  anyNA(dsub_snp_idx)
) {
  stop("Could not recover original SNP indexes from Dsub subset.")
}

cat(
  "Dsub-informative SNPs retained: ",
  format(length(dsub_snp_idx), big.mark = ","),
  "\n"
)

# Restrict the 44 D. suzukii population units to exactly the same
# SNP set before merging.
dsuz_pool <- pooldata.subset(
  pop_pool,
  pool.index = idx_dsuz,
  snp.index = dsub_snp_idx,
  return.snp.idx = TRUE,
  verbose = TRUE
)

if (dsuz_pool@nsnp != dsub_pool@nsnp) {
  stop(
    "Dsu/Dsub SNP counts remain inconsistent after harmonization: ",
    dsuz_pool@nsnp,
    " vs ",
    dsub_pool@nsnp
  )
}

if (!identical(
  rownames(dsuz_pool@snp.info),
  rownames(dsub_pool@snp.info)
)) {
  stop("Dsu and Dsub SNP index row names differ after harmonization.")
}

same_snp_content <- all(
  dsuz_pool@snp.info$Chromosome == dsub_pool@snp.info$Chromosome &
  dsuz_pool@snp.info$Position   == dsub_pool@snp.info$Position &
  dsuz_pool@snp.info$RefAllele  == dsub_pool@snp.info$RefAllele &
  dsuz_pool@snp.info$AltAllele  == dsub_pool@snp.info$AltAllele
)

if (!same_snp_content) {
  stop("Dsu and Dsub chromosome/position/allele information differs after harmonization.")
}

cat("Dsu/Dsub SNP identities and coordinates match: PASSED\n")

cat(
  "Harmonized SNP universe: ",
  format(dsuz_pool@nsnp, big.mark = ","),
  " SNPs\n\n"
)

# ============================================================
# 4. MERGE HYBRID OBJECT
# ============================================================

hybrid <- data.merge(
  dsuz_pool,
  dsub_pool,
  fake.pool.size = 1e6,
  verbose = FALSE
)

if (hybrid@npools != EXPECTED_TOTAL_UNITS) {
  stop(
    "Expected 55 hybrid analytical units; found ",
    hybrid@npools
  )
}

if (hybrid@nsnp != dsub_pool@nsnp) {
  stop(
    "Unexpected SNP count after harmonized merge: ",
    hybrid@nsnp,
    " expected ",
    dsub_pool@nsnp
  )
}

cat("Hybrid object:\n")
cat("  D. suzukii populations: ", length(dsuz_names), "\n")
cat("  Dsub individuals:       ", nrow(dsub_meta), "\n")
cat("  Total units:             ", hybrid@npools, "\n")
cat("  SNPs:                    ",
    format(hybrid@nsnp, big.mark = ","), "\n\n")

# ============================================================
# 5. HYBRID METADATA
# ============================================================

meta_dsuz <- plan[
  plan$Role == "Dsuz",
  ,
  drop = FALSE
]

hybrid_meta_dsuz <- data.frame(
  UnitID = meta_dsuz$PopID,
  SourcePopID = meta_dsuz$PopID,
  Species = meta_dsuz$Species,
  Group = meta_dsuz$Group,
  Type = meta_dsuz$Type,
  Role = "Dsuz",
  stringsAsFactors = FALSE
)

hybrid_meta_dsub <- data.frame(
  UnitID = dsub_meta$SampleLabel,
  SourcePopID = dsub_meta$PopID,
  Species = dsub_meta$Species,
  Group = dsub_meta$Group,
  Type = dsub_meta$Type,
  Role = "DsubIndividual",
  stringsAsFactors = FALSE
)

hybrid_meta <- rbind(
  hybrid_meta_dsuz,
  hybrid_meta_dsub
)

idx_meta <- match(
  hybrid@poolnames,
  hybrid_meta$UnitID
)

if (anyNA(idx_meta)) {
  stop("Hybrid metadata alignment failed.")
}

hybrid_meta <- hybrid_meta[
  idx_meta,
  ,
  drop = FALSE
]

if (!identical(
  hybrid_meta$UnitID,
  hybrid@poolnames
)) {
  stop("Hybrid metadata order mismatch.")
}

write.table(
  hybrid_meta,
  file.path(
    outdir,
    "F3_individual_Dsub_hybrid_metadata.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 6. COMPUTE F3 ONLY
# ============================================================

cat("Running compute.fstats() — F3 only...\n\n")

start_time <- Sys.time()

fs <- compute.fstats(
  x = hybrid,
  nsnp.per.bjack.block = NSNP_BLOCK,
  computeDstat = FALSE,
  computeF3 = TRUE,
  computeF4 = FALSE,
  output.pairwise.fst = FALSE,
  output.pairwise.div = FALSE,
  computeQmat = FALSE,
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

saveRDS(
  fs,
  file.path(
    outdir,
    "fstats_55units_Dsub_individual_F3_block200k.rds"
  )
)

# ============================================================
# 7. EXTRACT TARGETED RECIPROCAL F3
# ============================================================

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

f3 <- cbind(
  f3_comp,
  f3_val
)

role_map <- setNames(
  hybrid_meta$Role,
  hybrid_meta$UnitID
)

pop_map <- setNames(
  hybrid_meta$SourcePopID,
  hybrid_meta$UnitID
)

f3$TargetRole <- unname(
  role_map[f3$Px]
)

f3$P1Role <- unname(
  role_map[f3$P1]
)

f3$P2Role <- unname(
  role_map[f3$P2]
)

keep <- (
  f3$TargetRole == "DsubIndividual" &
  (
    (
      f3$P1Role == "DsubIndividual" &
      f3$P2Role == "Dsuz"
    ) |
    (
      f3$P1Role == "Dsuz" &
      f3$P2Role == "DsubIndividual"
    )
  )
)

targeted <- f3[
  keep,
  ,
  drop = FALSE
]

if (nrow(targeted) != EXPECTED_TARGETED_F3) {
  stop(
    "Expected 4,840 targeted reciprocal F3 tests; found ",
    nrow(targeted)
  )
}

targeted$TargetIndividual <- targeted$Px

targeted$DsubSourceIndividual <- ifelse(
  targeted$P1Role == "DsubIndividual",
  targeted$P1,
  targeted$P2
)

targeted$DsuzSourcePopulation <- ifelse(
  targeted$P1Role == "Dsuz",
  targeted$P1,
  targeted$P2
)

targeted$TargetPopulation <- unname(
  pop_map[targeted$TargetIndividual]
)

targeted$DsubSourcePopulation <- unname(
  pop_map[targeted$DsubSourceIndividual]
)

targeted$P_one_sided_negative <- pnorm(
  targeted[["Z-score"]]
)

targeted$FDR_BH <- p.adjust(
  targeted$P_one_sided_negative,
  method = "BH"
)

targeted$Significant_Zlt_minus1.65 <- (
  targeted[["Z-score"]] < Z_THRESHOLD
)

targeted$Significant_FDR05 <- (
  targeted$Estimate < 0 &
  targeted$FDR_BH < 0.05
)

tests_per_ind <- table(
  targeted$TargetIndividual
)

if (
  length(tests_per_ind) != 11 ||
  any(tests_per_ind != EXPECTED_TESTS_PER_DSUB_IND)
) {
  stop(
    "Each Dsub individual should have exactly 440 tests."
  )
}

# ============================================================
# 8. TARGET-INDIVIDUAL SUMMARY
# ============================================================

individual_summary <- do.call(
  rbind,
  lapply(
    sort(unique(targeted$TargetIndividual)),
    function(ind) {

      z <- targeted[
        targeted$TargetIndividual == ind,
        ,
        drop = FALSE
      ]

      best <- which.min(
        z[["Z-score"]]
      )

      data.frame(
        TargetIndividual = ind,
        TargetPopulation = z$TargetPopulation[1],
        N_tests = nrow(z),
        N_negative =
          sum(z$Estimate < 0, na.rm = TRUE),
        N_significant_Zlt_minus1.65 =
          sum(
            z$Significant_Zlt_minus1.65,
            na.rm = TRUE
          ),
        Percent_significant =
          100 * mean(
            z$Significant_Zlt_minus1.65,
            na.rm = TRUE
          ),
        N_significant_FDR05 =
          sum(
            z$Significant_FDR05,
            na.rm = TRUE
          ),
        Minimum_Z =
          z[["Z-score"]][best],
        F3_at_minimum_Z =
          z$Estimate[best],
        Dsub_source_at_minimum_Z =
          z$DsubSourceIndividual[best],
        Dsuz_source_at_minimum_Z =
          z$DsuzSourcePopulation[best],
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(individual_summary) <- NULL

individual_summary$Percent_significant <-
  round(
    individual_summary$Percent_significant,
    2
  )

# ============================================================
# 9. SUMMARY BY WILD Dsub POPULATION
# ============================================================

population_summary <- do.call(
  rbind,
  lapply(
    sort(unique(targeted$TargetPopulation)),
    function(popid) {

      z <- targeted[
        targeted$TargetPopulation == popid,
        ,
        drop = FALSE
      ]

      data.frame(
        TargetPopulation = popid,
        N_target_individuals =
          length(unique(z$TargetIndividual)),
        N_tests = nrow(z),
        N_negative =
          sum(z$Estimate < 0, na.rm = TRUE),
        N_significant_Zlt_minus1.65 =
          sum(
            z$Significant_Zlt_minus1.65,
            na.rm = TRUE
          ),
        Percent_significant =
          round(
            100 * mean(
              z$Significant_Zlt_minus1.65,
              na.rm = TRUE
            ),
            2
          ),
        N_individuals_with_signal =
          length(
            unique(
              z$TargetIndividual[
                z$Significant_Zlt_minus1.65
              ]
            )
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(population_summary) <- NULL

# ============================================================
# 10. OUTPUTS
# ============================================================

write.table(
  targeted,
  file.path(
    outdir,
    "F3_Dsub_individual_targets_Dsuz_affinity.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  individual_summary,
  file.path(
    outdir,
    "F3_Dsub_individual_target_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  population_summary,
  file.path(
    outdir,
    "F3_Dsub_individual_summary_by_population.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_df <- data.frame(
  Metric = c(
    "SNPs",
    "Dsuz_population_units",
    "Dsub_individual_units",
    "hybrid_units",
    "targeted_reciprocal_F3_tests",
    "significant_Z_lt_minus1.65",
    "Dsub_individuals_with_signal",
    "runtime_hours"
  ),
  Value = c(
    hybrid@nsnp,
    44,
    11,
    hybrid@npools,
    nrow(targeted),
    sum(
      targeted$Significant_Zlt_minus1.65,
      na.rm = TRUE
    ),
    length(
      unique(
        targeted$TargetIndividual[
          targeted$Significant_Zlt_minus1.65
        ]
      )
    ),
    runtime_hours
  ),
  stringsAsFactors = FALSE
)

write.table(
  summary_df,
  file.path(
    outdir,
    "F3_Dsub_individual_analysis_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n============================================================\n")
cat("INDIVIDUAL Dsub RECIPROCAL F3 COMPLETE\n")
cat("============================================================\n\n")

print(
  individual_summary,
  row.names = FALSE
)

cat("\nBY SOURCE POPULATION\n\n")

print(
  population_summary,
  row.names = FALSE
)

cat("\nOverall:\n")
print(summary_df, row.names = FALSE)

cat("\nDONE\n")
