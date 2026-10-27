#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(poolfstat))

# ============================================================
# FINAL POPULATION-LEVEL DATASET
#
# Input:
#   Frozen sample-level pooldata:
#     163 entries
#     8,004,305 SNPs
#
# Strategy:
#   - Actual Pool-Seq samples remain unchanged
#   - Dbiar remains unchanged
#   - Ind-Seq individuals are grouped by PopID with poolify()
#   - poolify() uses random-allele sampling
#   - objects are combined with data.merge()
#
# Expected final dataset:
#   44 D. suzukii populations
#    4 D. subpulchrella populations
#    2 D. biarmipes outgroups
#   --------------------------------
#   50 population units
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop(
    "Usage: Rscript 04_build_population_pooldata.R ",
    "<sample_pooldata.rds> <metadata_master.tsv> ",
    "<population_plan.tsv> <output.rds> <output_membership.tsv>"
  )
}

sample_rds   <- args[1]
meta_file    <- args[2]
plan_file    <- args[3]
out_rds      <- args[4]
out_members  <- args[5]

cat("============================================================\n")
cat("BUILD FINAL POPULATION-LEVEL POOLFSTAT OBJECT\n")
cat("============================================================\n")
cat("poolfstat version: ",
    as.character(packageVersion("poolfstat")), "\n\n")

# ------------------------------------------------------------
# Read frozen sample object
# ------------------------------------------------------------

pool <- readRDS(sample_rds)

if (!isS4(pool)) {
  stop("Input RDS is not a valid S4 pooldata object.")
}

if (pool@npools != 163) {
  stop("Expected 163 sample-level entries; found ", pool@npools)
}

if (pool@nsnp != 8004305) {
  stop("Expected 8,004,305 SNPs; found ", pool@nsnp)
}

# ------------------------------------------------------------
# Metadata
# ------------------------------------------------------------

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

if (!all(required_meta %in% colnames(meta))) {
  stop("Master metadata is incomplete.")
}

required_plan <- c(
  "PopID",
  "Species",
  "Group",
  "Type",
  "Role",
  "Nsamples"
)

if (!all(required_plan %in% colnames(plan))) {
  stop("Population plan is incomplete.")
}

# ------------------------------------------------------------
# Align metadata to frozen RDS
# ------------------------------------------------------------

idx <- match(pool@poolnames, meta$SampleLabel)

if (any(is.na(idx))) {
  stop(
    "RDS names absent from metadata: ",
    paste(pool@poolnames[is.na(idx)], collapse = ", ")
  )
}

meta <- meta[idx, , drop = FALSE]

if (!identical(pool@poolnames, meta$SampleLabel)) {
  stop("Metadata and sample-level pooldata do not align.")
}

# ------------------------------------------------------------
# Expected composition
# ------------------------------------------------------------

if (nrow(plan) != 50) {

  stop("Expected 50 population units; found ", nrow(plan))

}

if (sum(plan$Role == "Dsuz") != 44) {

  stop("Expected 44 D. suzukii populations.")

}

if (sum(plan$Role == "Dsub") != 4) {

  stop("Expected 4 D. subpulchrella populations.")

}

if (

  !setequal(

    plan$PopID[plan$Role == "Dsub"],

    c("Dsub_KO", "Dsub_CN", "DsubPool", "DsubPool2")

  )

) {

  stop("D. subpulchrella population names are incorrect.")

}

if (sum(plan$Role == "Outgroup") != 2) {

  stop("Expected 2 D. biarmipes outgroup representations.")

}

if (

  !setequal(

    plan$PopID[plan$Role == "Outgroup"],

    c("DbiarRef", "DbiarPool")

  )

) {

  stop("D. biarmipes outgroup population names are incorrect.")

}

# ============================================================
# REAL POOL-SEQ + OUTGROUP
# ============================================================

# Retain actual Pool-Seq and reference entries unchanged.
actual_idx <- which(
  meta$Type %in% c("PoolSeq", "Reference")
)

cat("Actual Pool-Seq/reference units: ",
    length(actual_idx), "\n")

if (length(actual_idx) != 34) {
  stop(
    "Expected 34 actual PoolSeq/reference units; found ",
    length(actual_idx)
  )
}

pop_final <- pooldata.subset(
  pool,
  pool.index = actual_idx
)

# Replace sample labels by authoritative population IDs.
pop_final@poolnames <- meta$PopID[actual_idx]

if (anyDuplicated(pop_final@poolnames)) {
  stop("Duplicated PopID among actual PoolSeq/reference units.")
}

# ============================================================
# IND-SEQ POPULATIONS
# ============================================================

ind_meta <- meta[
  meta$Type == "IndSeq",
  ,
  drop = FALSE
]

ind_popids <- unique(ind_meta$PopID)

cat("Ind-Seq populations to poolify: ",
    length(ind_popids), "\n")

if (length(ind_popids) != 16) {
  stop(
    "Expected 16 Ind-Seq populations ",
    "(14 Dsuz + 2 Dsub); found ",
    length(ind_popids)
  )
}

# Fixed seed for exact reproducibility.
BASE_SEED <- 20260911

for (k in seq_along(ind_popids)) {

  popid <- ind_popids[k]

  sample_names <- ind_meta$SampleLabel[
    ind_meta$PopID == popid
  ]

  sample_idx <- match(
    sample_names,
    pool@poolnames
  )

  if (any(is.na(sample_idx))) {
    stop(
      "Samples missing for population ",
      popid
    )
  }

  if (length(sample_idx) < 2) {
    stop(
      "poolify() requires >=2 samples: ",
      popid
    )
  }

  cat(
    "poolify ",
    popid,
    " : n=",
    length(sample_idx),
    " | seed=",
    BASE_SEED + k,
    "\n",
    sep = ""
  )

  set.seed(BASE_SEED + k)

  pseudo <- poolify(
    x = pool,
    sample.index = sample_idx,
    out.samplename = popid
  )

  # Merge the countdata pseudopool with the existing pooldata.
  # Explicitly use the poolfstat default fake.pool.size.
  pop_final <- data.merge(
    pop_final,
    pseudo,
    fake.pool.size = 1e6,
    verbose = FALSE
  )

  rm(pseudo)
  gc()
}

# ============================================================
# FINAL ORDER
# ============================================================

if (!isS4(pop_final)) {
  stop("Final population object is not S4.")
}

if (pop_final@npools != 50) {
  stop(
    "Expected 50 population units; found ",
    pop_final@npools
  )
}

if (pop_final@nsnp != 8004305) {
  stop(
    "SNP universe changed during poolification/merging: ",
    pop_final@nsnp
  )
}

final_order <- match(
  plan$PopID,
  pop_final@poolnames
)

if (any(is.na(final_order))) {
  stop(
    "Final object missing PopIDs: ",
    paste(
      plan$PopID[is.na(final_order)],
      collapse = ", "
    )
  )
}

pop_final <- pooldata.subset(
  pop_final,
  pool.index = final_order
)

if (!identical(
  pop_final@poolnames,
  plan$PopID
)) {
  stop("Final population ordering failed.")
}

# ============================================================
# MEMBERSHIP REPORT
# ============================================================

membership <- meta[
  ,
  c(
    "SampleVCF",
    "SampleLabel",
    "PopID",
    "Species",
    "Group",
    "Type",
    "Role"
  )
]

write.table(
  membership,
  out_members,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Save population metadata beside RDS
popmeta_file <- sub(
  "\\.rds$",
  ".metadata.tsv",
  out_rds
)

write.table(
  plan,
  popmeta_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  pop_final,
  out_rds
)

cat("\n============================================================\n")
cat("FINAL POPULATION DATASET\n")
cat("============================================================\n")
cat("Class:            ", class(pop_final), "\n")
cat("Population units: ", pop_final@npools, "\n")
cat("SNPs:             ",
    format(pop_final@nsnp, big.mark = ","), "\n")
cat("D. suzukii:       ", sum(plan$Role == "Dsuz"), "\n")
cat("D. subpulchrella: ", sum(plan$Role == "Dsub"), "\n")
cat("Outgroup:         ", sum(plan$Role == "Outgroup"), "\n")

cat("\nPopulation names:\n")
cat(paste(pop_final@poolnames, collapse = "\n"))
cat("\n\nSaved:\n", out_rds, "\n")
cat("\nDONE\n")
