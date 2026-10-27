#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(poolfstat))

# ============================================================
# Build definitive SAMPLE-LEVEL poolfstat object
#
# Input:
#   - final 161-entry autosomal VCF
#   - frozen metadata_master.tsv
#
# Steps:
#   1. Match metadata to VCF by SampleVCF (NOT by row order)
#   2. Use manuscript SampleLabel as pooldata names
#   3. Set haploid size:
#        IndSeq    = 2
#        PoolSeq   = 100
#        Reference = 100
#   4. Apply global MAF > 0.01 with vcf2pooldata()
#   5. Save definitive sample-level RDS
#   6. Export exact final SNP/sample summary
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(
    "Usage:\n",
    "Rscript 01_build_sample_pooldata.R ",
    "<final.vcf.gz> <metadata_master.tsv> ",
    "<output.rds> <summary.tsv>"
  )
}

vcf_file    <- args[1]
meta_file   <- args[2]
out_rds     <- args[3]
summary_tsv <- args[4]

cat("============================================================\n")
cat("FINAL SAMPLE-LEVEL POOLFSTAT IMPORT\n")
cat("============================================================\n")
cat("poolfstat version: ", as.character(packageVersion("poolfstat")), "\n")
cat("VCF:               ", vcf_file, "\n")
cat("Metadata:          ", meta_file, "\n")
cat("Output RDS:        ", out_rds, "\n")
cat("Summary:           ", summary_tsv, "\n\n")

# ------------------------------------------------------------
# Check input files
# ------------------------------------------------------------

if (!file.exists(vcf_file)) {
  stop("VCF not found: ", vcf_file)
}

if (!file.exists(meta_file)) {
  stop("Metadata not found: ", meta_file)
}

# ------------------------------------------------------------
# Read master metadata
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

required_cols <- c(
  "SampleVCF",
  "SampleLabel",
  "PopID",
  "Species",
  "Group",
  "Type",
  "Role",
  "TableRef"
)

if (!all(required_cols %in% colnames(meta))) {
  stop(
    "Metadata missing required columns: ",
    paste(
      required_cols[!required_cols %in% colnames(meta)],
      collapse = ", "
    )
  )
}

if (nrow(meta) != 163) {
  stop("Expected 163 metadata rows; found ", nrow(meta))
}

if (anyDuplicated(meta$SampleVCF)) {
  stop("Duplicated SampleVCF values in metadata.")
}

if (anyDuplicated(meta$SampleLabel)) {
  stop("Duplicated SampleLabel values in metadata.")
}

# ------------------------------------------------------------
# Read actual VCF sample names
# ------------------------------------------------------------

vcf_samples <- system2(
  "bcftools",
  args = c("query", "-l", vcf_file),
  stdout = TRUE,
  stderr = TRUE
)

if (length(vcf_samples) != 163) {
  stop(
    "Expected 163 samples in VCF; found ",
    length(vcf_samples)
  )
}

if (anyDuplicated(vcf_samples)) {
  stop("Duplicated sample IDs in VCF.")
}

# ------------------------------------------------------------
# Align metadata to VCF by SampleVCF
# ------------------------------------------------------------

idx <- match(vcf_samples, meta$SampleVCF)

if (any(is.na(idx))) {
  stop(
    "VCF samples absent from metadata: ",
    paste(vcf_samples[is.na(idx)], collapse = ", ")
  )
}

meta_vcf <- meta[idx, , drop = FALSE]

if (!identical(meta_vcf$SampleVCF, vcf_samples)) {
  stop("Metadata could not be aligned exactly to VCF order.")
}

cat("VCF/metadata identity check: PASSED\n")
cat("Samples: ", nrow(meta_vcf), "\n\n")

# ------------------------------------------------------------
# Haploid sample sizes
#
# Individual genomes represent diploid individuals -> 2
# PoolSeq/reference entries -> 100
#
# Actual PoolSeq read information remains in the VCF.
# randomallele.pca subsequently samples one allele/read per
# sample and SNP.
# ------------------------------------------------------------

poolsizes <- rep(100, nrow(meta_vcf))
poolsizes[meta_vcf$Type == "IndSeq"] <- 2

if (any(!meta_vcf$Type %in% c("PoolSeq", "IndSeq", "Reference"))) {
  stop(
    "Unexpected Type values: ",
    paste(
      unique(
        meta_vcf$Type[
          !meta_vcf$Type %in% c("PoolSeq", "IndSeq", "Reference")
        ]
      ),
      collapse = ", "
    )
  )
}

cat("Haploid-size assignments:\n")
cat("  IndSeq = 2:    ", sum(meta_vcf$Type == "IndSeq"), "\n")
cat("  PoolSeq = 100: ", sum(meta_vcf$Type == "PoolSeq"), "\n")
cat("  Reference=100: ", sum(meta_vcf$Type == "Reference"), "\n\n")

# ------------------------------------------------------------
# Number of variant records before final MAF recalculation
# ------------------------------------------------------------

input_n_snp <- suppressWarnings(
  as.numeric(
    system2(
      "bcftools",
      args = c("index", "-n", vcf_file),
      stdout = TRUE
    )
  )
)

cat("Input VCF variant records: ", format(input_n_snp, big.mark = ","), "\n")

# ------------------------------------------------------------
# Import VCF and apply FINAL global MAF > 0.01
#
# The VCF has already undergone:
#   QUAL
#   NS
#   depth
#   strand/read-position
#   biallelic SNP
#   autosomal
#   SNP-density masking
#
# Here we recalculate MAF on exactly the frozen 163 entries.
# ------------------------------------------------------------

cat("\nImporting VCF with vcf2pooldata()...\n")
cat("Applying global min.maf = 0.01\n\n")

pool_all <- vcf2pooldata(
  vcf.file = vcf_file,
  poolsizes = poolsizes,
  poolnames = meta_vcf$SampleLabel,
  min.cov.per.pool = -1,
  min.rc = 1,
  max.cov.per.pool = 1e6,
  min.maf = 0.01,
  remove.indels = FALSE,
  min.dist.from.indels = 0,
  nlines.per.readblock = 200000,
  verbose = TRUE
)

# ------------------------------------------------------------
# Integrity checks
# ------------------------------------------------------------

if (!isS4(pool_all)) {
  stop("vcf2pooldata did not return an S4 pooldata object.")
}

if (pool_all@npools != 163) {
  stop(
    "Expected pooldata npools = 163; found ",
    pool_all@npools
  )
}

if (!identical(pool_all@poolnames, meta_vcf$SampleLabel)) {
  stop("pooldata names do not match manuscript SampleLabel order.")
}

if (pool_all@nsnp == 0) {
  stop("Final object contains zero SNPs.")
}

# ------------------------------------------------------------
# Save definitive sample-level object
# ------------------------------------------------------------

dir.create(
  dirname(out_rds),
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(pool_all, out_rds)

# Save the exact metadata order represented in the RDS
meta_out <- paste0(
  sub("\\.rds$", "", out_rds),
  ".metadata.tsv"
)

write.table(
  meta_vcf,
  meta_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Save SNP information retained after MAF filtering
snpinfo_out <- paste0(
  sub("\\.rds$", "", out_rds),
  ".snp_info.tsv"
)

write.table(
  as.data.frame(pool_all@snp.info),
  snpinfo_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Dataset summary
# ------------------------------------------------------------

final_n_snp <- pool_all@nsnp

removed_maf <- input_n_snp - final_n_snp

retained_pct <- 100 * final_n_snp / input_n_snp

summary_df <- data.frame(
  Metric = c(
    "poolfstat_version",
    "input_vcf_entries",
    "final_entries",
    "Dsuz_entries",
    "Dsub_entries",
    "Dbiar_entries",
    "PoolSeq_entries",
    "IndSeq_entries",
    "Reference_entries",
    "input_SNPs_before_final_MAF",
    "final_SNPs_MAF_gt_0.01",
    "SNPs_removed_by_final_MAF",
    "percent_SNPs_retained"
  ),
  Value = c(
    as.character(packageVersion("poolfstat")),
    length(vcf_samples),
    pool_all@npools,
    sum(meta_vcf$Role == "Dsuz"),
    sum(meta_vcf$Role == "Dsub"),
    sum(meta_vcf$Role == "Outgroup"),
    sum(meta_vcf$Type == "PoolSeq"),
    sum(meta_vcf$Type == "IndSeq"),
    sum(meta_vcf$Type == "Reference"),
    input_n_snp,
    final_n_snp,
    removed_maf,
    retained_pct
  ),
  stringsAsFactors = FALSE
)

write.table(
  summary_df,
  summary_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n============================================================\n")
cat("FINAL DATASET SUMMARY\n")
cat("============================================================\n")

cat("Sample-level entries:       ", pool_all@npools, "\n")
cat("D. suzukii:                 ", sum(meta_vcf$Role == "Dsuz"), "\n")
cat("D. subpulchrella:           ", sum(meta_vcf$Role == "Dsub"), "\n")
cat("D. biarmipes:               ", sum(meta_vcf$Role == "Outgroup"), "\n")

cat("\nInput SNPs:                  ",
    format(input_n_snp, big.mark = ","), "\n")

cat("Final SNPs after MAF >0.01: ",
    format(final_n_snp, big.mark = ","), "\n")

cat("SNPs removed:               ",
    format(removed_maf, big.mark = ","), "\n")

cat("SNPs retained (%):          ",
    round(retained_pct, 3), "\n")

cat("\nSaved RDS:\n", out_rds, "\n")
cat("\nSaved metadata:\n", meta_out, "\n")
cat("\nSaved SNP information:\n", snpinfo_out, "\n")
cat("\nSaved summary:\n", summary_tsv, "\n")

cat("\nDONE\n")
