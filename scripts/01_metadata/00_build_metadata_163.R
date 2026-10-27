#!/usr/bin/env Rscript

base <- "manuscript_introgression_2026"

old_meta <- file.path(base, "01_metadata", "metadata_master.tsv")
removed_file <- file.path(base, "01_metadata", "metadata_removed_samples.tsv")
sample_file <- file.path(base, "publication_data_2026", "02_metadata",
                         "final_163_vcf_samples.txt")
outdir <- file.path(base, "publication_data_2026", "02_metadata")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta <- read.table(
  old_meta,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)

removed <- read.table(
  removed_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)

vcf_ids <- readLines(sample_file)

if (length(vcf_ids) != 163)
  stop("Expected 163 VCF entries; found ", length(vcf_ids))

# ------------------------------------------------------------
# Add DsubpRef -> DsubPool2
# ------------------------------------------------------------

dsub2 <- data.frame(
  SampleVCF   = "DsubpRef",
  SampleLabel = "DsubPool2",
  PopID       = "DsubPool2",
  Species     = "Drosophila subpulchrella",
  Group       = "Native",
  Type        = "PoolSeq",
  Role        = "Dsub",
  TableRef    = "TableS1",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Add D-bia -> DbiarPool
# ------------------------------------------------------------

dbia_pool <- data.frame(
  SampleVCF   = "D-bia",
  SampleLabel = "DbiarPool",
  PopID       = "DbiarPool",
  Species     = "Drosophila biarmipes",
  Group       = "Outgroup",
  Type        = "PoolSeq",
  Role        = "Outgroup",
  TableRef    = "TableS1",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Rename existing Dbiar reference representation explicitly
# ------------------------------------------------------------

ii <- meta$SampleVCF == "DbiarRef"

if (sum(ii) != 1)
  stop("Expected exactly one DbiarRef in existing metadata.")

meta$SampleLabel[ii] <- "DbiarRef"
meta$PopID[ii]       <- "DbiarRef"

# ------------------------------------------------------------
# Combine and order EXACTLY as VCF
# ------------------------------------------------------------

newmeta <- rbind(meta, dsub2, dbia_pool)

idx <- match(vcf_ids, newmeta$SampleVCF)

if (any(is.na(idx))) {
  stop(
    "VCF samples missing from metadata: ",
    paste(vcf_ids[is.na(idx)], collapse = ", ")
  )
}

newmeta <- newmeta[idx, , drop = FALSE]

if (!identical(newmeta$SampleVCF, vcf_ids))
  stop("Metadata ordering does not match VCF.")

if (anyDuplicated(newmeta$SampleVCF))
  stop("Duplicated SampleVCF.")

if (anyDuplicated(newmeta$SampleLabel))
  stop("Duplicated SampleLabel.")

# ------------------------------------------------------------
# Expected sample composition
# ------------------------------------------------------------

stopifnot(nrow(newmeta) == 163)
stopifnot(sum(newmeta$Species == "Drosophila suzukii") == 148)
stopifnot(sum(newmeta$Species == "Drosophila subpulchrella") == 13)
stopifnot(sum(newmeta$Species == "Drosophila biarmipes") == 2)

stopifnot(sum(newmeta$Type == "PoolSeq") == 33)
stopifnot(sum(newmeta$Type == "IndSeq") == 129)
stopifnot(sum(newmeta$Type == "Reference") == 1)

# ------------------------------------------------------------
# Population plan
# ------------------------------------------------------------

popids <- unique(newmeta$PopID)

popplan <- do.call(
  rbind,
  lapply(popids, function(p) {

    x <- newmeta[newmeta$PopID == p, , drop = FALSE]

    if (length(unique(x$Species)) != 1)
      stop("Multiple species in PopID: ", p)

    if (length(unique(x$Type)) != 1)
      stop("Multiple types in PopID: ", p)

    if (length(unique(x$Role)) != 1)
      stop("Multiple roles in PopID: ", p)

    if (length(unique(x$Group)) != 1)
      stop("Multiple groups in PopID: ", p)

    data.frame(
      PopID = p,
      Species = unique(x$Species),
      Group = unique(x$Group),
      Type = unique(x$Type),
      Role = unique(x$Role),
      Nsamples = nrow(x),
      stringsAsFactors = FALSE
    )
  })
)

dsuz_pops <- popplan$PopID[popplan$Role == "Dsuz"]
dsub_pops <- popplan$PopID[popplan$Role == "Dsub"]
out_pops  <- popplan$PopID[popplan$Role == "Outgroup"]

stopifnot(length(dsuz_pops) == 44)

stopifnot(
  setequal(
    dsub_pops,
    c("Dsub_CN", "Dsub_KO", "DsubPool", "DsubPool2")
  )
)

stopifnot(
  setequal(
    out_pops,
    c("DbiarRef", "DbiarPool")
  )
)

stopifnot(nrow(popplan) == 50)

# ------------------------------------------------------------
# Write publication metadata
# ------------------------------------------------------------

write.table(
  newmeta,
  file.path(outdir, "metadata_master_163.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  newmeta[, c("SampleVCF", "SampleLabel", "PopID")],
  file.path(outdir, "metadata_name_crosswalk_163.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  popplan,
  file.path(outdir, "metadata_population_plan_50pops.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n=============================================\n")
cat("FINAL 163-ENTRY METADATA\n")
cat("=============================================\n")
cat("VCF entries:          ", nrow(newmeta), "\n")
cat("D. suzukii entries:   ",
    sum(newmeta$Role == "Dsuz"), "\n")
cat("D. subpulchrella:     ",
    sum(newmeta$Role == "Dsub"), "\n")
cat("D. biarmipes:         ",
    sum(newmeta$Role == "Outgroup"), "\n")
cat("PoolSeq entries:      ",
    sum(newmeta$Type == "PoolSeq"), "\n")
cat("IndSeq entries:       ",
    sum(newmeta$Type == "IndSeq"), "\n")
cat("Reference entries:    ",
    sum(newmeta$Type == "Reference"), "\n")
cat("Population units:     ", nrow(popplan), "\n")
cat("Dsub populations:     ",
    paste(dsub_pops, collapse = ", "), "\n")
cat("Outgroup populations: ",
    paste(out_pops, collapse = ", "), "\n")
cat("=============================================\n")
