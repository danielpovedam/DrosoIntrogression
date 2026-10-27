#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage: Rscript 08_F3_primary_sensitivity_50pop.R ",
    "<F3_primary.tsv> <output_directory>"
  )
}

input_file <- args[1]
outdir     <- args[2]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

f3 <- read.delim(
  input_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_cols <- c(
  "Target",
  "SameSpeciesSource",
  "OtherSpeciesSource",
  "TargetGroup",
  "TargetType",
  "Estimate",
  "Z-score"
)

if (!all(required_cols %in% names(f3))) {
  stop("Missing required columns.")
}

if (nrow(f3) != 7568) {
  stop("Expected 7,568 primary F3 configurations; found ", nrow(f3))
}

DSUB_ALL <- c(
  "Dsub_CN",
  "Dsub_KO",
  "DsubPool",
  "DsubPool2"
)

if (!setequal(unique(f3$OtherSpeciesSource), DSUB_ALL)) {
  stop("Unexpected Dsub reference set.")
}

Z_THRESHOLD <- -1.65

f3$Significant <- f3[["Z-score"]] < Z_THRESHOLD

clean_seq <- function(x) {
  ifelse(
    x == "PoolSeq", "Pool-Seq",
    ifelse(x == "IndSeq", "Ind-Seq", x)
  )
}

clean_region <- function(x) {
  y <- x
  y[y %in% c("Native", "NativeAsia", "Asia")] <- "Native Asia"
  y[y %in% c("NorthAmerica", "North_America")] <- "North America"
  y[y %in% c("SouthAmerica", "South_America")] <- "South America"
  y
}

f3$Sequencing <- clean_seq(f3$TargetType)
f3$Region <- clean_region(f3$TargetGroup)

scenario_sets <- list(

  Full = list(
    remove_CN_Nin = FALSE,
    dsub = DSUB_ALL
  ),

  Exclude_CN_Nin = list(
    remove_CN_Nin = TRUE,
    dsub = DSUB_ALL
  ),

  Wild_Dsub_only = list(
    remove_CN_Nin = FALSE,
    dsub = c("Dsub_CN", "Dsub_KO")
  ),

  PoolSeq_Dsub_only = list(
    remove_CN_Nin = FALSE,
    dsub = c("DsubPool", "DsubPool2")
  ),

  Exclude_DsubPool2 = list(
    remove_CN_Nin = FALSE,
    dsub = c("Dsub_CN", "Dsub_KO", "DsubPool")
  ),

  DsubPool_only = list(
    remove_CN_Nin = FALSE,
    dsub = "DsubPool"
  ),

  DsubPool2_only = list(
    remove_CN_Nin = FALSE,
    dsub = "DsubPool2"
  ),

  Exclude_CN_Nin_Wild_Dsub = list(
    remove_CN_Nin = TRUE,
    dsub = c("Dsub_CN", "Dsub_KO")
  )
)

expected_n <- c(
  Full = 44 * 43 * 4,
  Exclude_CN_Nin = 43 * 42 * 4,
  Wild_Dsub_only = 44 * 43 * 2,
  PoolSeq_Dsub_only = 44 * 43 * 2,
  Exclude_DsubPool2 = 44 * 43 * 3,
  DsubPool_only = 44 * 43,
  DsubPool2_only = 44 * 43,
  Exclude_CN_Nin_Wild_Dsub = 43 * 42 * 2
)

summary_list <- list()
region_list  <- list()
seq_list     <- list()
dsub_list    <- list()
focal_list   <- list()

for (scenario in names(scenario_sets)) {

  pars <- scenario_sets[[scenario]]

  x <- f3[
    f3$OtherSpeciesSource %in% pars$dsub,
    ,
    drop = FALSE
  ]

  if (pars$remove_CN_Nin) {
    x <- x[
      x$Target != "CN-Nin" &
      x$SameSpeciesSource != "CN-Nin",
      ,
      drop = FALSE
    ]
  }

  if (nrow(x) != expected_n[[scenario]]) {
    stop(
      "Unexpected number of tests for ", scenario,
      ": expected ", expected_n[[scenario]],
      ", found ", nrow(x)
    )
  }

  focal_all <- unique(x$Target)
  focal_sig <- unique(x$Target[x$Significant])

  summary_list[[scenario]] <- data.frame(
    Scenario = scenario,
    N_Dsuz_focal_populations = length(focal_all),
    N_Dsub_references = length(unique(x$OtherSpeciesSource)),
    N_tests = nrow(x),
    N_significant = sum(x$Significant, na.rm = TRUE),
    Percent_significant =
      100 * mean(x$Significant, na.rm = TRUE),
    N_focal_populations_with_signal =
      length(focal_sig),
    Percent_focal_populations_with_signal =
      100 * length(focal_sig) / length(focal_all),
    stringsAsFactors = FALSE
  )

  region_list[[scenario]] <- do.call(
    rbind,
    lapply(sort(unique(x$Region)), function(reg) {

      z <- x[x$Region == reg, , drop = FALSE]

      data.frame(
        Scenario = scenario,
        Region = reg,
        N_focal_populations = length(unique(z$Target)),
        N_focal_populations_with_signal =
          length(unique(z$Target[z$Significant])),
        N_tests = nrow(z),
        N_significant = sum(z$Significant, na.rm = TRUE),
        Percent_significant =
          100 * mean(z$Significant, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )

  seq_list[[scenario]] <- do.call(
    rbind,
    lapply(sort(unique(x$Sequencing)), function(sq) {

      z <- x[x$Sequencing == sq, , drop = FALSE]

      data.frame(
        Scenario = scenario,
        Sequencing = sq,
        N_focal_populations = length(unique(z$Target)),
        N_focal_populations_with_signal =
          length(unique(z$Target[z$Significant])),
        N_tests = nrow(z),
        N_significant = sum(z$Significant, na.rm = TRUE),
        Percent_significant =
          100 * mean(z$Significant, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )

  dsub_list[[scenario]] <- do.call(
    rbind,
    lapply(sort(unique(x$OtherSpeciesSource)), function(ds) {

      z <- x[x$OtherSpeciesSource == ds, , drop = FALSE]

      data.frame(
        Scenario = scenario,
        Dsub_reference = ds,
        N_tests = nrow(z),
        N_significant = sum(z$Significant, na.rm = TRUE),
        Percent_significant =
          100 * mean(z$Significant, na.rm = TRUE),
        N_focal_populations_with_signal =
          length(unique(z$Target[z$Significant])),
        stringsAsFactors = FALSE
      )
    })
  )

  focal_list[[scenario]] <- do.call(
    rbind,
    lapply(sort(unique(x$Target)), function(target) {

      z <- x[x$Target == target, , drop = FALSE]
      best <- which.min(z[["Z-score"]])

      data.frame(
        Scenario = scenario,
        Target = target,
        Region = z$Region[1],
        Sequencing = z$Sequencing[1],
        N_tests = nrow(z),
        N_significant =
          sum(z$Significant, na.rm = TRUE),
        Percent_significant =
          100 * mean(z$Significant, na.rm = TRUE),
        Any_significant =
          any(z$Significant, na.rm = TRUE),
        Minimum_Z =
          z[["Z-score"]][best],
        F3_at_minimum_Z =
          z$Estimate[best],
        Dsuz_reference_at_minimum_Z =
          z$SameSpeciesSource[best],
        Dsub_reference_at_minimum_Z =
          z$OtherSpeciesSource[best],
        stringsAsFactors = FALSE
      )
    })
  )
}

summary_all <- do.call(rbind, summary_list)
region_all  <- do.call(rbind, region_list)
seq_all     <- do.call(rbind, seq_list)
dsub_all    <- do.call(rbind, dsub_list)
focal_all   <- do.call(rbind, focal_list)

rownames(summary_all) <- NULL
rownames(region_all)  <- NULL
rownames(seq_all)     <- NULL
rownames(dsub_all)    <- NULL
rownames(focal_all)   <- NULL

summary_all$Percent_significant <-
  round(summary_all$Percent_significant, 2)

summary_all$Percent_focal_populations_with_signal <-
  round(summary_all$Percent_focal_populations_with_signal, 2)

region_all$Percent_significant <-
  round(region_all$Percent_significant, 2)

seq_all$Percent_significant <-
  round(seq_all$Percent_significant, 2)

dsub_all$Percent_significant <-
  round(dsub_all$Percent_significant, 2)

focal_all$Percent_significant <-
  round(focal_all$Percent_significant, 2)

write.table(
  summary_all,
  file.path(outdir, "F3_primary_sensitivity_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  region_all,
  file.path(outdir, "F3_primary_sensitivity_by_region.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  seq_all,
  file.path(outdir, "F3_primary_sensitivity_by_sequencing.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  dsub_all,
  file.path(outdir, "F3_primary_sensitivity_by_Dsub_reference.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  focal_all,
  file.path(outdir, "F3_primary_sensitivity_by_focal_population.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  focal_all[focal_all$Any_significant, , drop = FALSE],
  file.path(
    outdir,
    "F3_primary_sensitivity_significant_focal_populations.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n============================================================\n")
cat("FINAL F3 PRIMARY SENSITIVITY SUMMARY\n")
cat("============================================================\n\n")

print(summary_all, row.names = FALSE)

cat("\nBY D. SUBPULCHRELLA REFERENCE\n\n")
print(dsub_all, row.names = FALSE)

cat("\nDONE\n")
