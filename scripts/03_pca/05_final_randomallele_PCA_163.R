#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(poolfstat)
  library(ggplot2)
  library(patchwork)
})

# ============================================================
# FINAL RANDOM-ALLELE PCA
# Drosophila suzukii / D. subpulchrella introgression project
#
# Frozen input:
#   163 entries
#   8,004,305 autosomal SNPs
#
# PCA A = all 163 entries
# PCA B = exclude both D. biarmipes outgroups (161)
# PCA C = D. suzukii only (148)
# PCA D = D. suzukii excluding CN-Nin (147)
#
# IMPORTANT:
# CN-Nin is excluded ONLY in PCA D as a sensitivity analysis.
# It remains in the frozen dataset and downstream analyses.
#
# randomallele.pca() is stochastic, so fixed seeds are used.
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(
    "Usage: Rscript 02_final_randomallele_PCA.R ",
    "<pooldata.rds> <metadata.tsv> <results_dir> <figure_dir>"
  )
}

pool_rds <- args[1]
meta_file <- args[2]
outdir <- args[3]
figdir <- args[4]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("FINAL RANDOM-ALLELE PCA\n")
cat("============================================================\n")
cat("poolfstat version: ",
    as.character(packageVersion("poolfstat")), "\n")
cat("Input pooldata:   ", pool_rds, "\n")
cat("Metadata:         ", meta_file, "\n\n")

# ============================================================
# 1. Load frozen pooldata and metadata
# ============================================================

pool_all <- readRDS(pool_rds)

if (!isS4(pool_all)) {
  stop("Input RDS is not a pooldata S4 object.")
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

required_cols <- c(
  "SampleVCF",
  "SampleLabel",
  "PopID",
  "Species",
  "Group",
  "Type",
  "Role"
)

if (!all(required_cols %in% colnames(meta))) {
  stop(
    "Metadata missing: ",
    paste(
      required_cols[!required_cols %in% colnames(meta)],
      collapse = ", "
    )
  )
}

# Match metadata to RDS using manuscript SampleLabel
idx <- match(pool_all@poolnames, meta$SampleLabel)

if (any(is.na(idx))) {
  stop(
    "Pooldata names absent from metadata: ",
    paste(pool_all@poolnames[is.na(idx)], collapse = ", ")
  )
}

meta <- meta[idx, , drop = FALSE]

if (!identical(meta$SampleLabel, pool_all@poolnames)) {
  stop("Metadata and pooldata could not be aligned.")
}

# Frozen-dataset checks
if (pool_all@npools != 163) {
  stop("Expected 163 entries; found ", pool_all@npools)
}

if (pool_all@nsnp != 8004305) {
  stop(
    "Expected 8,004,305 SNPs; found ",
    pool_all@nsnp
  )
}

cat("Frozen dataset checks: PASSED\n")
cat("Entries: ", pool_all@npools, "\n")
cat("SNPs:    ",
    format(pool_all@nsnp, big.mark = ","), "\n\n")

# ============================================================
# 2. Manuscript display metadata
# ============================================================

meta$SpeciesDisplay <- meta$Species

meta$SpeciesDisplay[
  meta$Species == "Drosophila suzukii"
] <- "D. suzukii"

meta$SpeciesDisplay[
  meta$Species == "Drosophila subpulchrella"
] <- "D. subpulchrella"

meta$SpeciesDisplay[
  meta$Species == "Drosophila biarmipes"
] <- "D. biarmipes"


meta$DataTypeDisplay <- meta$Type

meta$DataTypeDisplay[
  meta$Type == "PoolSeq"
] <- "Pool-Seq"

meta$DataTypeDisplay[
  meta$Type == "IndSeq"
] <- "Ind-Seq"

meta$DataTypeDisplay[
  meta$Type == "Reference"
] <- "Reference"


meta$RegionDisplay <- meta$Group

meta$RegionDisplay[
  meta$Group == "Native"
] <- "Native Asia"

meta$RegionDisplay[
  meta$Group == "NorthAmerica"
] <- "North America"

meta$RegionDisplay[
  meta$Group == "SouthAmerica"
] <- "South America"

# DsubPool lab strain isn't used in region-coloured Dsuz panels.
meta$RegionDisplay[
  meta$Group == "LabStrain"
] <- "Laboratory strain"

# ============================================================
# 3. Publication colour and shape definitions
# ============================================================

# Species
species_cols <- c(
  "D. suzukii"       = "#1F5A94",
  "D. subpulchrella" = "#D67A2D",
  "D. biarmipes"     = "#4D4D4D"
)

# Geographic regions used for D. suzukii
region_cols <- c(
  "Native Asia"   = "#1B9E77",
  "Europe"        = "#3B6FB6",
  "North America" = "#D89028",
  "South America" = "#8A5AA5"
)

# Sequencing strategy
datatype_shapes <- c(
  "Pool-Seq"  = 21,
  "Ind-Seq"   = 22,
  "Reference" = 24
)

# ============================================================
# 4. Common publication theme
# ============================================================

theme_manuscript <- theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(
      size = 11,
      colour = "black"
    ),
    axis.text = element_text(
      size = 9,
      colour = "black"
    ),
    axis.line = element_line(
      linewidth = 0.45,
      colour = "black"
    ),
    axis.ticks = element_line(
      linewidth = 0.4,
      colour = "black"
    ),
    legend.title = element_text(
      size = 10,
      face = "plain"
    ),
    legend.text = element_text(
      size = 9
    ),
    legend.key.height = grid::unit(0.42, "cm"),
    legend.key.width = grid::unit(0.42, "cm"),
    plot.margin = margin(7, 8, 7, 8)
  )

# ============================================================
# 5. PCA A
# ALL 163 entries
# ============================================================

cat("============================================================\n")
cat("PCA A: ALL ENTRIES\n")
cat("============================================================\n")

set.seed(20260911)

pca_A <- randomallele.pca(
  pool_all,
  plot.pcs = NULL
)

scores_A <- as.data.frame(pca_A$pop.loadings)
colnames(scores_A) <- paste0(
  "PC",
  seq_len(ncol(scores_A))
)

scores_A <- cbind(
  meta,
  scores_A
)

var_A <- pca_A$perc.var

write.table(
  scores_A,
  file.path(outdir, "PCA_A_ALL_163_scores.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    PC = paste0("PC", seq_along(var_A)),
    PercentVariance = var_A
  ),
  file.path(outdir, "PCA_A_ALL_163_variance.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  pca_A,
  file.path(outdir, "PCA_A_ALL_163.rds")
)

cat(
  "n = ", nrow(scores_A),
  " | PC1 = ", round(var_A[1], 3), "%",
  " | PC2 = ", round(var_A[2], 3), "%\n\n",
  sep = ""
)

# ============================================================
# 6. PCA B
# Exclude both D. biarmipes representations
# ============================================================

cat("============================================================\n")
cat("PCA B: NO OUTGROUP\n")
cat("============================================================\n")

idx_B <- which(meta$Role != "Outgroup")

if (length(idx_B) != 161) {
  stop(
    "PCA B expected 161 entries; found ",
    length(idx_B)
  )
}

pool_B <- pooldata.subset(
  pool_all,
  pool.index = idx_B
)

meta_B <- meta[idx_B, , drop = FALSE]

set.seed(20260912)

pca_B <- randomallele.pca(
  pool_B,
  plot.pcs = NULL
)

scores_B <- as.data.frame(pca_B$pop.loadings)
colnames(scores_B) <- paste0(
  "PC",
  seq_len(ncol(scores_B))
)

scores_B <- cbind(
  meta_B,
  scores_B
)

var_B <- pca_B$perc.var

write.table(
  scores_B,
  file.path(outdir, "PCA_B_NO_OUTGROUP_161_scores.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    PC = paste0("PC", seq_along(var_B)),
    PercentVariance = var_B
  ),
  file.path(outdir, "PCA_B_NO_OUTGROUP_161_variance.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  pca_B,
  file.path(outdir, "PCA_B_NO_OUTGROUP_161.rds")
)

cat(
  "n = ", nrow(scores_B),
  " | PC1 = ", round(var_B[1], 3), "%",
  " | PC2 = ", round(var_B[2], 3), "%\n\n",
  sep = ""
)

rm(pool_B)
gc()

# ============================================================
# 7. PCA C
# D. suzukii only, including CN-Nin
# ============================================================

cat("============================================================\n")
cat("PCA C: D. SUZUKII ONLY\n")
cat("============================================================\n")

idx_C <- which(meta$Role == "Dsuz")

if (length(idx_C) != 148) {
  stop(
    "PCA C expected 148 entries; found ",
    length(idx_C)
  )
}

pool_C <- pooldata.subset(
  pool_all,
  pool.index = idx_C
)

meta_C <- meta[idx_C, , drop = FALSE]

set.seed(20260913)

pca_C <- randomallele.pca(
  pool_C,
  plot.pcs = NULL
)

scores_C <- as.data.frame(pca_C$pop.loadings)
colnames(scores_C) <- paste0(
  "PC",
  seq_len(ncol(scores_C))
)

scores_C <- cbind(
  meta_C,
  scores_C
)

var_C <- pca_C$perc.var

write.table(
  scores_C,
  file.path(outdir, "PCA_C_DSUZ_148_scores.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    PC = paste0("PC", seq_along(var_C)),
    PercentVariance = var_C
  ),
  file.path(outdir, "PCA_C_DSUZ_148_variance.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  pca_C,
  file.path(outdir, "PCA_C_DSUZ_148.rds")
)

cat(
  "n = ", nrow(scores_C),
  " | PC1 = ", round(var_C[1], 3), "%",
  " | PC2 = ", round(var_C[2], 3), "%\n\n",
  sep = ""
)

rm(pool_C)
gc()

# ============================================================
# 8. PCA D
# D. suzukii only, excluding CN-Nin
# Sensitivity analysis ONLY
# ============================================================

cat("============================================================\n")
cat("PCA D: D. SUZUKII EXCLUDING CN-Nin\n")
cat("============================================================\n")

if (!("CN-Nin" %in% meta$SampleLabel)) {
  stop("CN-Nin not found in SampleLabel.")
}

idx_D <- which(
  meta$Role == "Dsuz" &
  meta$SampleLabel != "CN-Nin"
)

if (length(idx_D) != 147) {
  stop(
    "PCA D expected 147 entries; found ",
    length(idx_D)
  )
}

pool_D <- pooldata.subset(
  pool_all,
  pool.index = idx_D
)

meta_D <- meta[idx_D, , drop = FALSE]

set.seed(20260914)

pca_D <- randomallele.pca(
  pool_D,
  plot.pcs = NULL
)

scores_D <- as.data.frame(pca_D$pop.loadings)
colnames(scores_D) <- paste0(
  "PC",
  seq_len(ncol(scores_D))
)

scores_D <- cbind(
  meta_D,
  scores_D
)

var_D <- pca_D$perc.var

write.table(
  scores_D,
  file.path(
    outdir,
    "PCA_D_DSUZ_NO_CN-Nin_147_scores.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    PC = paste0("PC", seq_along(var_D)),
    PercentVariance = var_D
  ),
  file.path(
    outdir,
    "PCA_D_DSUZ_NO_CN-Nin_147_variance.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  pca_D,
  file.path(
    outdir,
    "PCA_D_DSUZ_NO_CN-Nin_147.rds"
  )
)

cat(
  "n = ", nrow(scores_D),
  " | PC1 = ", round(var_D[1], 3), "%",
  " | PC2 = ", round(var_D[2], 3), "%\n\n",
  sep = ""
)

rm(pool_D)
gc()

# ============================================================
# 9. Publication plots
# NO panel subtitles
# ============================================================

# ----------------------------
# Panel A
# ----------------------------

pA <- ggplot(
  scores_A,
  aes(
    x = PC1,
    y = PC2,
    fill = SpeciesDisplay,
    shape = DataTypeDisplay
  )
) +
  geom_point(
    size = 2.8,
    stroke = 0.55,
    colour = "black",
    alpha = 0.90
  ) +
  scale_fill_manual(
    values = species_cols,
    name = "Species"
  ) +
  scale_shape_manual(
    values = datatype_shapes,
    name = "Data type"
  ) +
  labs(
    x = sprintf(
      "PC1 (%.1f%%)",
      var_A[1]
    ),
    y = sprintf(
      "PC2 (%.1f%%)",
      var_A[2]
    )
  ) +
  theme_manuscript

# ----------------------------
# Panel B
# ----------------------------

pB <- ggplot(
  scores_B,
  aes(
    x = PC1,
    y = PC2,
    fill = SpeciesDisplay,
    shape = DataTypeDisplay
  )
) +
  geom_point(
    size = 2.8,
    stroke = 0.55,
    colour = "black",
    alpha = 0.90
  ) +
  scale_fill_manual(
    values = species_cols[
      c(
        "D. suzukii",
        "D. subpulchrella"
      )
    ],
    name = "Species"
  ) +
  scale_shape_manual(
    values = datatype_shapes[
      c(
        "Pool-Seq",
        "Ind-Seq"
      )
    ],
    name = "Data type"
  ) +
  labs(
    x = sprintf(
      "PC1 (%.1f%%)",
      var_B[1]
    ),
    y = sprintf(
      "PC2 (%.1f%%)",
      var_B[2]
    )
  ) +
  theme_manuscript

# ----------------------------
# Panel C
# D. suzukii including CN-Nin
# ----------------------------

pC <- ggplot(
  scores_C,
  aes(
    x = PC1,
    y = PC2,
    fill = RegionDisplay,
    shape = DataTypeDisplay
  )
) +
  geom_point(
    size = 2.8,
    stroke = 0.55,
    colour = "black",
    alpha = 0.90
  ) +
  scale_fill_manual(
    values = region_cols,
    name = "Region"
  ) +
  scale_shape_manual(
    values = datatype_shapes[
      c(
        "Pool-Seq",
        "Ind-Seq"
      )
    ],
    name = "Data type"
  ) +
  labs(
    x = sprintf(
      "PC1 (%.1f%%)",
      var_C[1]
    ),
    y = sprintf(
      "PC2 (%.1f%%)",
      var_C[2]
    )
  ) +
  theme_manuscript

# ----------------------------
# Panel D
# D. suzukii excluding CN-Nin
# ----------------------------

pD <- ggplot(
  scores_D,
  aes(
    x = PC1,
    y = PC2,
    fill = RegionDisplay,
    shape = DataTypeDisplay
  )
) +
  geom_point(
    size = 2.8,
    stroke = 0.55,
    colour = "black",
    alpha = 0.90
  ) +
  scale_fill_manual(
    values = region_cols,
    name = "Region"
  ) +
  scale_shape_manual(
    values = datatype_shapes[
      c(
        "Pool-Seq",
        "Ind-Seq"
      )
    ],
    name = "Data type"
  ) +
  labs(
    x = sprintf(
      "PC1 (%.1f%%)",
      var_D[1]
    ),
    y = sprintf(
      "PC2 (%.1f%%)",
      var_D[2]
    )
  ) +
  theme_manuscript

# ============================================================
# 10. Export individual panels
# PDF, SVG, PNG 600 dpi
# ============================================================

ggsave(
  file.path(figdir, "PCA_A_ALL.pdf"),
  pA,
  width = 6.3,
  height = 5.2,
  units = "in"
)

ggsave(
  file.path(figdir, "PCA_A_ALL.svg"),
  pA,
  width = 6.3,
  height = 5.2,
  units = "in",
  device = grDevices::svg
)

ggsave(
  file.path(figdir, "PCA_A_ALL.png"),
  pA,
  width = 6.3,
  height = 5.2,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggsave(
  file.path(figdir, "PCA_B_NO_OUTGROUP.pdf"),
  pB,
  width = 6.3,
  height = 5.2,
  units = "in"
)

ggsave(
  file.path(figdir, "PCA_B_NO_OUTGROUP.svg"),
  pB,
  width = 6.3,
  height = 5.2,
  units = "in",
  device = grDevices::svg
)

ggsave(
  file.path(figdir, "PCA_B_NO_OUTGROUP.png"),
  pB,
  width = 6.3,
  height = 5.2,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggsave(
  file.path(figdir, "PCA_C_DSUZ_ALL.pdf"),
  pC,
  width = 6.3,
  height = 5.2,
  units = "in"
)

ggsave(
  file.path(figdir, "PCA_C_DSUZ_ALL.svg"),
  pC,
  width = 6.3,
  height = 5.2,
  units = "in",
  device = grDevices::svg
)

ggsave(
  file.path(figdir, "PCA_C_DSUZ_ALL.png"),
  pC,
  width = 6.3,
  height = 5.2,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggsave(
  file.path(figdir, "PCA_D_DSUZ_NO_CN-Nin.pdf"),
  pD,
  width = 6.3,
  height = 5.2,
  units = "in"
)

ggsave(
  file.path(figdir, "PCA_D_DSUZ_NO_CN-Nin.svg"),
  pD,
  width = 6.3,
  height = 5.2,
  units = "in",
  device = grDevices::svg
)

ggsave(
  file.path(figdir, "PCA_D_DSUZ_NO_CN-Nin.png"),
  pD,
  width = 6.3,
  height = 5.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ============================================================
# 11. Combined four-panel publication figure
# 2 x 2
# ONLY A/B/C/D panel labels
# ============================================================

combined_ABCD <- (
  pA + pB
) / (
  pC + pD
) +
  plot_layout(
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 14,
        colour = "black"
      )
    )
  ) &
  theme(
    legend.position = "bottom"
  )

ggsave(
  file.path(figdir, "PCA_FINAL_ABCD.pdf"),
  combined_ABCD,
  width = 12.8,
  height = 10.0,
  units = "in"
)

ggsave(
  file.path(figdir, "PCA_FINAL_ABCD.svg"),
  combined_ABCD,
  width = 12.8,
  height = 10.0,
  units = "in",
  device = grDevices::svg
)

ggsave(
  file.path(figdir, "PCA_FINAL_ABCD.png"),
  combined_ABCD,
  width = 12.8,
  height = 10.0,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ============================================================
# 12. Also produce A-B-C version
# Gives us flexibility for main/supplement later
# ============================================================

combined_ABC <- pA + pB + pC +
  plot_layout(
    nrow = 1,
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 14,
        colour = "black"
      )
    )
  ) &
  theme(
    legend.position = "bottom"
  )

ggsave(
  file.path(figdir, "PCA_FINAL_ABC.pdf"),
  combined_ABC,
  width = 16.0,
  height = 5.5,
  units = "in"
)

ggsave(
  file.path(figdir, "PCA_FINAL_ABC.svg"),
  combined_ABC,
  width = 16.0,
  height = 5.5,
  units = "in",
  device = grDevices::svg
)

ggsave(
  file.path(figdir, "PCA_FINAL_ABC.png"),
  combined_ABC,
  width = 16.0,
  height = 5.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ============================================================
# 13. Summary report
# ============================================================

summary_pca <- data.frame(
  PCA = c(
    "A_ALL",
    "B_NO_OUTGROUP",
    "C_DSUZ_ALL",
    "D_DSUZ_NO_CN-Nin"
  ),
  N_entries = c(
    nrow(scores_A),
    nrow(scores_B),
    nrow(scores_C),
    nrow(scores_D)
  ),
  PC1_percent = c(
    var_A[1],
    var_B[1],
    var_C[1],
    var_D[1]
  ),
  PC2_percent = c(
    var_A[2],
    var_B[2],
    var_C[2],
    var_D[2]
  ),
  Exclusion = c(
    "None",
    "DbiarRef + DbiarPool",
    "All non-Dsuz",
    "All non-Dsuz + CN-Nin"
  ),
  stringsAsFactors = FALSE
)

write.table(
  summary_pca,
  file.path(outdir, "PCA_FINAL_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n============================================================\n")
cat("PCA ANALYSIS COMPLETE\n")
cat("============================================================\n")

print(summary_pca)

cat("\nResults:\n", outdir, "\n")
cat("\nFigures:\n", figdir, "\n")
cat("\nCN-Nin was excluded ONLY from PCA D.\n")
cat("CN-Nin remains in the frozen analytical dataset.\n")
cat("\nDONE\n")
