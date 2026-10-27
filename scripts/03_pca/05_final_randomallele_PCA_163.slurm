#!/bin/bash
#SBATCH --job-name=DSU_PCA163
#SBATCH --chdir=/home/dpoveda/work/micetral/DSU
#SBATCH --output=/home/dpoveda/work/micetral/DSU/manuscript_introgression_2026/final_results_2026/logs/PCA163_%j.out
#SBATCH --error=/home/dpoveda/work/micetral/DSU/manuscript_introgression_2026/final_results_2026/logs/PCA163_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G

set -euo pipefail

module purge
module load devel/Miniconda/Miniconda3
source activate bio-env

BASE=/home/dpoveda/work/micetral/DSU
PROJ=${BASE}/manuscript_introgression_2026

POOL=${PROJ}/publication_data_2026/03_poolfstat/sample_level/pooldata_163samples_autosomes_MAF001.rds

META=${PROJ}/publication_data_2026/02_metadata/metadata_master_163.tsv

SCRIPT=${PROJ}/final_results_2026/scripts/05_final_randomallele_PCA_163.R

OUT=${PROJ}/final_results_2026/PCA
FIG=${PROJ}/final_results_2026/PCA/figures

mkdir -p "${OUT}" "${FIG}" "${PROJ}/final_results_2026/logs"

echo "============================================================"
echo "FINAL RANDOM-ALLELE PCA — 163 ENTRY DATASET"
echo "============================================================"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo

Rscript -e '
cat("R version: ", R.version.string, "\n")
cat("poolfstat: ", as.character(packageVersion("poolfstat")), "\n")
cat("ggplot2: ", as.character(packageVersion("ggplot2")), "\n")
cat("patchwork: ", as.character(packageVersion("patchwork")), "\n")
'

echo
echo "Starting PCA..."
echo

Rscript "${SCRIPT}" \
  "${POOL}" \
  "${META}" \
  "${OUT}" \
  "${FIG}"

echo
echo "============================================================"
echo "PCA SUMMARY"
echo "============================================================"

column -t -s $'\t' \
  "${OUT}/PCA_FINAL_summary.tsv"

echo
echo "Generated figures:"
ls -lh "${FIG}"

echo
echo "Finished: $(date)"
