#!/bin/bash
#SBATCH --job-name=DSU_F3_DsubInd
#SBATCH --chdir=/home/dpoveda/work/micetral/DSU
#SBATCH --output=/home/dpoveda/work/micetral/DSU/manuscript_introgression_2026/final_results_2026/logs/F3_DsubInd_%j.out
#SBATCH --error=/home/dpoveda/work/micetral/DSU/manuscript_introgression_2026/final_results_2026/logs/F3_DsubInd_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=180G

set -euo pipefail

module purge
module load devel/Miniconda/Miniconda3
source activate bio-env

BASE=/home/dpoveda/work/micetral/DSU
PROJ=${BASE}/manuscript_introgression_2026

SAMPLE=${PROJ}/publication_data_2026/03_poolfstat/sample_level/pooldata_163samples_autosomes_MAF001.rds

POP=${PROJ}/publication_data_2026/03_poolfstat/population_level/pooldata_50pop_autosomes_MAF001.rds

META=${PROJ}/publication_data_2026/02_metadata/metadata_master_163.tsv

PLAN=${PROJ}/publication_data_2026/02_metadata/metadata_population_plan_50pops.tsv

SCRIPT=${PROJ}/final_results_2026/scripts/09_F3_Dsub_individual_reciprocal_sensitivity.R

OUT=${PROJ}/final_results_2026/Fstats/04_F3_Dsub_individual

mkdir -p "${OUT}" "${PROJ}/final_results_2026/logs"

Rscript "${SCRIPT}" \
  "${SAMPLE}" \
  "${POP}" \
  "${META}" \
  "${PLAN}" \
  "${OUT}"

echo
echo "========================================"
echo "FINAL OUTPUTS"
echo "========================================"

ls -lh "${OUT}"

echo
cat "${OUT}/F3_Dsub_individual_analysis_summary.tsv"

echo
echo "Finished: $(date)"
