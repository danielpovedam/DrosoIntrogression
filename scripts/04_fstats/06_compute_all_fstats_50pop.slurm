#!/bin/bash
#SBATCH --job-name=DSU_fstats50
#SBATCH --chdir=/home/dpoveda/work/micetral/DSU
#SBATCH --output=/home/dpoveda/work/micetral/DSU/manuscript_introgression_2026/final_results_2026/logs/fstats50_%j.out
#SBATCH --error=/home/dpoveda/work/micetral/DSU/manuscript_introgression_2026/final_results_2026/logs/fstats50_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=180G

set -euo pipefail

module purge
module load devel/Miniconda/Miniconda3
source activate bio-env

BASE=/home/dpoveda/work/micetral/DSU
PROJ=${BASE}/manuscript_introgression_2026

INPUT=${PROJ}/publication_data_2026/03_poolfstat/population_level/pooldata_50pop_autosomes_MAF001.rds
SCRIPT=${PROJ}/final_results_2026/scripts/06_compute_all_fstats_50pop.R
OUTDIR=${PROJ}/final_results_2026/Fstats/01_full_50pop

mkdir -p "${OUTDIR}" "${PROJ}/final_results_2026/logs"

OUTRDS=${OUTDIR}/fstats_50pop_block200k.rds
OUTF3=${OUTDIR}/all_F3_50pop_block200k.tsv.gz
OUTF4=${OUTDIR}/all_F4_50pop_block200k.tsv.gz
SUMMARY=${OUTDIR}/fstats_50pop_block200k_summary.tsv

echo "============================================================"
echo "COMPLETE F3/F4 ANALYSIS — 50 POPULATIONS"
echo "============================================================"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo

Rscript "${SCRIPT}" \
  "${INPUT}" \
  "${OUTRDS}" \
  "${OUTF3}" \
  "${OUTF4}" \
  "${SUMMARY}"

echo
column -t -s $'\t' "${SUMMARY}"

echo
echo "Output files:"
ls -lh "${OUTDIR}"

echo
echo "Finished: $(date)"
