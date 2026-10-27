#!/bin/bash
#SBATCH --job-name=DSU_pool163
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=180G
#SBATCH --output=/home/dpoveda/work/micetral/DSU/manuscript_introgression_2026/final_results_2026/logs/pool163_%j.out
#SBATCH --error=/home/dpoveda/work/micetral/DSU/manuscript_introgression_2026/final_results_2026/logs/pool163_%j.err

set -euo pipefail

module load devel/Miniconda/Miniconda3
source activate bio-env

BASE=/home/dpoveda/work/micetral/DSU/manuscript_introgression_2026

VCF=${BASE}/publication_data_2026/01_vcf/DSU_final_autosomes_163entries.vcf.gz
META=${BASE}/publication_data_2026/02_metadata/metadata_master_163.tsv

OUTDIR=${BASE}/publication_data_2026/03_poolfstat/sample_level
SCRIPT=${BASE}/final_results_2026/scripts/03_build_sample_pooldata_163.R

mkdir -p "${OUTDIR}"
mkdir -p "${BASE}/final_results_2026/logs"

OUTRDS=${OUTDIR}/pooldata_163samples_autosomes_MAF001.rds
SUMMARY=${OUTDIR}/final_dataset_summary_163.tsv

echo "============================================================"
echo "BUILD 163-SAMPLE AUTOSOMAL POOLFSTAT OBJECT"
echo "============================================================"
echo "Started: $(date)"
echo "VCF:     ${VCF}"
echo "Metadata:${META}"
echo

Rscript "${SCRIPT}" \
  "${VCF}" \
  "${META}" \
  "${OUTRDS}" \
  "${SUMMARY}"

echo
echo "Finished: $(date)"
