#!/bin/bash
#SBATCH --job-name=DSU_poolify50
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=120G
#SBATCH --output=/home/dpoveda/work/micetral/DSU/manuscript_introgression_2026/final_results_2026/logs/poolify50_%j.out
#SBATCH --error=/home/dpoveda/work/micetral/DSU/manuscript_introgression_2026/final_results_2026/logs/poolify50_%j.err

set -euo pipefail

module load devel/Miniconda/Miniconda3
source activate bio-env

BASE=/home/dpoveda/work/micetral/DSU/manuscript_introgression_2026

INRDS=${BASE}/publication_data_2026/03_poolfstat/sample_level/pooldata_163samples_autosomes_MAF001.rds
META=${BASE}/publication_data_2026/02_metadata/metadata_master_163.tsv
PLAN=${BASE}/publication_data_2026/02_metadata/metadata_population_plan_50pops.tsv

OUTDIR=${BASE}/publication_data_2026/03_poolfstat/population_level
SCRIPT=${BASE}/final_results_2026/scripts/04_build_population_pooldata_50.R

mkdir -p "${OUTDIR}"
mkdir -p "${BASE}/final_results_2026/logs"

OUTRDS=${OUTDIR}/pooldata_50pop_autosomes_MAF001.rds
MEMBERS=${OUTDIR}/population_membership_50pop.tsv

Rscript "${SCRIPT}" \
  "${INRDS}" \
  "${META}" \
  "${PLAN}" \
  "${OUTRDS}" \
  "${MEMBERS}"
