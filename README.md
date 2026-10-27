# DrosoIntrogression

Reproducible analysis workflow for the manuscript investigating historical
*Drosophila subpulchrella*-related allele sharing in *Drosophila suzukii*.

## Dataset

The revised autosomal dataset contains:

- 163 sequencing entries
- 148 *D. suzukii*
- 13 *D. subpulchrella*
- 2 *D. biarmipes* outgroup representations

Sequencing representations:

- 33 Pool-Seq
- 129 Ind-Seq
- 1 reference genome

The final autosomal input VCF contains 8,545,336 biallelic SNP positions on
chromosomes 2L, 2R, 3, and 4 after the genomic filtering and masking used for
the manuscript.

A global MAF > 0.01 filter was recalculated with poolfstat v3.1.0 on the
complete 163-entry dataset, retaining:

- 8,004,305 autosomal SNPs

## Population-level dataset

Individual-sequencing samples are combined by population using
`poolfstat::poolify()` with fixed random seeds.

The final population-level dataset contains 50 units:

- 44 *D. suzukii*
- 4 *D. subpulchrella*
  - Dsub_CN
  - Dsub_KO
  - DsubPool
  - DsubPool2
- 2 *D. biarmipes*
  - DbiarRef
  - DbiarPool

## Repository organization

- `scripts/01_metadata/` — metadata construction and validation
- `scripts/02_poolfstat/` — sample- and population-level poolfstat datasets
- `scripts/03_pca/` — population structure analyses
- `scripts/04_fstats/` — genome-wide F3/F4 analyses
- `scripts/05_f4ratio/` — F4-ratio analyses
- `scripts/06_localized_introgression/` — localized allele-sharing analyses
- `scripts/07_figures_tables/` — manuscript figures and tables
- `docs/` — workflow documentation

Large genomic datasets and intermediate files are not stored in this
repository.
