# Rat skeletal-muscle contusion transcriptome time series

This repository contains analysis and figure-generation code associated with the Scientific Data Data Descriptor:

**A 0–72 h transcriptomic time series of rat skeletal muscle after contusion**

## Primary dataset

The 2026 primary cohort contains 77 bulk RNA-seq samples:

| Sampling stage | n |
|---|---:|
| Baseline | 7 |
| 0 h | 9 |
| 1 h | 8 |
| 2 h | 7 |
| 6 h | 6 |
| 10 h | 6 |
| 14 h | 7 |
| 18 h | 6 |
| 36 h | 7 |
| 60 h | 7 |
| 72 h | 7 |

Baseline represents uninjured skeletal muscle.

The 0-h group represents tissue collected immediately following experimental contusion and is not an uninjured control.

## Historical-reference cohort

A previously published 60-sample rat skeletal-muscle contusion cohort was independently reprocessed as a historical reference.

GEO accession: **GSE171243**

The two cohorts were analysed independently and were not subjected to cross-cohort batch correction.

## Repository structure

- `analysis/` — core statistical-analysis workflow
- `figures/` — main Figure 1–5 generation scripts
- `supplementary/` — Supplementary Figure S1–S8 analysis and plotting scripts
- `submission_qc/` — package and manuscript consistency utilities
- `docs/` — method definitions, software versions and selected session information

## Core workflow

1. Matrix-integrity checking
2. Count filtering and variance stabilization
3. Internal expression-structure analysis
4. Temporal multivariate analysis
5. Baseline-relative differential-expression analysis
6. Continuous-time spline likelihood-ratio testing
7. Tau/SPM temporal-specificity annotation
8. Historical-reference cohort reprocessing
9. Cross-cohort response-concordance analysis
10. Sample-retention and leave-one-out sensitivity analysis

## Expression representations

Raw gene counts were used for count filtering and statistical modelling.

Processed TPM values were used for descriptive expression summaries and temporal-specificity annotation.

DESeq2 variance-stabilized expression values were used for multivariate expression-space analyses.

## Reproducibility

The recommended execution order is provided in `docs/run_order.md`.

## Software

See `docs/software_versions.txt`.

## Input data

See `docs/input_files.md`.

Technical-validation scripts that require external sequencing-QC files use
the `SCIDATA_DATA_ROOT` environment variable. See
`docs/path_configuration.md` for details.

Large sequencing and processed-data files are not stored in this repository.

The 2026 primary RNA-seq repository accession will be added after deposition.

## Citation

Citation information for the associated Scientific Data Data Descriptor will be added when available.

## License

See `LICENSE`.
