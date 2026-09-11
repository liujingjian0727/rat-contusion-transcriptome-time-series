# Analysis run order

## Core analysis

Run the principal analysis scripts in the following order:

1. `analysis/01_matrix_integrity_check.R`
2. `analysis/02_count_filtering_and_VST.R`
3. `analysis/03_2026_internal_expression_structure.R`
4. `analysis/04_2026_temporal_structure.R`
5. `analysis/05_2026_timepoint_vs_baseline_DESeq2.R`
6. `analysis/06_2026_continuous_time_spline_LRT.R`
7. `analysis/06B_define_high_confidence_dynamic_genes.R`
8. `analysis/07_2026_TPM_tau_SPM_peak_stage.R`
9. `analysis/08_2021_historical_reference_reanalysis.R`
10. `analysis/09_cross_cohort_temporal_response_concordance.R`
11. `analysis/10_2026_sample_exclusion_sensitivity.R`
12. `analysis/10B_single_sample_leave_one_out.R`

## Main figures

- `figures/11C_make_Figure1_resource_overview_v2.R`
- `figures/11D_make_Figure2_technical_validation.R`
- `figures/11B_make_final_Figures_3_4_5.R`
- `figures/11A_make_Figure4C_4D.R`

## Supplementary analyses

Scripts beginning with `12`, `13` and `14` implement technical-validation
analyses and generate Supplementary Figures S1-S8.

## Analysis conventions

- Raw gene counts are used for filtering and statistical modelling.
- Processed TPM values are used for descriptive expression summaries and
  temporal-specificity annotation.
- DESeq2 variance-stabilized expression values are used for multivariate analyses.
- Baseline represents uninjured skeletal muscle.
- 0 h represents tissue collected immediately after injury and is not Baseline.
- The 2026 primary cohort and 2021 historical-reference cohort are analysed
  independently.
- No cross-cohort batch correction is applied.
