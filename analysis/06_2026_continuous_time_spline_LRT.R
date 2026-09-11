#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  bitmapType = "cairo"
)

suppressPackageStartupMessages({

  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("R package 'DESeq2' is required.")
  }

  if (!requireNamespace("splines", quietly = TRUE)) {
    stop("R package 'splines' is required.")
  }

  library(DESeq2)
  library(splines)
})

cat("\n")
cat("============================================================\n")
cat("06 2026 CONTINUOUS POST-INJURY TEMPORAL MODELING\n")
cat("DESeq2 natural-spline likelihood-ratio test\n")
cat("Scientific Data rat skeletal muscle transcriptome project\n")
cat("============================================================\n\n")


# ============================================================
# 0. Fixed parameters
# ============================================================

EXPECTED_ALL_SAMPLES <- 77
EXPECTED_INJURY_SAMPLES <- 70
EXPECTED_INPUT_GENES <- 18364

MIN_COUNT <- 10
MIN_SAMPLES <- 6

SPLINE_DF <- 4
FDR_THRESHOLD <- 0.05

injury_group_levels <- c(
  "R0h",
  "R1h",
  "R2h",
  "R6h",
  "R10h",
  "R14h",
  "R18h",
  "R36h",
  "R60h",
  "R72h"
)

group_time_map <- c(
  R0h = 0,
  R1h = 1,
  R2h = 2,
  R6h = 6,
  R10h = 10,
  R14h = 14,
  R18h = 18,
  R36h = 36,
  R60h = 60,
  R72h = 72
)

group_label_map <- c(
  R0h = "0 h",
  R1h = "1 h",
  R2h = "2 h",
  R6h = "6 h",
  R10h = "10 h",
  R14h = "14 h",
  R18h = "18 h",
  R36h = "36 h",
  R60h = "60 h",
  R72h = "72 h"
)


# ============================================================
# 1. Input/output files
# ============================================================

count_file <- paste0(
  "02_count_filtering_and_VST_result/",
  "2026_primary/",
  "02_2026_filtered_counts.tsv"
)

meta_file <- paste0(
  "02_count_filtering_and_VST_result/",
  "2026_primary/",
  "02_2026_metadata_aligned.tsv"
)

outdir <- "06_2026_continuous_time_spline_LRT_result"

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 2. Check input files
# ============================================================

required_files <- c(
  count_file,
  meta_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  cat("ERROR: Missing input files:\n")

  for (f in missing_files) {
    cat("  -", f, "\n")
  }

  quit(status = 1)
}

cat("Input files found: PASS\n\n")


# ============================================================
# 3. Read filtered raw count matrix
# ============================================================

cat("Reading Step 02 filtered raw counts...\n")

count_df <- read.delim(
  count_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

colnames(count_df)[1] <- "GeneID"

if (anyDuplicated(count_df$GeneID)) {
  stop("Duplicated GeneIDs detected.")
}

if (any(
  is.na(count_df$GeneID) |
  trimws(count_df$GeneID) == ""
)) {
  stop("Missing GeneIDs detected.")
}

count_mat <- as.matrix(
  count_df[, -1, drop = FALSE]
)

storage.mode(count_mat) <- "numeric"

rownames(count_mat) <- count_df$GeneID

if (anyNA(count_mat)) {
  stop("NA values detected in count matrix.")
}

if (any(!is.finite(count_mat))) {
  stop("Non-finite values detected in count matrix.")
}

if (any(count_mat < 0)) {
  stop("Negative counts detected.")
}

if (any(
  abs(count_mat - round(count_mat)) > 1e-8
)) {
  stop("Noninteger counts detected.")
}

count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"


# ============================================================
# 4. Read metadata
# ============================================================

cat("Reading metadata...\n")

meta <- read.delim(
  meta_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

required_meta_cols <- c(
  "sampleID",
  "group",
  "time",
  "condition"
)

missing_meta_cols <- setdiff(
  required_meta_cols,
  colnames(meta)
)

if (length(missing_meta_cols) > 0) {

  stop(
    "Missing metadata columns: ",
    paste(
      missing_meta_cols,
      collapse = ", "
    )
  )
}

if (anyDuplicated(meta$sampleID)) {
  stop("Duplicated sampleIDs detected.")
}


# ============================================================
# 5. Align metadata to counts
# ============================================================

all_samples <- colnames(count_mat)

if (!setequal(
  all_samples,
  meta$sampleID
)) {

  stop(
    "Count matrix and metadata sample sets differ."
  )
}

meta <- meta[
  match(
    all_samples,
    meta$sampleID
  ),
  ,
  drop = FALSE
]

if (!identical(
  all_samples,
  meta$sampleID
)) {

  stop(
    "Metadata alignment failed."
  )
}

rownames(meta) <- meta$sampleID


# ============================================================
# 6. Validate Step 02 dimensions
# ============================================================

cat("\nStep 02 input dimensions:\n")

cat(
  "Genes       :",
  nrow(count_mat),
  "\n"
)

cat(
  "All samples :",
  ncol(count_mat),
  "\n\n"
)

if (nrow(count_mat) != EXPECTED_INPUT_GENES) {

  stop(
    "Expected ",
    EXPECTED_INPUT_GENES,
    " input genes, observed ",
    nrow(count_mat)
  )
}

if (ncol(count_mat) != EXPECTED_ALL_SAMPLES) {

  stop(
    "Expected ",
    EXPECTED_ALL_SAMPLES,
    " samples, observed ",
    ncol(count_mat)
  )
}

cat("Input dimension validation: PASS\n\n")


# ============================================================
# 7. Select injury samples only
# ============================================================

cat("Selecting injury samples only...\n")

injury_meta <- meta[
  meta$condition == "Injury",
  ,
  drop = FALSE
]

if (nrow(injury_meta) != EXPECTED_INJURY_SAMPLES) {

  stop(
    "Expected ",
    EXPECTED_INJURY_SAMPLES,
    " injury samples, observed ",
    nrow(injury_meta)
  )
}

unexpected_groups <- setdiff(
  unique(injury_meta$group),
  injury_group_levels
)

if (length(unexpected_groups) > 0) {

  stop(
    "Unexpected injury groups: ",
    paste(
      unexpected_groups,
      collapse = ", "
    )
  )
}

missing_groups <- setdiff(
  injury_group_levels,
  unique(injury_meta$group)
)

if (length(missing_groups) > 0) {

  stop(
    "Missing expected injury groups: ",
    paste(
      missing_groups,
      collapse = ", "
    )
  )
}

injury_meta$group <- factor(
  injury_meta$group,
  levels = injury_group_levels
)

# Explicit numeric post-injury time
injury_meta$time_h <- as.numeric(
  group_time_map[
    as.character(injury_meta$group)
  ]
)

if (anyNA(injury_meta$time_h)) {
  stop("Failed to assign numeric injury times.")
}

injury_samples <- injury_meta$sampleID

injury_counts <- count_mat[
  ,
  injury_samples,
  drop = FALSE
]

if (!identical(
  colnames(injury_counts),
  injury_meta$sampleID
)) {

  stop(
    "Injury count matrix and metadata order differ."
  )
}

cat(
  "Injury samples selected :",
  ncol(injury_counts),
  "\n"
)

cat(
  "Baseline samples excluded:",
  sum(meta$condition == "Control"),
  "\n\n"
)


# ============================================================
# 8. Injury-stage sample summary
# ============================================================

sample_summary <- data.frame(
  Group = injury_group_levels,
  Stage = unname(
    group_label_map[
      injury_group_levels
    ]
  ),
  Time_h = unname(
    group_time_map[
      injury_group_levels
    ]
  ),
  N = as.integer(
    table(
      factor(
        injury_meta$group,
        levels = injury_group_levels
      )
    )
  ),
  stringsAsFactors = FALSE
)

write.table(
  sample_summary,
  file = file.path(
    outdir,
    "06_2026_injury_sample_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Post-injury sample distribution:\n")
print(
  sample_summary,
  row.names = FALSE
)
cat("\n")


# ============================================================
# 9. Analysis-specific injury-sample expression filter
# ============================================================

cat(
  "Applying injury-specific expression eligibility filter...\n"
)

samples_ge_threshold <- rowSums(
  injury_counts >= MIN_COUNT
)

keep_injury <- samples_ge_threshold >= MIN_SAMPLES

genes_before <- nrow(injury_counts)
genes_after <- sum(keep_injury)
genes_removed <- genes_before - genes_after

all_zero_injury <- sum(
  rowSums(injury_counts) == 0
)

filter_metrics <- data.frame(
  GeneID = rownames(injury_counts),
  Injury_samples_with_count_GE_10 =
    samples_ge_threshold,
  Total_injury_raw_count =
    rowSums(injury_counts),
  Keep_for_continuous_time_model =
    keep_injury,
  stringsAsFactors = FALSE
)

write.table(
  filter_metrics,
  file = file.path(
    outdir,
    "06_2026_injury_time_model_filtering_metrics.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

eligible_counts <- injury_counts[
  keep_injury,
  ,
  drop = FALSE
]

write.table(
  data.frame(
    GeneID = rownames(eligible_counts),
    stringsAsFactors = FALSE
  ),
  file = file.path(
    outdir,
    "06_2026_time_model_eligible_gene_list.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("\nInjury-specific gene filtering:\n")

cat(
  "Genes entering Step 06      :",
  genes_before,
  "\n"
)

cat(
  "Genes retained for modeling :",
  genes_after,
  "\n"
)

cat(
  "Genes removed               :",
  genes_removed,
  "\n"
)

cat(
  "All-zero injury genes       :",
  all_zero_injury,
  "\n"
)

cat(
  "Criterion                   : count >=",
  MIN_COUNT,
  "in >=",
  MIN_SAMPLES,
  "injury samples\n\n"
)

if (genes_after < 1000) {
  stop("Unexpectedly few genes retained for time modeling.")
}


# ============================================================
# 10. Define spline model
# ============================================================

full_formula <- as.formula(
  paste0(
    "~ splines::ns(time_h, df = ",
    SPLINE_DF,
    ")"
  )
)

reduced_formula <- ~ 1

cat("Temporal model:\n")
cat(
  "Full    : ",
  deparse(full_formula),
  "\n",
  sep = ""
)

cat(
  "Reduced : ",
  deparse(reduced_formula),
  "\n",
  sep = ""
)

cat(
  "Spline df: ",
  SPLINE_DF,
  "\n\n",
  sep = ""
)


# ============================================================
# 11. Explicitly inspect spline basis and knots
# ============================================================

spline_basis <- splines::ns(
  injury_meta$time_h,
  df = SPLINE_DF
)

spline_knots <- attr(
  spline_basis,
  "knots"
)

boundary_knots <- attr(
  spline_basis,
  "Boundary.knots"
)

spline_parameters <- data.frame(
  Parameter = c(
    "Spline_df",
    paste0(
      "Internal_knot_",
      seq_along(spline_knots)
    ),
    "Boundary_knot_min",
    "Boundary_knot_max"
  ),
  Value = c(
    SPLINE_DF,
    as.numeric(spline_knots),
    boundary_knots[1],
    boundary_knots[2]
  ),
  stringsAsFactors = FALSE
)

write.table(
  spline_parameters,
  file = file.path(
    outdir,
    "06_2026_natural_spline_parameters.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Natural spline parameters:\n")
print(
  spline_parameters,
  row.names = FALSE
)
cat("\n")


# ============================================================
# 12. Validate model matrices
# ============================================================

full_model_matrix <- model.matrix(
  full_formula,
  data = injury_meta
)

reduced_model_matrix <- model.matrix(
  reduced_formula,
  data = injury_meta
)

full_rank <- qr(
  full_model_matrix
)$rank

reduced_rank <- qr(
  reduced_model_matrix
)$rank

cat("Model-matrix validation:\n")

cat(
  "Full matrix dimensions    :",
  nrow(full_model_matrix),
  "x",
  ncol(full_model_matrix),
  "\n"
)

cat(
  "Full matrix rank          :",
  full_rank,
  "\n"
)

cat(
  "Reduced matrix dimensions :",
  nrow(reduced_model_matrix),
  "x",
  ncol(reduced_model_matrix),
  "\n"
)

cat(
  "Reduced matrix rank       :",
  reduced_rank,
  "\n\n"
)

if (full_rank != ncol(full_model_matrix)) {
  stop("Full model matrix is not full rank.")
}

if (reduced_rank != ncol(reduced_model_matrix)) {
  stop("Reduced model matrix is not full rank.")
}

cat("Model matrix rank validation: PASS\n\n")


# ============================================================
# 13. Save model matrix
# ============================================================

full_model_output <- data.frame(
  sampleID = injury_meta$sampleID,
  time_h = injury_meta$time_h,
  full_model_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.table(
  full_model_output,
  file = file.path(
    outdir,
    "06_2026_full_spline_model_matrix.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 14. Create DESeq2 object
# ============================================================

cat("Creating DESeq2 time-course dataset...\n")

dds_time <- DESeqDataSetFromMatrix(
  countData = eligible_counts,
  colData = injury_meta,
  design = full_formula
)


# ============================================================
# 15. Run likelihood-ratio test
# ============================================================

cat(
  "Running DESeq2 likelihood-ratio test...\n"
)

cat(
  "Automatic Cook's count replacement: DISABLED\n\n"
)

dds_time <- DESeq(
  dds_time,
  test = "LRT",
  reduced = reduced_formula,
  minReplicatesForReplace = Inf,
  quiet = FALSE
)

cat("\nDESeq2 spline LRT fitting completed.\n\n")


# ============================================================
# 16. Extract LRT results
# ============================================================

# Independent filtering is disabled because an explicit
# count-based analysis eligibility filter was already applied.
#
# Cook's-distance diagnostics remain active. Genes considered
# unsuitable by DESeq2 may therefore have NA test results.

lrt_res <- results(
  dds_time,
  alpha = FDR_THRESHOLD,
  independentFiltering = FALSE,
  cooksCutoff = TRUE
)

lrt_df_native <- as.data.frame(
  lrt_res
)

lrt_df_native$GeneID <- rownames(
  lrt_df_native
)

# IMPORTANT:
# For an LRT, the statistical test evaluates the full spline
# model against the intercept-only model. The reported
# log2FoldChange does NOT represent the omnibus temporal effect.
# Therefore only baseMean/stat/pvalue/padj are retained in the
# publication-facing temporal statistics table.

lrt_df <- lrt_df_native[
  ,
  c(
    "GeneID",
    "baseMean",
    "stat",
    "pvalue",
    "padj"
  ),
  drop = FALSE
]

lrt_df$Dynamic_FDR005 <- (
  !is.na(lrt_df$padj) &
  lrt_df$padj < FDR_THRESHOLD
)


# ============================================================
# 17. Cook's-distance diagnostics
# ============================================================

if ("cooks" %in% assayNames(dds_time)) {

  cooks_mat <- assays(
    dds_time
  )[["cooks"]]

  max_cooks <- apply(
    cooks_mat,
    1,
    function(x) {

      if (all(is.na(x))) {
        return(NA_real_)
      }

      max(
        x,
        na.rm = TRUE
      )
    }
  )

  lrt_df$Maximum_Cooks_distance <-
    max_cooks[
      lrt_df$GeneID
    ]

} else {

  lrt_df$Maximum_Cooks_distance <-
    NA_real_
}


# ============================================================
# 18. Sort and save complete results
# ============================================================

ord <- order(
  is.na(lrt_df$padj),
  lrt_df$padj,
  lrt_df$pvalue
)

lrt_df <- lrt_df[
  ord,
  ,
  drop = FALSE
]

write.table(
  lrt_df,
  file = file.path(
    outdir,
    "06_2026_spline_LRT_full_statistics.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 19. Dynamic gene list
# ============================================================

dynamic_df <- lrt_df[
  lrt_df$Dynamic_FDR005,
  ,
  drop = FALSE
]

write.table(
  dynamic_df,
  file = file.path(
    outdir,
    "06_2026_dynamic_genes_FDR005.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  data.frame(
    GeneID = dynamic_df$GeneID,
    stringsAsFactors = FALSE
  ),
  file = file.path(
    outdir,
    "06_2026_dynamic_gene_IDs_FDR005.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 20. Result diagnostics
# ============================================================

n_tested <- nrow(lrt_df)

n_pvalue_na <- sum(
  is.na(lrt_df$pvalue)
)

n_padj_na <- sum(
  is.na(lrt_df$padj)
)

n_dynamic <- nrow(dynamic_df)

n_fdr001 <- sum(
  !is.na(lrt_df$padj) &
  lrt_df$padj < 0.01
)

n_fdr0001 <- sum(
  !is.na(lrt_df$padj) &
  lrt_df$padj < 0.001
)

dynamic_percent <- 100 *
  n_dynamic /
  n_tested


# ============================================================
# 21. Timepoint-level descriptive expression summaries
# using DESeq2 VST only for trajectory support
# ============================================================

cat(
  "Generating VST stage-level expression summaries ",
  "for downstream trajectory visualization...\n"
)

vsd_time <- varianceStabilizingTransformation(
  dds_time,
  blind = FALSE
)

vst_time_mat <- assay(
  vsd_time
)

stage_mean_vst <- sapply(
  injury_group_levels,
  function(g) {

    samples <- injury_meta$sampleID[
      injury_meta$group == g
    ]

    rowMeans(
      vst_time_mat[
        ,
        samples,
        drop = FALSE
      ]
    )
  }
)

colnames(stage_mean_vst) <-
  injury_group_levels

stage_median_vst <- sapply(
  injury_group_levels,
  function(g) {

    samples <- injury_meta$sampleID[
      injury_meta$group == g
    ]

    apply(
      vst_time_mat[
        ,
        samples,
        drop = FALSE
      ],
      1,
      median
    )
  }
)

colnames(stage_median_vst) <-
  injury_group_levels


# ============================================================
# 22. Calculate descriptive temporal amplitude
# ============================================================

vst_range <- apply(
  stage_mean_vst,
  1,
  function(x) {
    max(x) - min(x)
  }
)

vst_sd <- apply(
  stage_mean_vst,
  1,
  sd
)

peak_index <- max.col(
  stage_mean_vst,
  ties.method = "first"
)

trough_index <- max.col(
  -stage_mean_vst,
  ties.method = "first"
)

descriptive_metrics <- data.frame(
  GeneID = rownames(stage_mean_vst),
  Mean_VST_temporal_range =
    vst_range,
  SD_of_stage_mean_VST =
    vst_sd,
  Highest_mean_VST_stage =
    injury_group_levels[
      peak_index
    ],
  Lowest_mean_VST_stage =
    injury_group_levels[
      trough_index
    ],
  stringsAsFactors = FALSE
)

write.table(
  descriptive_metrics,
  file = file.path(
    outdir,
    "06_2026_time_model_descriptive_VST_metrics.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 23. Append descriptive metrics to dynamic gene table
# ============================================================

dynamic_annotated <- merge(
  dynamic_df,
  descriptive_metrics,
  by = "GeneID",
  all.x = TRUE,
  sort = FALSE
)

dynamic_annotated <- dynamic_annotated[
  match(
    dynamic_df$GeneID,
    dynamic_annotated$GeneID
  ),
  ,
  drop = FALSE
]

write.table(
  dynamic_annotated,
  file = file.path(
    outdir,
    "06_2026_dynamic_genes_FDR005_with_VST_metrics.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 24. Save stage-mean and stage-median VST matrices
# ============================================================

write.table(
  data.frame(
    GeneID = rownames(stage_mean_vst),
    stage_mean_vst,
    check.names = FALSE
  ),
  file = file.path(
    outdir,
    "06_2026_stage_mean_VST.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  data.frame(
    GeneID = rownames(stage_median_vst),
    stage_median_vst,
    check.names = FALSE
  ),
  file = file.path(
    outdir,
    "06_2026_stage_median_VST.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 25. Summary table
# ============================================================

summary_table <- data.frame(
  Metric = c(
    "All_2026_samples",
    "Baseline_samples_excluded",
    "Injury_samples_modeled",
    "Step02_input_genes",
    "Genes_retained_for_time_model",
    "Genes_removed_by_injury_specific_filter",
    "Minimum_count_threshold",
    "Minimum_sample_threshold",
    "Spline_df",
    "Full_model_columns",
    "Full_model_rank",
    "Genes_tested",
    "Genes_with_Pvalue_NA",
    "Genes_with_PADJ_NA",
    "Dynamic_genes_FDR_LT_0.05",
    "Dynamic_genes_FDR_LT_0.01",
    "Dynamic_genes_FDR_LT_0.001",
    "Dynamic_gene_percent"
  ),
  Value = c(
    EXPECTED_ALL_SAMPLES,
    EXPECTED_ALL_SAMPLES -
      EXPECTED_INJURY_SAMPLES,
    EXPECTED_INJURY_SAMPLES,
    genes_before,
    genes_after,
    genes_removed,
    MIN_COUNT,
    MIN_SAMPLES,
    SPLINE_DF,
    ncol(full_model_matrix),
    full_rank,
    n_tested,
    n_pvalue_na,
    n_padj_na,
    n_dynamic,
    n_fdr001,
    n_fdr0001,
    dynamic_percent
  ),
  stringsAsFactors = FALSE
)

write.table(
  summary_table,
  file = file.path(
    outdir,
    "06_2026_continuous_time_model_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 26. Save DESeq2 objects
# ============================================================

saveRDS(
  dds_time,
  file = file.path(
    outdir,
    "06_2026_spline_LRT_DESeq2_object.rds"
  )
)

saveRDS(
  vsd_time,
  file = file.path(
    outdir,
    "06_2026_time_model_VST_object.rds"
  )
)


# ============================================================
# 27. Save analysis metadata
# ============================================================

analysis_meta <- as.data.frame(
  injury_meta
)

analysis_meta$group <- as.character(
  analysis_meta$group
)

write.table(
  analysis_meta,
  file = file.path(
    outdir,
    "06_2026_injury_metadata_used.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 28. Save session information
# ============================================================

sink(
  file.path(
    outdir,
    "06_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 29. Final console report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("06 CONTINUOUS TIME SPLINE-LRT COMPLETED\n")
cat("============================================================\n")

cat(
  "All 2026 samples                : ",
  EXPECTED_ALL_SAMPLES,
  "\n",
  sep = ""
)

cat(
  "Baseline samples excluded       : ",
  EXPECTED_ALL_SAMPLES -
    EXPECTED_INJURY_SAMPLES,
  "\n",
  sep = ""
)

cat(
  "Injury samples modeled          : ",
  EXPECTED_INJURY_SAMPLES,
  "\n",
  sep = ""
)

cat("------------------------------------------------------------\n")

cat(
  "Step 02 genes                   : ",
  genes_before,
  "\n",
  sep = ""
)

cat(
  "Genes eligible for time model   : ",
  genes_after,
  "\n",
  sep = ""
)

cat(
  "Genes removed for Step 06       : ",
  genes_removed,
  "\n",
  sep = ""
)

cat("------------------------------------------------------------\n")

cat(
  "Spline degrees of freedom       : ",
  SPLINE_DF,
  "\n",
  sep = ""
)

cat(
  "Full model matrix rank          : ",
  full_rank,
  "/",
  ncol(full_model_matrix),
  "\n",
  sep = ""
)

cat(
  "FDR threshold                   : ",
  FDR_THRESHOLD,
  "\n",
  sep = ""
)

cat("------------------------------------------------------------\n")

cat(
  "Genes tested                    : ",
  n_tested,
  "\n",
  sep = ""
)

cat(
  "Genes with p-value NA           : ",
  n_pvalue_na,
  "\n",
  sep = ""
)

cat(
  "Genes with adjusted P NA        : ",
  n_padj_na,
  "\n",
  sep = ""
)

cat(
  "Dynamic genes FDR < 0.05        : ",
  n_dynamic,
  "\n",
  sep = ""
)

cat(
  "Dynamic genes FDR < 0.01        : ",
  n_fdr001,
  "\n",
  sep = ""
)

cat(
  "Dynamic genes FDR < 0.001       : ",
  n_fdr0001,
  "\n",
  sep = ""
)

cat(
  "Dynamic genes (%)               : ",
  sprintf(
    "%.2f",
    dynamic_percent
  ),
  "\n",
  sep = ""
)

cat("------------------------------------------------------------\n")
cat("Baseline was NOT included in the continuous-time model.\n")
cat("No CPM or logCPM was used.\n")
cat("Raw counts were analyzed with DESeq2 LRT.\n")
cat("Automatic Cook's count replacement was disabled.\n")
cat("DESeq2 Cook's-distance diagnostics remained active.\n")
cat("Independent filtering was disabled after explicit count filtering.\n")
cat("The LRT tests the overall temporal contribution of the spline.\n")
cat("LRT log2FoldChange is NOT interpreted as an omnibus effect size.\n")
cat("No pathway enrichment analysis was performed.\n")
cat("============================================================\n\n")

cat("Output directory:\n")
cat(outdir, "\n\n")

cat("PASS\n")
cat(
  "Step 06 is complete. ",
  "The next step is Step 07: ",
  "2026 TMM.TPM group-median expression, ",
  "Tau, SPM and peak-stage annotations.\n"
)
