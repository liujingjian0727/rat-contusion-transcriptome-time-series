#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1
)

cat("\n")
cat("============================================================\n")
cat("07 2026 TMM.TPM TEMPORAL SPECIFICITY ANNOTATION\n")
cat("Group median expression, Tau, SPM and peak-stage annotation\n")
cat("Scientific Data rat skeletal muscle transcriptome project\n")
cat("============================================================\n\n")


# ============================================================
# 0. Fixed parameters
# ============================================================

EXPECTED_GENES <- 18364
EXPECTED_SAMPLES <- 77

TAU_THRESHOLD <- 0.85
SPM_THRESHOLD <- 0.50

all_group_levels <- c(
  "Baseline",
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

group_labels <- c(
  Baseline = "Baseline",
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
# 1. Input/output
# ============================================================

tpm_file <- paste0(
  "02_count_filtering_and_VST_result/",
  "2026_primary/",
  "02_2026_filtered_TMM_TPM.tsv"
)

meta_file <- paste0(
  "02_count_filtering_and_VST_result/",
  "2026_primary/",
  "02_2026_metadata_aligned.tsv"
)

dynamic_file <- paste0(
  "06_2026_continuous_time_spline_LRT_result/",
  "06_2026_dynamic_gene_IDs_FDR005.tsv"
)

high_conf_file <- paste0(
  "06_2026_continuous_time_spline_LRT_result/",
  "06B_high_confidence_dynamic_gene_IDs.tsv"
)

outdir <- "07_2026_TPM_tau_SPM_peak_stage_result"

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 2. Check input files
# ============================================================

required_files <- c(
  tpm_file,
  meta_file,
  dynamic_file,
  high_conf_file
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

cat("All input files found: PASS\n\n")


# ============================================================
# 3. Read TPM matrix
# ============================================================

cat("Reading filtered TMM.TPM matrix...\n")

tpm_df <- read.delim(
  tpm_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

colnames(tpm_df)[1] <- "GeneID"

if (anyDuplicated(tpm_df$GeneID)) {
  stop("Duplicated GeneIDs detected in TPM matrix.")
}

if (any(
  is.na(tpm_df$GeneID) |
  trimws(tpm_df$GeneID) == ""
)) {
  stop("Missing GeneIDs detected.")
}

tpm_mat <- as.matrix(
  tpm_df[, -1, drop = FALSE]
)

storage.mode(tpm_mat) <- "numeric"

rownames(tpm_mat) <- tpm_df$GeneID

if (anyNA(tpm_mat)) {
  stop("NA values detected in TPM matrix.")
}

if (any(!is.finite(tpm_mat))) {
  stop("Non-finite values detected in TPM matrix.")
}

if (any(tpm_mat < 0)) {
  stop("Negative TPM values detected.")
}


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

required_meta <- c(
  "sampleID",
  "group",
  "time",
  "condition"
)

missing_meta <- setdiff(
  required_meta,
  colnames(meta)
)

if (length(missing_meta) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_meta, collapse = ", ")
  )
}

if (anyDuplicated(meta$sampleID)) {
  stop("Duplicated sampleIDs detected.")
}


# ============================================================
# 5. Explicit sample matching
# ============================================================

sample_ids <- colnames(tpm_mat)

if (!setequal(
  sample_ids,
  meta$sampleID
)) {
  stop("TPM and metadata sample sets differ.")
}

meta <- meta[
  match(
    sample_ids,
    meta$sampleID
  ),
  ,
  drop = FALSE
]

if (!identical(
  sample_ids,
  meta$sampleID
)) {
  stop("Metadata sample alignment failed.")
}

rownames(meta) <- meta$sampleID

meta$group <- factor(
  meta$group,
  levels = all_group_levels
)

if (anyNA(meta$group)) {
  stop("Unexpected group detected.")
}


# ============================================================
# 6. Dimension checks
# ============================================================

cat("\nInput dimensions:\n")

cat(
  "TMM.TPM  :",
  nrow(tpm_mat),
  "genes x",
  ncol(tpm_mat),
  "samples\n"
)

cat(
  "Metadata :",
  nrow(meta),
  "samples\n\n"
)

if (nrow(tpm_mat) != EXPECTED_GENES) {
  stop(
    "Expected ",
    EXPECTED_GENES,
    " genes, observed ",
    nrow(tpm_mat)
  )
}

if (ncol(tpm_mat) != EXPECTED_SAMPLES) {
  stop(
    "Expected ",
    EXPECTED_SAMPLES,
    " samples, observed ",
    ncol(tpm_mat)
  )
}

cat("Dimension validation: PASS\n\n")


# ============================================================
# 7. Read dynamic gene memberships
# ============================================================

dynamic_df <- read.delim(
  dynamic_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

high_conf_df <- read.delim(
  high_conf_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"GeneID" %in% colnames(dynamic_df)) {
  stop("GeneID column missing from dynamic gene file.")
}

if (!"GeneID" %in% colnames(high_conf_df)) {
  stop("GeneID column missing from high-confidence gene file.")
}

dynamic_ids <- unique(
  dynamic_df$GeneID
)

high_conf_ids <- unique(
  high_conf_df$GeneID
)

cat(
  "LRT dynamic gene IDs              :",
  length(dynamic_ids),
  "\n"
)

cat(
  "High-confidence dynamic gene IDs  :",
  length(high_conf_ids),
  "\n\n"
)


# ============================================================
# 8. Group sample summary
# ============================================================

group_summary <- data.frame(
  Group = all_group_levels,
  Stage = unname(
    group_labels[
      all_group_levels
    ]
  ),
  N = as.integer(
    table(
      factor(
        meta$group,
        levels = all_group_levels
      )
    )
  ),
  stringsAsFactors = FALSE
)

write.table(
  group_summary,
  file = file.path(
    outdir,
    "07_2026_group_sample_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Sample distribution:\n")
print(
  group_summary,
  row.names = FALSE
)
cat("\n")


# ============================================================
# 9. Calculate group-median TPM
# ============================================================

cat("Calculating group-median TMM.TPM...\n")

group_median_tpm <- sapply(
  all_group_levels,
  function(g) {

    samples <- meta$sampleID[
      meta$group == g
    ]

    apply(
      tpm_mat[
        ,
        samples,
        drop = FALSE
      ],
      1,
      median
    )
  }
)

rownames(group_median_tpm) <- rownames(
  tpm_mat
)

colnames(group_median_tpm) <- all_group_levels


# ============================================================
# 10. Group mean TPM also retained as descriptive resource
# ============================================================

group_mean_tpm <- sapply(
  all_group_levels,
  function(g) {

    samples <- meta$sampleID[
      meta$group == g
    ]

    rowMeans(
      tpm_mat[
        ,
        samples,
        drop = FALSE
      ]
    )
  }
)

rownames(group_mean_tpm) <- rownames(
  tpm_mat
)

colnames(group_mean_tpm) <- all_group_levels


# ============================================================
# 11. Save median and mean TPM matrices
# ============================================================

write.table(
  data.frame(
    GeneID = rownames(group_median_tpm),
    group_median_tpm,
    check.names = FALSE
  ),
  file = file.path(
    outdir,
    "07_2026_group_median_TMM_TPM_all11stages.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  data.frame(
    GeneID = rownames(group_mean_tpm),
    group_mean_tpm,
    check.names = FALSE
  ),
  file = file.path(
    outdir,
    "07_2026_group_mean_TMM_TPM_all11stages.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 12. log2(TPM + 1) median matrix for visualization only
# ============================================================

group_median_log2tpm <- log2(
  group_median_tpm + 1
)

write.table(
  data.frame(
    GeneID = rownames(group_median_log2tpm),
    group_median_log2tpm,
    check.names = FALSE
  ),
  file = file.path(
    outdir,
    "07_2026_group_median_log2TPMplus1_all11stages.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 13. Functions for Tau and SPM
# ============================================================

calculate_tau <- function(x) {

  if (anyNA(x)) {
    return(NA_real_)
  }

  xmax <- max(x)

  if (!is.finite(xmax) || xmax <= 0) {
    return(NA_real_)
  }

  n <- length(x)

  if (n < 2) {
    return(NA_real_)
  }

  sum(
    1 - x / xmax
  ) / (n - 1)
}


calculate_spm <- function(x) {

  if (anyNA(x)) {
    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }

  denominator <- sum(
    x^2
  )

  if (!is.finite(denominator) ||
      denominator <= 0) {

    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }

  (x^2) / denominator
}


# ============================================================
# 14. Generic temporal-specificity function
# ============================================================

calculate_specificity_set <- function(
    expression_matrix,
    stage_names,
    prefix) {

  cat(
    "Calculating ",
    prefix,
    " Tau/SPM annotations...\n",
    sep = ""
  )

  x <- expression_matrix[
    ,
    stage_names,
    drop = FALSE
  ]

  # Tau
  tau <- apply(
    x,
    1,
    calculate_tau
  )

  # SPM
  spm <- t(
    apply(
      x,
      1,
      calculate_spm
    )
  )

  rownames(spm) <- rownames(x)

  colnames(spm) <- paste0(
    stage_names,
    "_SPM"
  )

  # Peak median TPM
  peak_index <- max.col(
    x,
    ties.method = "first"
  )

  peak_stage <- stage_names[
    peak_index
  ]

  peak_tpm <- x[
    cbind(
      seq_len(nrow(x)),
      peak_index
    )
  ]

  # Peak SPM
  peak_spm <- spm[
    cbind(
      seq_len(nrow(spm)),
      peak_index
    )
  ]

  # Number of tied maximum stages
  number_peak_ties <- apply(
    x,
    1,
    function(v) {
      sum(v == max(v))
    }
  )

  all_zero <- rowSums(x) == 0

  tau[
    all_zero
  ] <- NA_real_

  peak_stage[
    all_zero
  ] <- NA_character_

  peak_tpm[
    all_zero
  ] <- NA_real_

  peak_spm[
    all_zero
  ] <- NA_real_

  number_peak_ties[
    all_zero
  ] <- NA_integer_

  # Specificity flags
  tau_flag <- (
    !is.na(tau) &
    tau >= TAU_THRESHOLD
  )

  peak_spm_flag <- (
    !is.na(peak_spm) &
    peak_spm >= SPM_THRESHOLD
  )

  specific_flag <- (
    tau_flag &
    peak_spm_flag
  )

  # Summary annotation
  annotation <- data.frame(
    GeneID = rownames(x),
    Tau = tau,
    Peak_stage = peak_stage,
    Peak_stage_label = unname(
      group_labels[
        peak_stage
      ]
    ),
    Peak_median_TPM = peak_tpm,
    Peak_SPM = peak_spm,
    Number_of_peak_ties = number_peak_ties,
    Tau_GE_0.85 = tau_flag,
    Peak_SPM_GE_0.50 = peak_spm_flag,
    Tau085_and_SPM050 = specific_flag,
    stringsAsFactors = FALSE
  )

  # Add all SPM columns
  annotation <- cbind(
    annotation,
    as.data.frame(
      spm,
      check.names = FALSE
    )
  )

  # Save full annotation
  write.table(
    annotation,
    file = file.path(
      outdir,
      paste0(
        "07_2026_",
        prefix,
        "_Tau_SPM_peak_annotation.tsv"
      )
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  # Save Tau table
  write.table(
    annotation[
      ,
      c(
        "GeneID",
        "Tau",
        "Peak_stage",
        "Peak_stage_label",
        "Peak_median_TPM",
        "Peak_SPM"
      )
    ],
    file = file.path(
      outdir,
      paste0(
        "07_2026_",
        prefix,
        "_Tau_peak_summary.tsv"
      )
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  # Save SPM only
  write.table(
    data.frame(
      GeneID = rownames(spm),
      spm,
      check.names = FALSE
    ),
    file = file.path(
      outdir,
      paste0(
        "07_2026_",
        prefix,
        "_SPM_matrix.tsv"
      )
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  # Save specific genes
  specific_df <- annotation[
    annotation$Tau085_and_SPM050,
    ,
    drop = FALSE
  ]

  write.table(
    specific_df,
    file = file.path(
      outdir,
      paste0(
        "07_2026_",
        prefix,
        "_time_specific_genes_tau085_spm050.tsv"
      )
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  return(
    list(
      annotation = annotation,
      specific = specific_df
    )
  )
}


# ============================================================
# 15. All 11 sampling stages
# ============================================================

result_all11 <- calculate_specificity_set(
  expression_matrix = group_median_tpm,
  stage_names = all_group_levels,
  prefix = "all11stages"
)


# ============================================================
# 16. Injury-only 10 sampling stages
# ============================================================

result_injury10 <- calculate_specificity_set(
  expression_matrix = group_median_tpm,
  stage_names = injury_group_levels,
  prefix = "injury10stages"
)


# ============================================================
# 17. Integrate LRT/high-confidence membership
# ============================================================

cat(
  "Integrating temporal-model membership annotations...\n"
)

master_annotation <- result_all11$annotation[
  ,
  c(
    "GeneID",
    "Tau",
    "Peak_stage",
    "Peak_stage_label",
    "Peak_median_TPM",
    "Peak_SPM",
    "Number_of_peak_ties",
    "Tau_GE_0.85",
    "Peak_SPM_GE_0.50",
    "Tau085_and_SPM050"
  ),
  drop = FALSE
]

colnames(master_annotation)[
  colnames(master_annotation) == "Tau"
] <- "Tau_all11"

colnames(master_annotation)[
  colnames(master_annotation) == "Peak_stage"
] <- "Peak_stage_all11"

colnames(master_annotation)[
  colnames(master_annotation) == "Peak_stage_label"
] <- "Peak_stage_label_all11"

colnames(master_annotation)[
  colnames(master_annotation) == "Peak_median_TPM"
] <- "Peak_median_TPM_all11"

colnames(master_annotation)[
  colnames(master_annotation) == "Peak_SPM"
] <- "Peak_SPM_all11"

colnames(master_annotation)[
  colnames(master_annotation) == "Number_of_peak_ties"
] <- "Number_of_peak_ties_all11"

colnames(master_annotation)[
  colnames(master_annotation) == "Tau_GE_0.85"
] <- "Tau_GE_0.85_all11"

colnames(master_annotation)[
  colnames(master_annotation) == "Peak_SPM_GE_0.50"
] <- "Peak_SPM_GE_0.50_all11"

colnames(master_annotation)[
  colnames(master_annotation) == "Tau085_and_SPM050"
] <- "Tau085_and_SPM050_all11"


injury_annotation_small <- result_injury10$annotation[
  ,
  c(
    "GeneID",
    "Tau",
    "Peak_stage",
    "Peak_stage_label",
    "Peak_median_TPM",
    "Peak_SPM",
    "Number_of_peak_ties",
    "Tau_GE_0.85",
    "Peak_SPM_GE_0.50",
    "Tau085_and_SPM050"
  ),
  drop = FALSE
]

colnames(
  injury_annotation_small
)[-1] <- paste0(
  colnames(
    injury_annotation_small
  )[-1],
  "_injury10"
)


master_annotation <- merge(
  master_annotation,
  injury_annotation_small,
  by = "GeneID",
  all.x = TRUE,
  sort = FALSE
)

master_annotation <- master_annotation[
  match(
    rownames(group_median_tpm),
    master_annotation$GeneID
  ),
  ,
  drop = FALSE
]

master_annotation$LRT_dynamic_FDR005 <-
  master_annotation$GeneID %in%
  dynamic_ids

master_annotation$High_confidence_dynamic <-
  master_annotation$GeneID %in%
  high_conf_ids


# ============================================================
# 18. Add median TPM for every stage to master table
# ============================================================

median_stage_table <- data.frame(
  GeneID = rownames(group_median_tpm),
  group_median_tpm,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

master_annotation <- merge(
  master_annotation,
  median_stage_table,
  by = "GeneID",
  all.x = TRUE,
  sort = FALSE
)

master_annotation <- master_annotation[
  match(
    rownames(group_median_tpm),
    master_annotation$GeneID
  ),
  ,
  drop = FALSE
]

write.table(
  master_annotation,
  file = file.path(
    outdir,
    "07_2026_master_temporal_specificity_annotation.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 19. High-confidence dynamic subset with Tau/SPM
# ============================================================

hc_annotation <- master_annotation[
  master_annotation$High_confidence_dynamic,
  ,
  drop = FALSE
]

write.table(
  hc_annotation,
  file = file.path(
    outdir,
    "07_2026_high_confidence_dynamic_Tau_SPM_annotation.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 20. High-confidence + injury-stage-specific subset
# ============================================================

hc_injury_specific <- master_annotation[
  master_annotation$High_confidence_dynamic &
  !is.na(
    master_annotation$Tau085_and_SPM050_injury10
  ) &
  master_annotation$Tau085_and_SPM050_injury10,
  ,
  drop = FALSE
]

write.table(
  hc_injury_specific,
  file = file.path(
    outdir,
    paste0(
      "07_2026_high_confidence_dynamic_",
      "injury_time_specific_tau085_spm050.tsv"
    )
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 21. Stage-specific gene counts
# ============================================================

count_by_peak_stage <- function(
    df,
    peak_col,
    flag_col,
    stage_levels,
    output_name) {

  x <- df[
    !is.na(df[[flag_col]]) &
    df[[flag_col]],
    ,
    drop = FALSE
  ]

  tab <- table(
    factor(
      x[[peak_col]],
      levels = stage_levels
    )
  )

  out <- data.frame(
    Peak_stage = stage_levels,
    Peak_stage_label = unname(
      group_labels[
        stage_levels
      ]
    ),
    Gene_number = as.integer(tab),
    stringsAsFactors = FALSE
  )

  write.table(
    out,
    file = file.path(
      outdir,
      output_name
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  return(out)
}


all11_stage_counts <- count_by_peak_stage(
  master_annotation,
  "Peak_stage_all11",
  "Tau085_and_SPM050_all11",
  all_group_levels,
  "07_2026_all11_time_specific_gene_numbers.tsv"
)

injury10_stage_counts <- count_by_peak_stage(
  master_annotation,
  "Peak_stage_injury10",
  "Tau085_and_SPM050_injury10",
  injury_group_levels,
  "07_2026_injury10_time_specific_gene_numbers.tsv"
)


# ============================================================
# 22. High-confidence injury-specific gene counts
# ============================================================

hc_peak_tab <- table(
  factor(
    hc_injury_specific$Peak_stage_injury10,
    levels = injury_group_levels
  )
)

hc_stage_counts <- data.frame(
  Peak_stage = injury_group_levels,
  Peak_stage_label = unname(
    group_labels[
      injury_group_levels
    ]
  ),
  Gene_number = as.integer(
    hc_peak_tab
  ),
  stringsAsFactors = FALSE
)

write.table(
  hc_stage_counts,
  file = file.path(
    outdir,
    paste0(
      "07_2026_high_confidence_dynamic_",
      "injury_time_specific_gene_numbers.tsv"
    )
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 23. Summary statistics
# ============================================================

n_allzero_all11 <- sum(
  is.na(
    result_all11$annotation$Tau
  )
)

n_allzero_injury10 <- sum(
  is.na(
    result_injury10$annotation$Tau
  )
)

n_tau_all11 <- sum(
  result_all11$annotation$Tau_GE_0.85,
  na.rm = TRUE
)

n_spm_all11 <- sum(
  result_all11$annotation$Peak_SPM_GE_0.50,
  na.rm = TRUE
)

n_specific_all11 <- nrow(
  result_all11$specific
)

n_tau_injury10 <- sum(
  result_injury10$annotation$Tau_GE_0.85,
  na.rm = TRUE
)

n_spm_injury10 <- sum(
  result_injury10$annotation$Peak_SPM_GE_0.50,
  na.rm = TRUE
)

n_specific_injury10 <- nrow(
  result_injury10$specific
)

summary_table <- data.frame(
  Metric = c(
    "Filtered_genes",
    "Samples",
    "All11_stages",
    "Injury_only_stages",
    "LRT_dynamic_FDR005",
    "High_confidence_dynamic",
    "Genes_with_undefined_Tau_all11",
    "Genes_with_undefined_Tau_injury10",
    "Tau_GE_0.85_all11",
    "Peak_SPM_GE_0.50_all11",
    "Tau085_and_SPM050_all11",
    "Tau_GE_0.85_injury10",
    "Peak_SPM_GE_0.50_injury10",
    "Tau085_and_SPM050_injury10",
    "High_confidence_and_injury_specific"
  ),
  Value = c(
    nrow(tpm_mat),
    ncol(tpm_mat),
    length(all_group_levels),
    length(injury_group_levels),
    length(dynamic_ids),
    length(high_conf_ids),
    n_allzero_all11,
    n_allzero_injury10,
    n_tau_all11,
    n_spm_all11,
    n_specific_all11,
    n_tau_injury10,
    n_spm_injury10,
    n_specific_injury10,
    nrow(hc_injury_specific)
  ),
  stringsAsFactors = FALSE
)

write.table(
  summary_table,
  file = file.path(
    outdir,
    "07_2026_Tau_SPM_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 24. Save method-definition file
# ============================================================

method_lines <- c(
  "2026 temporal specificity annotation",
  "",
  paste0(
    "Input genes: ",
    EXPECTED_GENES
  ),
  paste0(
    "Input samples: ",
    EXPECTED_SAMPLES
  ),
  "",
  "Expression statistic:",
  "Median TMM.TPM across biological replicates within each sampling stage.",
  "",
  "Tau definition:",
  "Tau = sum(1 - x_i / max(x)) / (n - 1)",
  "where x_i is the median TMM.TPM at sampling stage i.",
  "",
  "SPM definition:",
  "SPM_i = x_i^2 / sum_j(x_j^2)",
  "The SPM values across stages sum to 1 for genes with non-zero expression.",
  "",
  paste0(
    "Tau threshold: ",
    TAU_THRESHOLD
  ),
  paste0(
    "Peak-stage SPM threshold: ",
    SPM_THRESHOLD
  ),
  "",
  "Tau/SPM are descriptive temporal-specificity annotations.",
  "They are not differential-expression significance tests.",
  "",
  "Two annotation systems were generated:",
  "1. all11stages: Baseline plus 10 post-injury sampling stages.",
  "2. injury10stages: 10 post-injury sampling stages only.",
  "",
  "log2(TPM + 1) was generated for visualization only.",
  "Tau and SPM were calculated from untransformed median TMM.TPM."
)

writeLines(
  method_lines,
  con = file.path(
    outdir,
    "07_Tau_SPM_method_definition.txt"
  )
)


# ============================================================
# 25. Session information
# ============================================================

sink(
  file.path(
    outdir,
    "07_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 26. Final report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("07 TMM.TPM / TAU / SPM ANNOTATION COMPLETED\n")
cat("============================================================\n")

cat(
  "Genes analyzed                         : ",
  nrow(tpm_mat),
  "\n",
  sep = ""
)

cat(
  "Samples analyzed                       : ",
  ncol(tpm_mat),
  "\n",
  sep = ""
)

cat("------------------------------------------------------------\n")

cat(
  "LRT dynamic genes                      : ",
  length(dynamic_ids),
  "\n",
  sep = ""
)

cat(
  "High-confidence dynamic genes          : ",
  length(high_conf_ids),
  "\n",
  sep = ""
)

cat("------------------------------------------------------------\n")

cat("ALL 11 STAGES:\n")

cat(
  "  Tau >= 0.85                          : ",
  n_tau_all11,
  "\n",
  sep = ""
)

cat(
  "  Peak SPM >= 0.50                     : ",
  n_spm_all11,
  "\n",
  sep = ""
)

cat(
  "  Tau >= 0.85 AND Peak SPM >= 0.50    : ",
  n_specific_all11,
  "\n",
  sep = ""
)

cat("------------------------------------------------------------\n")

cat("INJURY-ONLY 10 STAGES:\n")

cat(
  "  Tau >= 0.85                          : ",
  n_tau_injury10,
  "\n",
  sep = ""
)

cat(
  "  Peak SPM >= 0.50                     : ",
  n_spm_injury10,
  "\n",
  sep = ""
)

cat(
  "  Tau >= 0.85 AND Peak SPM >= 0.50    : ",
  n_specific_injury10,
  "\n",
  sep = ""
)

cat(
  "  + high-confidence dynamic            : ",
  nrow(hc_injury_specific),
  "\n",
  sep = ""
)

cat("------------------------------------------------------------\n")
cat("Tau/SPM calculated from untransformed group-median TMM.TPM.\n")
cat("log2(TPM+1) is visualization-only.\n")
cat("No CPM or logCPM was used.\n")
cat("No differential-expression testing was performed in Step 07.\n")
cat("============================================================\n\n")

cat("High-confidence injury-specific genes by peak stage:\n")

print(
  hc_stage_counts,
  row.names = FALSE
)

cat("\nOutput directory:\n")
cat(outdir, "\n\n")

cat("PASS\n")
cat(
  "Step 07 is complete. ",
  "The next phase will reanalyze the 2021 historical-reference cohort ",
  "using the same frozen count/VST/DESeq2 framework before ",
  "cross-cohort temporal concordance analysis.\n"
)
