#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1
)

suppressPackageStartupMessages({

  pkgs <- c(
    "DESeq2",
    "vegan"
  )

  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Required R package missing: ", p)
    }
  }

  library(DESeq2)
  library(vegan)
})

cat("\n")
cat("============================================================\n")
cat("10 2026 SAMPLE-EXCLUSION SENSITIVITY ANALYSIS\n")
cat("Robustness of DE and cross-cohort response concordance\n")
cat("Scientific Data rat skeletal muscle transcriptome project\n")
cat("============================================================\n\n")


# ============================================================
# 0. Fixed parameters
# ============================================================

PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD <- 1

EXPECTED_GENES_2026 <- 18364
EXPECTED_COMMON_GENES <- 16918

N_PERMUTATIONS <- 9999
RANDOM_SEED <- 20260910

group_levels <- c(
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

injury_groups_2026 <- group_levels[
  group_levels != "Baseline"
]

groups_2021 <- c(
  "R4h",
  "R8h",
  "R12h",
  "R16h",
  "R20h",
  "R24h",
  "R48h"
)


# ============================================================
# 1. Input/output
# ============================================================

count_file <- paste0(
  "02_count_filtering_and_VST_result/",
  "2026_primary/",
  "02_2026_filtered_counts.tsv"
)

vst_file <- paste0(
  "02_count_filtering_and_VST_result/",
  "2026_primary/",
  "02_2026_VST_matrix.tsv"
)

meta_file <- paste0(
  "02_count_filtering_and_VST_result/",
  "2026_primary/",
  "02_2026_metadata_aligned.tsv"
)

review_file <- paste0(
  "03_2026_internal_expression_structure_result/",
  "03_2026_expression_outlier_diagnostics.tsv"
)

primary_lfc_file <- paste0(
  "05_2026_timepoint_vs_baseline_DESeq2_result/",
  "05_2026_genomewide_log2FC_matrix.tsv"
)

primary_padj_file <- paste0(
  "05_2026_timepoint_vs_baseline_DESeq2_result/",
  "05_2026_genomewide_PADJ_matrix.tsv"
)

lfc_2021_file <- paste0(
  "08_2021_historical_reference_reanalysis_result/",
  "08_2021_genomewide_log2FC_matrix.tsv"
)

primary_cross_file <- paste0(
  "09_cross_cohort_temporal_response_concordance_result/",
  "09_cross_cohort_Spearman_matrix.tsv"
)

outdir <- "10_2026_sample_exclusion_sensitivity_result"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Check files
# ============================================================

required_files <- c(
  count_file,
  vst_file,
  meta_file,
  review_file,
  primary_lfc_file,
  primary_padj_file,
  lfc_2021_file,
  primary_cross_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  cat("Missing input files:\n")

  for (f in missing_files) {
    cat("  -", f, "\n")
  }

  quit(status = 1)
}

cat("All required input files found: PASS\n\n")


# ============================================================
# 3. Generic gene-matrix reader
# ============================================================

read_gene_matrix <- function(file) {

  x <- read.delim(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )

  colnames(x)[1] <- "GeneID"

  if (anyDuplicated(x$GeneID)) {
    stop("Duplicated GeneID in ", file)
  }

  x
}


# ============================================================
# 4. Read 2026 counts
# ============================================================

cat("Reading 2026 filtered counts...\n")

count_df <- read_gene_matrix(
  count_file
)

count_mat <- as.matrix(
  count_df[, -1, drop = FALSE]
)

storage.mode(count_mat) <- "numeric"

rownames(count_mat) <- count_df$GeneID

if (anyNA(count_mat) ||
    any(!is.finite(count_mat)) ||
    any(count_mat < 0)) {

  stop("Invalid count values.")
}

if (any(
  abs(count_mat - round(count_mat)) > 1e-8
)) {
  stop("Noninteger counts detected.")
}

count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"


# ============================================================
# 5. Read VST
# ============================================================

cat("Reading 2026 VST matrix...\n")

vst_df <- read_gene_matrix(
  vst_file
)

vst_mat <- as.matrix(
  vst_df[, -1, drop = FALSE]
)

storage.mode(vst_mat) <- "numeric"

rownames(vst_mat) <- vst_df$GeneID


# ============================================================
# 6. Explicit count/VST alignment
# ============================================================

if (!setequal(
  rownames(count_mat),
  rownames(vst_mat)
)) {
  stop("Count/VST GeneID sets differ.")
}

if (!setequal(
  colnames(count_mat),
  colnames(vst_mat)
)) {
  stop("Count/VST sample sets differ.")
}

vst_mat <- vst_mat[
  rownames(count_mat),
  colnames(count_mat),
  drop = FALSE
]


# ============================================================
# 7. Read metadata
# ============================================================

meta <- read.delim(
  meta_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!setequal(
  colnames(count_mat),
  meta$sampleID
)) {
  stop("Count/metadata sample sets differ.")
}

meta <- meta[
  match(
    colnames(count_mat),
    meta$sampleID
  ),
  ,
  drop = FALSE
]

rownames(meta) <- meta$sampleID

meta$group <- factor(
  meta$group,
  levels = group_levels
)


# ============================================================
# 8. Validate gene number
# ============================================================

if (nrow(count_mat) != EXPECTED_GENES_2026) {
  stop(
    "Expected ",
    EXPECTED_GENES_2026,
    " genes, observed ",
    nrow(count_mat)
  )
}

cat(
  "2026 fixed gene universe : ",
  nrow(count_mat),
  "\n\n",
  sep = ""
)


# ============================================================
# 9. Read Step 03 expression-review results
# ============================================================

review <- read.delim(
  review_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_review_cols <- c(
  "sampleID",
  "Review_level"
)

if (!all(
  required_review_cols %in%
    colnames(review)
)) {
  stop(
    "Expected Step 03 review columns missing."
  )
}

multi_metric_samples <- review$sampleID[
  review$Review_level ==
    "Multi_metric_review"
]

multi_metric_samples <- unique(
  multi_metric_samples
)

cat(
  "Step 03 multi-metric review samples: ",
  length(multi_metric_samples),
  "\n",
  sep = ""
)

if (length(multi_metric_samples) > 0) {

  cat(
    "  ",
    paste(
      multi_metric_samples,
      collapse = ", "
    ),
    "\n\n",
    sep = ""
  )
}


# ============================================================
# 10. Prior raw-QC review sample
# ============================================================

prior_qc_review <- "R2h_4"

if (!prior_qc_review %in%
    meta$sampleID) {

  warning(
    prior_qc_review,
    " was not found and will be ignored."
  )

  prior_qc_samples <- character(0)

} else {

  prior_qc_samples <- prior_qc_review
}

cat(
  "Prior raw-QC review sample     : ",
  ifelse(
    length(prior_qc_samples) == 0,
    "None",
    paste(
      prior_qc_samples,
      collapse = ", "
    )
  ),
  "\n\n",
  sep = ""
)


# ============================================================
# 11. Define scenarios
# ============================================================

scenario_exclusions <- list(

  Primary_all77 =
    character(0),

  Exclude_multi_metric =
    multi_metric_samples,

  Exclude_R2h4 =
    prior_qc_samples,

  Exclude_union =
    unique(
      c(
        multi_metric_samples,
        prior_qc_samples
      )
    )
)


# ============================================================
# 12. Save scenario definitions
# ============================================================

scenario_definition <- do.call(
  rbind,
  lapply(
    names(scenario_exclusions),
    function(s) {

      exc <- scenario_exclusions[[s]]

      data.frame(
        Scenario = s,
        Excluded_sample_number =
          length(exc),
        Excluded_samples =
          ifelse(
            length(exc) == 0,
            "None",
            paste(
              exc,
              collapse = ";"
            )
          ),
        Retained_samples =
          nrow(meta) - length(exc),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.table(
  scenario_definition,
  file = file.path(
    outdir,
    "10_sensitivity_scenario_definitions.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Sensitivity scenarios:\n")
print(
  scenario_definition,
  row.names = FALSE
)
cat("\n")


# ============================================================
# 13. Read frozen primary response matrices
# ============================================================

primary_lfc_df <- read_gene_matrix(
  primary_lfc_file
)

primary_padj_df <- read_gene_matrix(
  primary_padj_file
)

primary_lfc <- as.matrix(
  primary_lfc_df[
    ,
    injury_groups_2026,
    drop = FALSE
  ]
)

primary_padj <- as.matrix(
  primary_padj_df[
    ,
    injury_groups_2026,
    drop = FALSE
  ]
)

storage.mode(primary_lfc) <- "numeric"
storage.mode(primary_padj) <- "numeric"

rownames(primary_lfc) <-
  primary_lfc_df$GeneID

rownames(primary_padj) <-
  primary_padj_df$GeneID

primary_lfc <- primary_lfc[
  rownames(count_mat),
  ,
  drop = FALSE
]

primary_padj <- primary_padj[
  rownames(count_mat),
  ,
  drop = FALSE
]


# ============================================================
# 14. Read 2021 response matrix
# ============================================================

lfc21_df <- read_gene_matrix(
  lfc_2021_file
)

common_genes <- rownames(count_mat)[
  rownames(count_mat) %in%
    lfc21_df$GeneID
]

if (length(common_genes) !=
    EXPECTED_COMMON_GENES) {

  stop(
    "Expected ",
    EXPECTED_COMMON_GENES,
    " common genes, observed ",
    length(common_genes)
  )
}

lfc21 <- as.matrix(
  lfc21_df[
    match(
      common_genes,
      lfc21_df$GeneID
    ),
    groups_2021,
    drop = FALSE
  ]
)

storage.mode(lfc21) <- "numeric"

rownames(lfc21) <- common_genes


# ============================================================
# 15. Read primary Step 09 cross-cohort matrix
# ============================================================

primary_cross_df <- read.delim(
  primary_cross_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rownames(primary_cross_df) <-
  primary_cross_df[[1]]

primary_cross <- as.matrix(
  primary_cross_df[
    ,
    -1,
    drop = FALSE
  ]
)

storage.mode(primary_cross) <- "numeric"

primary_cross <- primary_cross[
  injury_groups_2026,
  groups_2021,
  drop = FALSE
]


# ============================================================
# 16. Storage objects
# ============================================================

scenario_lfc_list <- list()
scenario_padj_list <- list()

de_robustness_list <- list()
cross_robustness_list <- list()
global_structure_list <- list()


# ============================================================
# 17. Analyze each sensitivity scenario
# ============================================================

for (scenario in
     names(scenario_exclusions)) {

  cat(
    "\n============================================================\n"
  )

  cat(
    "Scenario: ",
    scenario,
    "\n",
    sep = ""
  )

  cat(
    "============================================================\n"
  )

  excluded <- scenario_exclusions[[scenario]]

  retained_samples <- setdiff(
    meta$sampleID,
    excluded
  )

  meta_s <- meta[
    retained_samples,
    ,
    drop = FALSE
  ]

  count_s <- count_mat[
    ,
    retained_samples,
    drop = FALSE
  ]

  vst_s <- vst_mat[
    ,
    retained_samples,
    drop = FALSE
  ]

  meta_s$group <- droplevels(
    factor(
      meta_s$group,
      levels = group_levels
    )
  )


  # ----------------------------------------------------------
  # Check all groups retained
  # ----------------------------------------------------------

  group_counts <- table(
    factor(
      meta_s$group,
      levels = group_levels
    )
  )

  if (any(group_counts == 0)) {
    stop(
      "A biological stage has no samples in scenario ",
      scenario
    )
  }

  cat(
    "Retained samples : ",
    ncol(count_s),
    "\n",
    sep = ""
  )

  cat(
    "Group sizes      : ",
    paste(
      paste0(
        names(group_counts),
        "=",
        as.integer(group_counts)
      ),
      collapse = ", "
    ),
    "\n",
    sep = ""
  )


  # ==========================================================
  # 17A. Global temporal structure
  # ==========================================================

  sample_distance <- dist(
    t(vst_s),
    method = "euclidean"
  )

  set.seed(RANDOM_SEED)

  perm <- vegan::adonis2(
    sample_distance ~ group,
    data = meta_s,
    permutations = N_PERMUTATIONS
  )

  bd <- vegan::betadisper(
    sample_distance,
    group = meta_s$group,
    type = "median",
    bias.adjust = TRUE
  )

  set.seed(RANDOM_SEED)

  bd_perm <- vegan::permutest(
    bd,
    permutations = N_PERMUTATIONS
  )

  perm_r2 <- as.numeric(
    perm[
      "Model",
      "R2"
    ]
  )

  perm_f <- as.numeric(
    perm[
      "Model",
      "F"
    ]
  )

  perm_p <- as.numeric(
    perm[
      "Model",
      "Pr(>F)"
    ]
  )

  permdisp_p <- as.numeric(
    bd_perm$tab[
      "Groups",
      "Pr(>F)"
    ]
  )

  global_structure_list[[scenario]] <- data.frame(
    Scenario = scenario,
    Retained_samples = ncol(count_s),
    Excluded_samples = length(excluded),
    PERMANOVA_R2 = perm_r2,
    PERMANOVA_F = perm_f,
    PERMANOVA_P = perm_p,
    PERMDISP_P = permdisp_p,
    stringsAsFactors = FALSE
  )


  # ==========================================================
  # 17B. DESeq2
  # ==========================================================

  cat("Running DESeq2...\n")

  dds <- DESeqDataSetFromMatrix(
    countData = count_s,
    colData = meta_s,
    design = ~ group
  )

  dds$group <- relevel(
    dds$group,
    ref = "Baseline"
  )

  dds <- DESeq(
    dds,
    test = "Wald",
    minReplicatesForReplace = Inf,
    quiet = TRUE
  )


  scenario_lfc <- matrix(
    NA_real_,
    nrow = nrow(count_s),
    ncol = length(
      injury_groups_2026
    ),
    dimnames = list(
      rownames(count_s),
      injury_groups_2026
    )
  )

  scenario_padj <- scenario_lfc


  # ==========================================================
  # 17C. Per-stage DE robustness
  # ==========================================================

  for (g in injury_groups_2026) {

    res <- results(
      dds,
      contrast = c(
        "group",
        g,
        "Baseline"
      ),
      alpha = PADJ_THRESHOLD,
      independentFiltering = FALSE
    )

    x <- as.data.frame(res)

    scenario_lfc[, g] <-
      x[
        rownames(count_s),
        "log2FoldChange"
      ]

    scenario_padj[, g] <-
      x[
        rownames(count_s),
        "padj"
      ]


    # Primary vs sensitivity response
    a <- primary_lfc[, g]
    b <- scenario_lfc[, g]

    complete <- (
      is.finite(a) &
      is.finite(b)
    )

    response_rho <- cor(
      a[complete],
      b[complete],
      method = "spearman"
    )

    response_r <- cor(
      a[complete],
      b[complete],
      method = "pearson"
    )

    median_abs_lfc_delta <- median(
      abs(
        a[complete] -
        b[complete]
      ),
      na.rm = TRUE
    )

    primary_deg <- (
      !is.na(primary_padj[, g]) &
      primary_padj[, g] <
        PADJ_THRESHOLD &
      abs(primary_lfc[, g]) >=
        LFC_THRESHOLD
    )

    sensitivity_deg <- (
      !is.na(scenario_padj[, g]) &
      scenario_padj[, g] <
        PADJ_THRESHOLD &
      abs(scenario_lfc[, g]) >=
        LFC_THRESHOLD
    )

    shared_deg <- (
      primary_deg &
      sensitivity_deg
    )

    n_primary <- sum(primary_deg)
    n_sens <- sum(sensitivity_deg)
    n_shared <- sum(shared_deg)

    if (n_shared > 0) {

      same_direction <- sum(
        sign(
          primary_lfc[
            shared_deg,
            g
          ]
        ) ==
        sign(
          scenario_lfc[
            shared_deg,
            g
          ]
        )
      )

      direction_concordance <-
        same_direction /
        n_shared

    } else {

      direction_concordance <-
        NA_real_
    }

    de_robustness_list[[length(de_robustness_list) + 1]] <- data.frame(
      Scenario = scenario,
      Group = g,
      Primary_DEG = n_primary,
      Sensitivity_DEG = n_sens,
      DEG_number_difference =
        n_sens - n_primary,
      Shared_DEG = n_shared,
      Shared_DEG_direction_concordance =
        direction_concordance,
      Genomewide_log2FC_Spearman =
        response_rho,
      Genomewide_log2FC_Pearson =
        response_r,
      Median_absolute_log2FC_difference =
        median_abs_lfc_delta,
      stringsAsFactors = FALSE
    )
  }


  scenario_lfc_list[[scenario]] <- scenario_lfc

  scenario_padj_list[[scenario]] <- scenario_padj


  # ==========================================================
  # 17D. Cross-cohort concordance matrix
  # ==========================================================

  scenario_common_lfc <-
    scenario_lfc[
      common_genes,
      ,
      drop = FALSE
    ]

  cross_mat <- matrix(
    NA_real_,
    nrow = length(
      injury_groups_2026
    ),
    ncol = length(
      groups_2021
    ),
    dimnames = list(
      injury_groups_2026,
      groups_2021
    )
  )

  for (g26 in
       injury_groups_2026) {

    for (g21 in
         groups_2021) {

      a <- scenario_common_lfc[
        ,
        g26
      ]

      b <- lfc21[
        ,
        g21
      ]

      complete <- (
        is.finite(a) &
        is.finite(b)
      )

      cross_mat[
        g26,
        g21
      ] <- cor(
        a[complete],
        b[complete],
        method = "spearman"
      )
    }
  }


  # ----------------------------------------------------------
  # Save each scenario matrix
  # ----------------------------------------------------------

  write.table(
    data.frame(
      Cohort2026_stage =
        rownames(cross_mat),
      cross_mat,
      check.names = FALSE
    ),
    file = file.path(
      outdir,
      paste0(
        "10_",
        scenario,
        "_cross_cohort_Spearman_matrix.tsv"
      )
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )


  # ----------------------------------------------------------
  # Compare 70-cell matrix against frozen primary Step 09
  # ----------------------------------------------------------

  primary_vec <- as.vector(
    primary_cross
  )

  scenario_vec <- as.vector(
    cross_mat
  )

  matrix_cell_correlation <- cor(
    primary_vec,
    scenario_vec,
    method = "spearman"
  )

  matrix_mae <- mean(
    abs(
      primary_vec -
      scenario_vec
    ),
    na.rm = TRUE
  )

  matrix_max_delta <- max(
    abs(
      primary_vec -
      scenario_vec
    ),
    na.rm = TRUE
  )

  cross_robustness_list[[scenario]] <- data.frame(
    Scenario = scenario,

    Mean_cross_cohort_rho =
      mean(
        scenario_vec,
        na.rm = TRUE
      ),

    Median_cross_cohort_rho =
      median(
        scenario_vec,
        na.rm = TRUE
      ),

    Maximum_cross_cohort_rho =
      max(
        scenario_vec,
        na.rm = TRUE
      ),

    Correlation_with_primary_70cell_matrix =
      matrix_cell_correlation,

    Mean_absolute_change_from_primary =
      matrix_mae,

    Maximum_absolute_change_from_primary =
      matrix_max_delta,

    stringsAsFactors = FALSE
  )

  cat(
    "PERMANOVA R2 : ",
    sprintf(
      "%.4f",
      perm_r2
    ),
    "\n",
    sep = ""
  )

  cat(
    "PERMDISP P   : ",
    format(
      permdisp_p,
      scientific = TRUE
    ),
    "\n",
    sep = ""
  )

  cat(
    "Cross-cohort median rho : ",
    sprintf(
      "%.4f",
      median(
        scenario_vec,
        na.rm = TRUE
      )
    ),
    "\n",
    sep = ""
  )
}


# ============================================================
# 18. Combine outputs
# ============================================================

global_structure <- do.call(
  rbind,
  global_structure_list
)

de_robustness <- do.call(
  rbind,
  de_robustness_list
)

cross_robustness <- do.call(
  rbind,
  cross_robustness_list
)

rownames(global_structure) <- NULL
rownames(de_robustness) <- NULL
rownames(cross_robustness) <- NULL


# ============================================================
# 19. Save outputs
# ============================================================

write.table(
  global_structure,
  file = file.path(
    outdir,
    "10_global_structure_sensitivity.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  de_robustness,
  file = file.path(
    outdir,
    "10_DESeq2_response_sensitivity.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  cross_robustness,
  file = file.path(
    outdir,
    "10_cross_cohort_concordance_sensitivity.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 20. Nearest-time descriptive comparison
# ============================================================

cross_long_file <- paste0(
  "09_cross_cohort_temporal_response_concordance_result/",
  "09_cross_cohort_pairwise_concordance_long.tsv"
)

cross_long <- read.delim(
  cross_long_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

nearest_list <- list()

k <- 1

for (g21 in groups_2021) {

  x <- cross_long[
    cross_long$Cohort2021_group ==
      g21,
    ,
    drop = FALSE
  ]

  min_gap <- min(
    x$Absolute_time_difference_h
  )

  nearest <- x[
    x$Absolute_time_difference_h ==
      min_gap,
    ,
    drop = FALSE
  ]

  nearest_list[[k]] <- nearest

  k <- k + 1
}

nearest_time_pairs <- do.call(
  rbind,
  nearest_list
)

rownames(nearest_time_pairs) <- NULL

write.table(
  nearest_time_pairs,
  file = file.path(
    outdir,
    "10_nearest_time_cross_cohort_pairs.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 21. Time-gap bins
# ============================================================

cross_long$Time_gap_bin <- cut(
  cross_long$Absolute_time_difference_h,
  breaks = c(
    -Inf,
    2,
    6,
    12,
    Inf
  ),
  labels = c(
    "<=2 h",
    ">2-6 h",
    ">6-12 h",
    ">12 h"
  ),
  right = TRUE
)

time_gap_summary <- do.call(
  rbind,
  lapply(
    levels(
      cross_long$Time_gap_bin
    ),
    function(bin) {

      x <- cross_long[
        cross_long$Time_gap_bin ==
          bin,
        ,
        drop = FALSE
      ]

      data.frame(
        Time_gap_bin = bin,
        Pair_number = nrow(x),
        Mean_Spearman_rho =
          mean(
            x$Spearman_rho,
            na.rm = TRUE
          ),
        Median_Spearman_rho =
          median(
            x$Spearman_rho,
            na.rm = TRUE
          ),
        Mean_direction_concordance =
          mean(
            x$Shared_DEG_direction_concordance,
            na.rm = TRUE
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.table(
  time_gap_summary,
  file = file.path(
    outdir,
    "10_cross_cohort_time_gap_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 22. Save method note
# ============================================================

method_note <- c(
  "Step 10 sample-exclusion sensitivity analysis",
  "",
  "The Step 02 gene universe was held fixed.",
  "Genes were not re-filtered after sample exclusion.",
  "",
  "This ensures that sensitivity comparisons isolate",
  "the effect of sample inclusion/exclusion rather than",
  "changes in the analyzed gene universe.",
  "",
  "Primary scenario:",
  "all 77 samples.",
  "",
  "Sensitivity A:",
  "exclude Step 03 Multi_metric_review samples.",
  "",
  "Sensitivity B:",
  "exclude the previously reviewed raw-QC sample R2h_4.",
  "",
  "Sensitivity C:",
  "exclude the union of the above samples.",
  "",
  "For each scenario:",
  "1. DESeq2 Wald models were refitted.",
  "2. Automatic Cook's count replacement remained disabled.",
  "3. Baseline-relative genome-wide log2FC values were recalculated.",
  "4. Global VST PERMANOVA/PERMDISP was recalculated.",
  "5. 2026-vs-2021 genome-wide response concordance was recalculated.",
  "",
  "Nearest-time cross-cohort comparisons are descriptive.",
  "They are not treated as exact experimental time equivalents."
)

writeLines(
  method_note,
  con = file.path(
    outdir,
    "10_sensitivity_method_definition.txt"
  )
)


# ============================================================
# 23. Session info
# ============================================================

sink(
  file.path(
    outdir,
    "10_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 24. Final report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("10 SAMPLE-EXCLUSION SENSITIVITY COMPLETED\n")
cat("============================================================\n\n")

cat("Global temporal structure:\n")

print(
  global_structure,
  row.names = FALSE
)

cat("\n")
cat("Cross-cohort concordance robustness:\n")

print(
  cross_robustness,
  row.names = FALSE
)

cat("\n")
cat("Nearest-time cross-cohort pairs:\n")

print(
  nearest_time_pairs[
    ,
    c(
      "Cohort2026_stage",
      "Cohort2021_stage",
      "Absolute_time_difference_h",
      "Spearman_rho",
      "Shared_DEGs",
      "Shared_DEG_direction_concordance"
    )
  ],
  row.names = FALSE
)

cat("\n")
cat("Time-gap summary:\n")

print(
  time_gap_summary,
  row.names = FALSE
)

cat("\n")
cat("Output directory:\n")
cat(outdir, "\n\n")

cat("PASS\n")

cat(
  "Step 10 is complete. ",
  "If sensitivity results remain stable, ",
  "the analytical framework can be frozen for ",
  "final Technical Validation figure/table assembly.\n"
)
