#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(DESeq2)
  library(vegan)
})

cat("\n")
cat("============================================================\n")
cat("10B SINGLE-SAMPLE LEAVE-ONE-OUT INFLUENCE ANALYSIS\n")
cat("============================================================\n\n")

PADJ_CUTOFF <- 0.05
LFC_CUTOFF <- 1
N_PERM <- 9999
SEED <- 20260910

groups26 <- c(
  "R0h","R1h","R2h","R6h","R10h",
  "R14h","R18h","R36h","R60h","R72h"
)

groups21 <- c(
  "R4h","R8h","R12h","R16h",
  "R20h","R24h","R48h"
)

group_levels <- c("Baseline", groups26)

count_file <- paste0(
  "02_count_filtering_and_VST_result/2026_primary/",
  "02_2026_filtered_counts.tsv"
)

vst_file <- paste0(
  "02_count_filtering_and_VST_result/2026_primary/",
  "02_2026_VST_matrix.tsv"
)

meta_file <- paste0(
  "02_count_filtering_and_VST_result/2026_primary/",
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

lfc21_file <- paste0(
  "08_2021_historical_reference_reanalysis_result/",
  "08_2021_genomewide_log2FC_matrix.tsv"
)

primary_cross_file <- paste0(
  "09_cross_cohort_temporal_response_concordance_result/",
  "09_cross_cohort_Spearman_matrix.tsv"
)

outdir <- "10B_single_sample_leave_one_out_result"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


read_gene_matrix <- function(file) {

  x <- read.delim(
    file,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  colnames(x)[1] <- "GeneID"

  if (anyDuplicated(x$GeneID)) {
    stop("Duplicated GeneID: ", file)
  }

  x
}


# ============================================================
# Read counts
# ============================================================

count_df <- read_gene_matrix(count_file)

count_mat <- as.matrix(
  count_df[, -1, drop = FALSE]
)

storage.mode(count_mat) <- "numeric"
rownames(count_mat) <- count_df$GeneID

count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"


# ============================================================
# Read fixed VST
# ============================================================

vst_df <- read_gene_matrix(vst_file)

vst_mat <- as.matrix(
  vst_df[, -1, drop = FALSE]
)

storage.mode(vst_mat) <- "numeric"
rownames(vst_mat) <- vst_df$GeneID

vst_mat <- vst_mat[
  rownames(count_mat),
  colnames(count_mat),
  drop = FALSE
]


# ============================================================
# Metadata
# ============================================================

meta <- read.delim(
  meta_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

meta <- meta[
  match(colnames(count_mat), meta$sampleID),
  ,
  drop = FALSE
]

rownames(meta) <- meta$sampleID

meta$group <- factor(
  meta$group,
  levels = group_levels
)


# ============================================================
# Review samples
# ============================================================

review <- read.delim(
  review_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

review_samples <- unique(
  review$sampleID[
    review$Review_level == "Multi_metric_review"
  ]
)

test_samples <- unique(
  c(
    review_samples,
    "R2h_4"
  )
)

test_samples <- test_samples[
  test_samples %in% meta$sampleID
]

cat(
  "Samples tested:\n",
  paste(test_samples, collapse = ", "),
  "\n\n"
)


# ============================================================
# Primary 2026 matrices
# ============================================================

primary_lfc_df <- read_gene_matrix(
  primary_lfc_file
)

primary_padj_df <- read_gene_matrix(
  primary_padj_file
)

idx <- match(
  rownames(count_mat),
  primary_lfc_df$GeneID
)

primary_lfc <- as.matrix(
  primary_lfc_df[
    idx,
    groups26,
    drop = FALSE
  ]
)

storage.mode(primary_lfc) <- "numeric"
rownames(primary_lfc) <- rownames(count_mat)

idx <- match(
  rownames(count_mat),
  primary_padj_df$GeneID
)

primary_padj <- as.matrix(
  primary_padj_df[
    idx,
    groups26,
    drop = FALSE
  ]
)

storage.mode(primary_padj) <- "numeric"
rownames(primary_padj) <- rownames(count_mat)


# ============================================================
# 2021 matrix
# ============================================================

lfc21_df <- read_gene_matrix(lfc21_file)

common_genes <- intersect(
  rownames(count_mat),
  lfc21_df$GeneID
)

if (length(common_genes) != 16918) {
  stop(
    "Expected 16918 common genes; observed ",
    length(common_genes)
  )
}

lfc21 <- as.matrix(
  lfc21_df[
    match(common_genes, lfc21_df$GeneID),
    groups21,
    drop = FALSE
  ]
)

storage.mode(lfc21) <- "numeric"
rownames(lfc21) <- common_genes


# ============================================================
# Primary cross-cohort matrix
# ============================================================

primary_cross_df <- read.delim(
  primary_cross_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rownames(primary_cross_df) <-
  primary_cross_df[[1]]

primary_cross <- as.matrix(
  primary_cross_df[
    groups26,
    groups21,
    drop = FALSE
  ]
)

storage.mode(primary_cross) <- "numeric"


# ============================================================
# Storage
# ============================================================

overall_list <- list()
detail_list <- list()


# ============================================================
# Leave one sample out
# ============================================================

for (removed_sample in test_samples) {

  cat(
    "============================================================\n"
  )

  cat(
    "Removing sample: ",
    removed_sample,
    "\n",
    sep = ""
  )

  removed_group <- as.character(
    meta[
      removed_sample,
      "group"
    ]
  )

  retained <- setdiff(
    meta$sampleID,
    removed_sample
  )

  meta_s <- meta[
    retained,
    ,
    drop = FALSE
  ]

  counts_s <- count_mat[
    ,
    retained,
    drop = FALSE
  ]

  vst_s <- vst_mat[
    ,
    retained,
    drop = FALSE
  ]

  meta_s$group <- factor(
    meta_s$group,
    levels = group_levels
  )


  # ==========================================================
  # Global temporal structure
  # ==========================================================

  d <- dist(
    t(vst_s),
    method = "euclidean"
  )

  set.seed(SEED)

  perm <- vegan::adonis2(
    d ~ group,
    data = meta_s,
    permutations = N_PERM
  )

  bd <- vegan::betadisper(
    d,
    group = meta_s$group,
    type = "median",
    bias.adjust = TRUE
  )

  set.seed(SEED)

  bd_perm <- vegan::permutest(
    bd,
    permutations = N_PERM
  )

  perm_r2 <- as.numeric(
    perm["Model", "R2"]
  )

  perm_p <- as.numeric(
    perm["Model", "Pr(>F)"]
  )

  permdisp_p <- as.numeric(
    bd_perm$tab["Groups", "Pr(>F)"]
  )


  # ==========================================================
  # Refit DESeq2
  # ==========================================================

  dds <- DESeqDataSetFromMatrix(
    countData = counts_s,
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


  loo_lfc <- matrix(
    NA_real_,
    nrow = nrow(count_mat),
    ncol = length(groups26),
    dimnames = list(
      rownames(count_mat),
      groups26
    )
  )

  loo_padj <- loo_lfc


  # ==========================================================
  # Individual timepoint influence
  # ==========================================================

  for (g in groups26) {

    res <- results(
      dds,
      contrast = c(
        "group",
        g,
        "Baseline"
      ),
      alpha = PADJ_CUTOFF,
      independentFiltering = FALSE
    )

    rdf <- as.data.frame(res)

    loo_lfc[, g] <-
      rdf[
        rownames(count_mat),
        "log2FoldChange"
      ]

    loo_padj[, g] <-
      rdf[
        rownames(count_mat),
        "padj"
      ]


    complete <- (
      is.finite(primary_lfc[, g]) &
      is.finite(loo_lfc[, g])
    )

    lfc_rho <- cor(
      primary_lfc[complete, g],
      loo_lfc[complete, g],
      method = "spearman"
    )

    lfc_r <- cor(
      primary_lfc[complete, g],
      loo_lfc[complete, g],
      method = "pearson"
    )

    median_delta <- median(
      abs(
        primary_lfc[complete, g] -
        loo_lfc[complete, g]
      )
    )

    primary_deg <- (
      !is.na(primary_padj[, g]) &
      primary_padj[, g] < PADJ_CUTOFF &
      abs(primary_lfc[, g]) >= LFC_CUTOFF
    )

    loo_deg <- (
      !is.na(loo_padj[, g]) &
      loo_padj[, g] < PADJ_CUTOFF &
      abs(loo_lfc[, g]) >= LFC_CUTOFF
    )

    shared <- primary_deg & loo_deg

    if (sum(shared) > 0) {

      direction <- mean(
        sign(primary_lfc[shared, g]) ==
        sign(loo_lfc[shared, g])
      )

    } else {

      direction <- NA_real_
    }

    detail_list[[length(detail_list) + 1]] <-
      data.frame(
        Removed_sample = removed_sample,
        Removed_group = removed_group,
        Contrast = g,
        Primary_DEG = sum(primary_deg),
        Leave_one_out_DEG = sum(loo_deg),
        DEG_number_difference =
          sum(loo_deg) - sum(primary_deg),
        Shared_DEG = sum(shared),
        Shared_DEG_direction_concordance =
          direction,
        Genomewide_log2FC_Spearman =
          lfc_rho,
        Genomewide_log2FC_Pearson =
          lfc_r,
        Median_absolute_log2FC_difference =
          median_delta,
        stringsAsFactors = FALSE
      )
  }


  # ==========================================================
  # Recalculate cross-cohort 10 x 7 matrix
  # ==========================================================

  loo_common <- loo_lfc[
    common_genes,
    ,
    drop = FALSE
  ]

  cross_mat <- matrix(
    NA_real_,
    nrow = length(groups26),
    ncol = length(groups21),
    dimnames = list(
      groups26,
      groups21
    )
  )

  for (g26 in groups26) {

    for (g21 in groups21) {

      a <- loo_common[, g26]
      b <- lfc21[, g21]

      ok <- is.finite(a) & is.finite(b)

      cross_mat[g26, g21] <- cor(
        a[ok],
        b[ok],
        method = "spearman"
      )
    }
  }


  primary_vec <- as.vector(
    primary_cross
  )

  loo_vec <- as.vector(
    cross_mat
  )

  cross_matrix_rho <- cor(
    primary_vec,
    loo_vec,
    method = "spearman"
  )

  cross_mae <- mean(
    abs(
      loo_vec -
      primary_vec
    )
  )

  cross_max_delta <- max(
    abs(
      loo_vec -
      primary_vec
    )
  )


  # ==========================================================
  # Overall per-sample influence summary
  # ==========================================================

  this_detail <- do.call(
    rbind,
    detail_list
  )

  this_detail <- this_detail[
    this_detail$Removed_sample ==
      removed_sample,
    ,
    drop = FALSE
  ]

  overall_list[[length(overall_list) + 1]] <-
    data.frame(
      Removed_sample = removed_sample,
      Removed_group = removed_group,
      Retained_samples = length(retained),
      Remaining_samples_in_removed_group =
        sum(meta_s$group == removed_group),

      PERMANOVA_R2 = perm_r2,
      PERMANOVA_P = perm_p,
      PERMDISP_P = permdisp_p,

      Minimum_genomewide_log2FC_Spearman =
        min(
          this_detail$Genomewide_log2FC_Spearman
        ),

      Median_genomewide_log2FC_Spearman =
        median(
          this_detail$Genomewide_log2FC_Spearman
        ),

      Maximum_median_absolute_log2FC_difference =
        max(
          this_detail$Median_absolute_log2FC_difference
        ),

      Cross_cohort_70cell_matrix_Spearman =
        cross_matrix_rho,

      Cross_cohort_matrix_MAE =
        cross_mae,

      Cross_cohort_matrix_maximum_change =
        cross_max_delta,

      stringsAsFactors = FALSE
    )

  cat(
    "PERMANOVA R2                       : ",
    sprintf("%.4f", perm_r2),
    "\n",
    sep = ""
  )

  cat(
    "PERMDISP P                         : ",
    format(permdisp_p, scientific = TRUE),
    "\n",
    sep = ""
  )

  cat(
    "Minimum primary-vs-LOO log2FC rho  : ",
    sprintf(
      "%.4f",
      min(this_detail$Genomewide_log2FC_Spearman)
    ),
    "\n",
    sep = ""
  )

  cat(
    "70-cell cross-matrix correlation   : ",
    sprintf("%.4f", cross_matrix_rho),
    "\n\n",
    sep = ""
  )
}


# ============================================================
# Combine and save
# ============================================================

overall <- do.call(
  rbind,
  overall_list
)

detail <- do.call(
  rbind,
  detail_list
)

rownames(overall) <- NULL
rownames(detail) <- NULL

write.table(
  overall,
  file = file.path(
    outdir,
    "10B_single_sample_influence_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  detail,
  file = file.path(
    outdir,
    "10B_single_sample_DESeq2_influence_detail.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# Session info
# ============================================================

sink(
  file.path(
    outdir,
    "10B_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# Final report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("10B SINGLE-SAMPLE INFLUENCE COMPLETED\n")
cat("============================================================\n\n")

print(
  overall,
  row.names = FALSE
)

cat("\nOutput directory:\n")
cat(outdir, "\n\n")

cat("PASS\n")
