#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  bitmapType = "cairo"
)

suppressPackageStartupMessages({

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("R package 'pheatmap' is required.")
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("R package 'ggplot2' is required.")
  }

  library(pheatmap)
  library(ggplot2)
})

cat("\n")
cat("============================================================\n")
cat("09 CROSS-COHORT TEMPORAL RESPONSE CONCORDANCE\n")
cat("2026 primary cohort versus 2021 historical-reference cohort\n")
cat("Baseline-relative genome-wide log2FC comparison\n")
cat("Scientific Data rat skeletal muscle transcriptome project\n")
cat("============================================================\n\n")


# ============================================================
# 0. Fixed parameters
# ============================================================

EXPECTED_COMMON_GENES <- 16918

PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD <- 1

groups_2026 <- c(
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

groups_2021 <- c(
  "R4h",
  "R8h",
  "R12h",
  "R16h",
  "R20h",
  "R24h",
  "R48h"
)

times_2026 <- c(
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

times_2021 <- c(
  R4h = 4,
  R8h = 8,
  R12h = 12,
  R16h = 16,
  R20h = 20,
  R24h = 24,
  R48h = 48
)

labels_2026 <- c(
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

labels_2021 <- c(
  R4h = "4 h",
  R8h = "8 h",
  R12h = "12 h",
  R16h = "16 h",
  R20h = "20 h",
  R24h = "24 h",
  R48h = "48 h"
)


# ============================================================
# 1. Input/output
# ============================================================

logfc_2026_file <- paste0(
  "05_2026_timepoint_vs_baseline_DESeq2_result/",
  "05_2026_genomewide_log2FC_matrix.tsv"
)

padj_2026_file <- paste0(
  "05_2026_timepoint_vs_baseline_DESeq2_result/",
  "05_2026_genomewide_PADJ_matrix.tsv"
)

logfc_2021_file <- paste0(
  "08_2021_historical_reference_reanalysis_result/",
  "08_2021_genomewide_log2FC_matrix.tsv"
)

padj_2021_file <- paste0(
  "08_2021_historical_reference_reanalysis_result/",
  "08_2021_genomewide_PADJ_matrix.tsv"
)

outdir <- "09_cross_cohort_temporal_response_concordance_result"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. File checks
# ============================================================

required_files <- c(
  logfc_2026_file,
  padj_2026_file,
  logfc_2021_file,
  padj_2021_file
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

cat("All four input files found: PASS\n\n")


# ============================================================
# 3. Generic matrix reader
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
    stop(
      "Duplicated GeneIDs in: ",
      file
    )
  }

  if (any(
    is.na(x$GeneID) |
    trimws(x$GeneID) == ""
  )) {
    stop(
      "Missing GeneIDs in: ",
      file
    )
  }

  return(x)
}


# ============================================================
# 4. Read matrices
# ============================================================

cat("Reading 2026 response matrices...\n")

lfc26_df <- read_gene_matrix(
  logfc_2026_file
)

padj26_df <- read_gene_matrix(
  padj_2026_file
)

cat("Reading 2021 response matrices...\n")

lfc21_df <- read_gene_matrix(
  logfc_2021_file
)

padj21_df <- read_gene_matrix(
  padj_2021_file
)


# ============================================================
# 5. Validate expected timepoint columns
# ============================================================

missing_26_lfc <- setdiff(
  groups_2026,
  colnames(lfc26_df)
)

missing_26_padj <- setdiff(
  groups_2026,
  colnames(padj26_df)
)

missing_21_lfc <- setdiff(
  groups_2021,
  colnames(lfc21_df)
)

missing_21_padj <- setdiff(
  groups_2021,
  colnames(padj21_df)
)

if (length(c(
  missing_26_lfc,
  missing_26_padj,
  missing_21_lfc,
  missing_21_padj
)) > 0) {

  stop(
    "Expected timepoint columns are missing."
  )
}

cat("Timepoint-column validation: PASS\n\n")


# ============================================================
# 6. Validate within-cohort GeneID identity
# ============================================================

if (!setequal(
  lfc26_df$GeneID,
  padj26_df$GeneID
)) {
  stop(
    "2026 log2FC and PADJ GeneID sets differ."
  )
}

if (!setequal(
  lfc21_df$GeneID,
  padj21_df$GeneID
)) {
  stop(
    "2021 log2FC and PADJ GeneID sets differ."
  )
}


# ============================================================
# 7. Define common genes
# ============================================================

common_genes <- lfc26_df$GeneID[
  lfc26_df$GeneID %in%
    lfc21_df$GeneID
]

cat(
  "2026 genes                    : ",
  nrow(lfc26_df),
  "\n",
  sep = ""
)

cat(
  "2021 genes                    : ",
  nrow(lfc21_df),
  "\n",
  sep = ""
)

cat(
  "Common filtered genes         : ",
  length(common_genes),
  "\n\n",
  sep = ""
)

if (length(common_genes) !=
    EXPECTED_COMMON_GENES) {

  stop(
    "Expected ",
    EXPECTED_COMMON_GENES,
    " common genes, observed ",
    length(common_genes),
    "."
  )
}

write.table(
  data.frame(
    GeneID = common_genes,
    stringsAsFactors = FALSE
  ),
  file = file.path(
    outdir,
    "09_common_filtered_gene_list.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat(
  "Common-gene validation: ",
  EXPECTED_COMMON_GENES,
  " genes: PASS\n\n",
  sep = ""
)


# ============================================================
# 8. Align all matrices by common GeneID
# ============================================================

align_by_gene <- function(
    df,
    genes,
    value_columns) {

  idx <- match(
    genes,
    df$GeneID
  )

  if (anyNA(idx)) {
    stop(
      "Gene alignment failure."
    )
  }

  mat <- as.matrix(
    df[
      idx,
      value_columns,
      drop = FALSE
    ]
  )

  storage.mode(mat) <- "numeric"

  rownames(mat) <- genes

  return(mat)
}


lfc26 <- align_by_gene(
  lfc26_df,
  common_genes,
  groups_2026
)

padj26 <- align_by_gene(
  padj26_df,
  common_genes,
  groups_2026
)

lfc21 <- align_by_gene(
  lfc21_df,
  common_genes,
  groups_2021
)

padj21 <- align_by_gene(
  padj21_df,
  common_genes,
  groups_2021
)


# ============================================================
# 9. Validate log2FC matrices
# ============================================================

if (any(!is.finite(
  lfc26[
    !is.na(lfc26)
  ]
))) {
  stop("Non-finite 2026 log2FC detected.")
}

if (any(!is.finite(
  lfc21[
    !is.na(lfc21)
  ]
))) {
  stop("Non-finite 2021 log2FC detected.")
}

cat(
  "Gene-aligned response matrices:\n"
)

cat(
  "2026 : ",
  nrow(lfc26),
  " genes x ",
  ncol(lfc26),
  " stages\n",
  sep = ""
)

cat(
  "2021 : ",
  nrow(lfc21),
  " genes x ",
  ncol(lfc21),
  " stages\n\n",
  sep = ""
)


# ============================================================
# 10. Pairwise cross-cohort response concordance
# ============================================================

cat(
  "Calculating 10 x 7 cross-cohort response concordance...\n"
)

spearman_mat <- matrix(
  NA_real_,
  nrow = length(groups_2026),
  ncol = length(groups_2021),
  dimnames = list(
    groups_2026,
    groups_2021
  )
)

pearson_mat <- spearman_mat

n_complete_mat <- matrix(
  NA_integer_,
  nrow = length(groups_2026),
  ncol = length(groups_2021),
  dimnames = list(
    groups_2026,
    groups_2021
  )
)

shared_deg_mat <- n_complete_mat
same_direction_mat <- n_complete_mat

direction_concordance_mat <- matrix(
  NA_real_,
  nrow = length(groups_2026),
  ncol = length(groups_2021),
  dimnames = list(
    groups_2026,
    groups_2021
  )
)

long_results <- list()

counter <- 1

for (g26 in groups_2026) {

  for (g21 in groups_2021) {

    x <- lfc26[, g26]
    y <- lfc21[, g21]

    complete <- (
      is.finite(x) &
      is.finite(y)
    )

    n_complete <- sum(complete)

    if (n_complete < 100) {
      stop(
        "Unexpectedly few complete genes for ",
        g26,
        " vs ",
        g21
      )
    }

    rho <- cor(
      x[complete],
      y[complete],
      method = "spearman"
    )

    r <- cor(
      x[complete],
      y[complete],
      method = "pearson"
    )

    spearman_mat[
      g26,
      g21
    ] <- rho

    pearson_mat[
      g26,
      g21
    ] <- r

    n_complete_mat[
      g26,
      g21
    ] <- n_complete


    # --------------------------------------------------------
    # Shared DEG definition
    # --------------------------------------------------------

    sig26 <- (
      !is.na(padj26[, g26]) &
      padj26[, g26] <
        PADJ_THRESHOLD &
      abs(lfc26[, g26]) >=
        LFC_THRESHOLD
    )

    sig21 <- (
      !is.na(padj21[, g21]) &
      padj21[, g21] <
        PADJ_THRESHOLD &
      abs(lfc21[, g21]) >=
        LFC_THRESHOLD
    )

    shared_deg <- (
      sig26 &
      sig21 &
      complete
    )

    n_shared <- sum(
      shared_deg
    )

    if (n_shared > 0) {

      same_direction <- sum(
        sign(
          lfc26[
            shared_deg,
            g26
          ]
        ) ==
        sign(
          lfc21[
            shared_deg,
            g21
          ]
        )
      )

      direction_concordance <-
        same_direction /
        n_shared

    } else {

      same_direction <- 0
      direction_concordance <- NA_real_
    }

    shared_deg_mat[
      g26,
      g21
    ] <- n_shared

    same_direction_mat[
      g26,
      g21
    ] <- same_direction

    direction_concordance_mat[
      g26,
      g21
    ] <- direction_concordance


    # --------------------------------------------------------
    # Long-format record
    # --------------------------------------------------------

    long_results[[counter]] <- data.frame(
      Cohort2026_group = g26,
      Cohort2026_stage =
        unname(labels_2026[g26]),
      Cohort2026_time_h =
        unname(times_2026[g26]),

      Cohort2021_group = g21,
      Cohort2021_stage =
        unname(labels_2021[g21]),
      Cohort2021_time_h =
        unname(times_2021[g21]),

      Absolute_time_difference_h =
        abs(
          unname(times_2026[g26]) -
          unname(times_2021[g21])
        ),

      Complete_common_genes =
        n_complete,

      Spearman_rho = rho,

      Pearson_r = r,

      Shared_DEGs = n_shared,

      Shared_DEGs_same_direction =
        same_direction,

      Shared_DEG_direction_concordance =
        direction_concordance,

      stringsAsFactors = FALSE
    )

    counter <- counter + 1
  }
}

long_df <- do.call(
  rbind,
  long_results
)

rownames(long_df) <- NULL


# ============================================================
# 11. Save concordance matrices
# ============================================================

write_matrix <- function(
    mat,
    filename) {

  write.table(
    data.frame(
      Cohort2026_stage =
        rownames(mat),
      mat,
      check.names = FALSE
    ),
    file = file.path(
      outdir,
      filename
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}


write_matrix(
  spearman_mat,
  "09_cross_cohort_Spearman_matrix.tsv"
)

write_matrix(
  pearson_mat,
  "09_cross_cohort_Pearson_matrix.tsv"
)

write_matrix(
  n_complete_mat,
  "09_cross_cohort_complete_gene_number_matrix.tsv"
)

write_matrix(
  shared_deg_mat,
  "09_cross_cohort_shared_DEG_number_matrix.tsv"
)

write_matrix(
  direction_concordance_mat,
  "09_cross_cohort_shared_DEG_direction_concordance_matrix.tsv"
)

write.table(
  long_df,
  file = file.path(
    outdir,
    "09_cross_cohort_pairwise_concordance_long.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 12. Human-readable heatmap labels
# ============================================================

spearman_plot <- spearman_mat

rownames(
  spearman_plot
) <- unname(
  labels_2026[
    rownames(spearman_plot)
  ]
)

colnames(
  spearman_plot
) <- unname(
  labels_2021[
    colnames(spearman_plot)
  ]
)

pearson_plot <- pearson_mat

rownames(
  pearson_plot
) <- unname(
  labels_2026[
    rownames(pearson_plot)
  ]
)

colnames(
  pearson_plot
) <- unname(
  labels_2021[
    colnames(pearson_plot)
  ]
)

direction_plot <- direction_concordance_mat

rownames(
  direction_plot
) <- unname(
  labels_2026[
    rownames(direction_plot)
  ]
)

colnames(
  direction_plot
) <- unname(
  labels_2021[
    colnames(direction_plot)
  ]
)


# ============================================================
# 13. Spearman heatmap
# ============================================================

cat("Generating Spearman concordance heatmap...\n")

pdf(
  file.path(
    outdir,
    "09_cross_cohort_Spearman_response_concordance_heatmap.pdf"
  ),
  width = 8,
  height = 8
)

pheatmap(
  spearman_plot,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 9,
  main = paste0(
    "Genome-wide Baseline-relative response concordance\n",
    "2026 primary cohort vs 2021 historical reference"
  ),
  angle_col = 0,
  border_color = "white"
)

dev.off()


# ============================================================
# 14. TIFF Spearman heatmap
# ============================================================

tifile <- file.path(
  outdir,
  "09_cross_cohort_Spearman_response_concordance_heatmap_600dpi.tiff"
)

if (requireNamespace(
  "ragg",
  quietly = TRUE
)) {

  ragg::agg_tiff(
    filename = tifile,
    width = 8,
    height = 8,
    units = "in",
    res = 600,
    compression = "lzw"
  )

} else {

  grDevices::tiff(
    filename = tifile,
    width = 8,
    height = 8,
    units = "in",
    res = 600,
    type = "cairo",
    compression = "lzw"
  )
}

pheatmap(
  spearman_plot,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 9,
  main = paste0(
    "Genome-wide Baseline-relative response concordance\n",
    "2026 primary cohort vs 2021 historical reference"
  ),
  angle_col = 0,
  border_color = "white"
)

dev.off()


# ============================================================
# 15. Pearson heatmap
# ============================================================

pdf(
  file.path(
    outdir,
    "09_cross_cohort_Pearson_response_concordance_heatmap.pdf"
  ),
  width = 8,
  height = 8
)

pheatmap(
  pearson_plot,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 9,
  main = "Pearson response concordance",
  angle_col = 0,
  border_color = "white"
)

dev.off()


# ============================================================
# 16. Shared-DEG direction-concordance heatmap
# ============================================================

pdf(
  file.path(
    outdir,
    "09_shared_DEG_direction_concordance_heatmap.pdf"
  ),
  width = 8,
  height = 8
)

pheatmap(
  direction_plot,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 9,
  main = paste0(
    "Direction concordance among shared DEGs\n",
    "padj < 0.05 and |log2FC| >= 1 in both cohorts"
  ),
  angle_col = 0,
  border_color = "white"
)

dev.off()


# ============================================================
# 17. Best 2026 match for each 2021 stage
# ============================================================

best_for_2021 <- do.call(
  rbind,
  lapply(
    groups_2021,
    function(g21) {

      vals <- spearman_mat[
        ,
        g21
      ]

      best_idx <- which.max(
        vals
      )

      g26 <- names(vals)[
        best_idx
      ]

      row <- long_df[
        long_df$Cohort2026_group ==
          g26 &
        long_df$Cohort2021_group ==
          g21,
        ,
        drop = FALSE
      ]

      data.frame(
        Cohort2021_group = g21,
        Cohort2021_stage =
          unname(labels_2021[g21]),
        Cohort2021_time_h =
          unname(times_2021[g21]),

        Best_matching_2026_group =
          g26,

        Best_matching_2026_stage =
          unname(labels_2026[g26]),

        Best_matching_2026_time_h =
          unname(times_2026[g26]),

        Time_difference_h =
          abs(
            unname(times_2021[g21]) -
            unname(times_2026[g26])
          ),

        Spearman_rho =
          vals[best_idx],

        Pearson_r =
          row$Pearson_r,

        Shared_DEGs =
          row$Shared_DEGs,

        Shared_DEG_direction_concordance =
          row$Shared_DEG_direction_concordance,

        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(best_for_2021) <- NULL

write.table(
  best_for_2021,
  file = file.path(
    outdir,
    "09_best_2026_match_for_each_2021_stage.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 18. Best 2021 match for each 2026 stage
# ============================================================

best_for_2026 <- do.call(
  rbind,
  lapply(
    groups_2026,
    function(g26) {

      vals <- spearman_mat[
        g26,
        ,
        drop = TRUE
      ]

      best_idx <- which.max(
        vals
      )

      g21 <- names(vals)[
        best_idx
      ]

      row <- long_df[
        long_df$Cohort2026_group ==
          g26 &
        long_df$Cohort2021_group ==
          g21,
        ,
        drop = FALSE
      ]

      data.frame(
        Cohort2026_group = g26,
        Cohort2026_stage =
          unname(labels_2026[g26]),
        Cohort2026_time_h =
          unname(times_2026[g26]),

        Best_matching_2021_group =
          g21,

        Best_matching_2021_stage =
          unname(labels_2021[g21]),

        Best_matching_2021_time_h =
          unname(times_2021[g21]),

        Time_difference_h =
          abs(
            unname(times_2026[g26]) -
            unname(times_2021[g21])
          ),

        Spearman_rho =
          vals[best_idx],

        Pearson_r =
          row$Pearson_r,

        Shared_DEGs =
          row$Shared_DEGs,

        Shared_DEG_direction_concordance =
          row$Shared_DEG_direction_concordance,

        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(best_for_2026) <- NULL

write.table(
  best_for_2026,
  file = file.path(
    outdir,
    "09_best_2021_match_for_each_2026_stage.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 19. Overall concordance summaries
# ============================================================

max_index <- which(
  spearman_mat ==
    max(
      spearman_mat,
      na.rm = TRUE
    ),
  arr.ind = TRUE
)[1, ]

min_index <- which(
  spearman_mat ==
    min(
      spearman_mat,
      na.rm = TRUE
    ),
  arr.ind = TRUE
)[1, ]

max_g26 <- rownames(
  spearman_mat
)[
  max_index[1]
]

max_g21 <- colnames(
  spearman_mat
)[
  max_index[2]
]

min_g26 <- rownames(
  spearman_mat
)[
  min_index[1]
]

min_g21 <- colnames(
  spearman_mat
)[
  min_index[2]
]


# ============================================================
# 20. Descriptive relation with temporal distance
# ============================================================

time_distance_similarity_rho <- cor(
  long_df$Absolute_time_difference_h,
  long_df$Spearman_rho,
  method = "spearman"
)


# ============================================================
# 21. Scatter plot: absolute time gap vs response rho
# Descriptive only because 70 pairwise cells are not independent
# ============================================================

p_time <- ggplot(
  long_df,
  aes(
    x = Absolute_time_difference_h,
    y = Spearman_rho
  )
) +
  geom_point(
    size = 2.4,
    alpha = 0.75
  ) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    linewidth = 0.8
  ) +
  labs(
    x = "Absolute difference in post-injury time (h)",
    y = "Genome-wide Spearman response concordance"
  ) +
  theme_classic(
    base_size = 12
  )

ggsave(
  file.path(
    outdir,
    "09_temporal_gap_vs_response_concordance.pdf"
  ),
  plot = p_time,
  width = 7.5,
  height = 5.5,
  units = "in"
)


# ============================================================
# 22. Summary table
# ============================================================

summary_table <- data.frame(
  Metric = c(
    "2026_genes_before_intersection",
    "2021_genes_before_intersection",
    "Common_filtered_genes",
    "2026_timepoints",
    "2021_timepoints",
    "Cross_cohort_pairwise_comparisons",
    "Mean_Spearman_rho",
    "Median_Spearman_rho",
    "Maximum_Spearman_rho",
    "Maximum_rho_2026_stage",
    "Maximum_rho_2021_stage",
    "Minimum_Spearman_rho",
    "Minimum_rho_2026_stage",
    "Minimum_rho_2021_stage",
    "Mean_Pearson_r",
    "Median_Pearson_r",
    "Spearman_between_time_gap_and_response_similarity"
  ),

  Value = c(
    nrow(lfc26_df),
    nrow(lfc21_df),
    length(common_genes),
    length(groups_2026),
    length(groups_2021),
    nrow(long_df),

    mean(
      long_df$Spearman_rho,
      na.rm = TRUE
    ),

    median(
      long_df$Spearman_rho,
      na.rm = TRUE
    ),

    max(
      long_df$Spearman_rho,
      na.rm = TRUE
    ),

    max_g26,
    max_g21,

    min(
      long_df$Spearman_rho,
      na.rm = TRUE
    ),

    min_g26,
    min_g21,

    mean(
      long_df$Pearson_r,
      na.rm = TRUE
    ),

    median(
      long_df$Pearson_r,
      na.rm = TRUE
    ),

    time_distance_similarity_rho
  ),

  stringsAsFactors = FALSE
)

write.table(
  summary_table,
  file = file.path(
    outdir,
    "09_cross_cohort_concordance_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 23. Save aligned response matrices
# ============================================================

write.table(
  data.frame(
    GeneID = common_genes,
    lfc26,
    check.names = FALSE
  ),
  file = file.path(
    outdir,
    "09_2026_common_gene_log2FC_matrix.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  data.frame(
    GeneID = common_genes,
    lfc21,
    check.names = FALSE
  ),
  file = file.path(
    outdir,
    "09_2021_common_gene_log2FC_matrix.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 24. Method note
# ============================================================

method_note <- c(
  "Cross-cohort temporal-response concordance analysis",
  "",
  "The 2026 cohort and 2021 cohort were not pooled.",
  "No batch correction was performed between cohorts.",
  "",
  paste0(
    "Common independently filtered genes: ",
    length(common_genes)
  ),
  "",
  "For each cohort, each injury time point was compared with",
  "that cohort's own uninjured Baseline using DESeq2.",
  "",
  "Primary cross-cohort similarity metric:",
  "Spearman correlation between genome-wide Baseline-relative",
  "log2FoldChange vectors across common genes.",
  "",
  "Pearson correlation was retained as a secondary metric.",
  "",
  "Shared DEG definition:",
  paste0(
    "padj < ",
    PADJ_THRESHOLD,
    " and |log2FoldChange| >= ",
    LFC_THRESHOLD,
    " in both compared time points."
  ),
  "",
  "Shared-DEG direction concordance was calculated as the",
  "fraction of jointly significant genes with the same",
  "log2FoldChange sign in the two cohorts.",
  "",
  "Best-matching time points are descriptive trajectory",
  "annotations and are not interpreted as evidence that the",
  "two injury paradigms are experimentally identical.",
  "",
  "The relation between absolute time difference and response",
  "similarity is descriptive because pairwise matrix cells",
  "are statistically non-independent."
)

writeLines(
  method_note,
  con = file.path(
    outdir,
    "09_cross_cohort_method_definition.txt"
  )
)


# ============================================================
# 25. Session info
# ============================================================

sink(
  file.path(
    outdir,
    "09_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 26. Final console report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("09 CROSS-COHORT CONCORDANCE COMPLETED\n")
cat("============================================================\n")

cat(
  "Common genes                         : ",
  length(common_genes),
  "\n",
  sep = ""
)

cat(
  "2026 stages                          : ",
  length(groups_2026),
  "\n",
  sep = ""
)

cat(
  "2021 stages                          : ",
  length(groups_2021),
  "\n",
  sep = ""
)

cat(
  "Pairwise stage comparisons           : ",
  nrow(long_df),
  "\n",
  sep = ""
)

cat("------------------------------------------------------------\n")

cat(
  "Mean genome-wide Spearman rho        : ",
  sprintf(
    "%.4f",
    mean(
      long_df$Spearman_rho,
      na.rm = TRUE
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Median genome-wide Spearman rho      : ",
  sprintf(
    "%.4f",
    median(
      long_df$Spearman_rho,
      na.rm = TRUE
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Maximum genome-wide Spearman rho     : ",
  sprintf(
    "%.4f",
    max(
      long_df$Spearman_rho,
      na.rm = TRUE
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Maximum-rho stage pair               : 2026 ",
  unname(labels_2026[max_g26]),
  " <-> 2021 ",
  unname(labels_2021[max_g21]),
  "\n",
  sep = ""
)

cat(
  "Minimum genome-wide Spearman rho     : ",
  sprintf(
    "%.4f",
    min(
      long_df$Spearman_rho,
      na.rm = TRUE
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Time-gap vs response-rho Spearman    : ",
  sprintf(
    "%.4f",
    time_distance_similarity_rho
  ),
  "\n",
  sep = ""
)

cat("------------------------------------------------------------\n")

cat(
  "No cohort pooling was performed.\n"
)

cat(
  "No cross-cohort batch correction was performed.\n"
)

cat(
  "Each cohort was normalized and contrasted against its own Baseline.\n"
)

cat(
  "Genome-wide Spearman concordance is the primary validation metric.\n"
)

cat(
  "Pearson and shared-DEG direction concordance are secondary metrics.\n"
)

cat("============================================================\n\n")


cat(
  "Best 2026 match for each 2021 historical-reference stage:\n"
)

print(
  best_for_2021[
    ,
    c(
      "Cohort2021_stage",
      "Best_matching_2026_stage",
      "Time_difference_h",
      "Spearman_rho",
      "Shared_DEGs",
      "Shared_DEG_direction_concordance"
    )
  ],
  row.names = FALSE
)

cat("\nOutput directory:\n")
cat(outdir, "\n\n")

cat("PASS\n")

cat(
  "Step 09 is complete. ",
  "Review the cross-cohort concordance structure before ",
  "proceeding to sensitivity analyses and final figure assembly.\n"
)
