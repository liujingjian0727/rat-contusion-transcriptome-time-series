#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, bitmapType = "cairo")

suppressPackageStartupMessages({

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("R package 'ggplot2' is required.")
  }

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("R package 'pheatmap' is required.")
  }

  library(ggplot2)
  library(pheatmap)
})

cat("\n")
cat("============================================================\n")
cat("03 2026 INTERNAL EXPRESSION STRUCTURE\n")
cat("PCA, sample correlation, hierarchical clustering,\n")
cat("and quantitative outlier diagnostics\n")
cat("Scientific Data rat skeletal muscle transcriptome project\n")
cat("============================================================\n\n")


# ============================================================
# 0. Fixed parameters
# ============================================================

EXPECTED_SAMPLES <- 77
EXPECTED_GENES <- 18364

ROBUST_Z_THRESHOLD <- 3.5
N_PCS_OUTLIER <- 10

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

group_labels <- c(
  "Baseline",
  "0 h",
  "1 h",
  "2 h",
  "6 h",
  "10 h",
  "14 h",
  "18 h",
  "36 h",
  "60 h",
  "72 h"
)

names(group_labels) <- group_levels


# ============================================================
# 1. Input/output
# ============================================================

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

outdir <- "03_2026_internal_expression_structure_result"

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 2. Check input files
# ============================================================

required_files <- c(
  vst_file,
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
# 3. Read VST matrix
# ============================================================

cat("Reading VST matrix...\n")

vst_df <- read.delim(
  vst_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

colnames(vst_df)[1] <- "GeneID"

if (anyDuplicated(vst_df$GeneID)) {
  stop("Duplicated GeneIDs detected in VST matrix.")
}

vst_mat <- as.matrix(
  vst_df[, -1, drop = FALSE]
)

storage.mode(vst_mat) <- "numeric"

rownames(vst_mat) <- vst_df$GeneID

if (anyNA(vst_mat)) {
  stop("NA detected in VST matrix.")
}

if (any(!is.finite(vst_mat))) {
  stop("Non-finite values detected in VST matrix.")
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
    paste(missing_meta_cols, collapse = ", ")
  )
}

if (anyDuplicated(meta$sampleID)) {
  stop("Duplicated sampleIDs detected in metadata.")
}


# ============================================================
# 5. Explicit sample matching
# ============================================================

vst_samples <- colnames(vst_mat)

if (!setequal(vst_samples, meta$sampleID)) {
  stop("VST and metadata sample sets are different.")
}

meta <- meta[
  match(vst_samples, meta$sampleID),
  ,
  drop = FALSE
]

if (!identical(meta$sampleID, vst_samples)) {
  stop("Failed to align metadata to VST sample order.")
}

rownames(meta) <- meta$sampleID

meta$group <- factor(
  meta$group,
  levels = group_levels
)

meta$condition <- factor(
  meta$condition,
  levels = c("Control", "Injury")
)

if (anyNA(meta$group)) {
  stop("Unexpected group detected.")
}

meta$group_label <- factor(
  group_labels[as.character(meta$group)],
  levels = unname(group_labels)
)


# ============================================================
# 6. Dimension validation
# ============================================================

cat("\nInput dimensions:\n")

cat(
  "VST matrix :",
  nrow(vst_mat),
  "genes x",
  ncol(vst_mat),
  "samples\n"
)

cat(
  "Metadata   :",
  nrow(meta),
  "samples\n\n"
)

if (ncol(vst_mat) != EXPECTED_SAMPLES) {
  stop(
    "Unexpected sample number. Expected ",
    EXPECTED_SAMPLES,
    ", observed ",
    ncol(vst_mat),
    "."
  )
}

if (nrow(vst_mat) != EXPECTED_GENES) {
  stop(
    "Unexpected gene number. Expected ",
    EXPECTED_GENES,
    ", observed ",
    nrow(vst_mat),
    "."
  )
}

cat("Dimension validation: PASS\n\n")


# ============================================================
# 7. Group sample summary
# ============================================================

group_summary <- as.data.frame(
  table(meta$group)
)

colnames(group_summary) <- c(
  "Group",
  "Sample_number"
)

group_summary <- group_summary[
  group_summary$Sample_number > 0,
  ,
  drop = FALSE
]

group_summary$Time_label <- group_labels[
  as.character(group_summary$Group)
]

write.table(
  group_summary,
  file = file.path(
    outdir,
    "03_2026_group_sample_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("2026 sample distribution:\n")
print(group_summary, row.names = FALSE)
cat("\n")


# ============================================================
# 8. PCA
# ============================================================

cat("Running PCA on VST-transformed expression matrix...\n")

pca <- prcomp(
  t(vst_mat),
  center = TRUE,
  scale. = FALSE
)

variance <- pca$sdev^2

variance_percent <- 100 * variance / sum(variance)

pca_variance <- data.frame(
  PC = paste0(
    "PC",
    seq_along(variance_percent)
  ),
  Variance_percent = variance_percent,
  Cumulative_variance_percent = cumsum(
    variance_percent
  ),
  stringsAsFactors = FALSE
)

write.table(
  pca_variance,
  file = file.path(
    outdir,
    "03_2026_PCA_variance.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

pca_coordinates <- data.frame(
  sampleID = rownames(pca$x),
  pca$x,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pca_coordinates <- merge(
  meta[, c(
    "sampleID",
    "group",
    "group_label",
    "time",
    "condition"
  )],
  pca_coordinates,
  by = "sampleID",
  sort = FALSE
)

pca_coordinates <- pca_coordinates[
  match(vst_samples, pca_coordinates$sampleID),
  ,
  drop = FALSE
]

write.table(
  pca_coordinates,
  file = file.path(
    outdir,
    "03_2026_PCA_coordinates.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat(
  "PCA variance:\n",
  "PC1 = ",
  sprintf("%.2f", variance_percent[1]),
  "%\n",
  "PC2 = ",
  sprintf("%.2f", variance_percent[2]),
  "%\n\n",
  sep = ""
)


# ============================================================
# 9. PCA plot
# ============================================================

group_palette <- setNames(
  grDevices::hcl.colors(
    length(group_levels),
    palette = "Dynamic"
  ),
  unname(group_labels)
)

p_pca <- ggplot(
  pca_coordinates,
  aes(
    x = PC1,
    y = PC2,
    colour = group_label
  )
) +
  geom_point(
    size = 3.2,
    alpha = 0.88
  ) +
  labs(
    x = paste0(
      "PC1 (",
      sprintf("%.2f", variance_percent[1]),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      sprintf("%.2f", variance_percent[2]),
      "%)"
    ),
    colour = "Sampling stage"
  ) +
  scale_colour_manual(
    values = group_palette,
    drop = FALSE
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(
      face = "bold"
    ),
    axis.title = element_text(
      face = "bold"
    )
  )

ggsave(
  filename = file.path(
    outdir,
    "03_2026_PCA_PC1_PC2.pdf"
  ),
  plot = p_pca,
  width = 8.5,
  height = 6.5,
  units = "in"
)

ggsave(
  filename = file.path(
    outdir,
    "03_2026_PCA_PC1_PC2_600dpi.tiff"
  ),
  plot = p_pca,
  width = 8.5,
  height = 6.5,
  units = "in",
  dpi = 600,
  type = "cairo",
  compression = "lzw"
)


# ============================================================
# 10. PCA plot with sample labels
# ============================================================

p_pca_label <- p_pca +
  geom_text(
    aes(label = sampleID),
    size = 2.2,
    vjust = -0.7,
    check_overlap = TRUE,
    show.legend = FALSE
  )

ggsave(
  filename = file.path(
    outdir,
    "03_2026_PCA_PC1_PC2_with_sample_labels.pdf"
  ),
  plot = p_pca_label,
  width = 10,
  height = 8,
  units = "in"
)


# ============================================================
# 11. PCA scree plot
# ============================================================

n_scree <- min(
  20,
  nrow(pca_variance)
)

scree_data <- pca_variance[
  seq_len(n_scree),
  ,
  drop = FALSE
]

scree_data$PC_number <- seq_len(n_scree)

p_scree <- ggplot(
  scree_data,
  aes(
    x = PC_number,
    y = Variance_percent
  )
) +
  geom_col(
    width = 0.75
  ) +
  geom_point(
    size = 2
  ) +
  scale_x_continuous(
    breaks = seq_len(n_scree)
  ) +
  labs(
    x = "Principal component",
    y = "Explained variance (%)"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    axis.title = element_text(
      face = "bold"
    )
  )

ggsave(
  filename = file.path(
    outdir,
    "03_2026_PCA_scree_first20PCs.pdf"
  ),
  plot = p_scree,
  width = 8,
  height = 5.5,
  units = "in"
)


# ============================================================
# 12. Spearman sample correlation
# ============================================================

cat("Calculating sample-to-sample Spearman correlations...\n")

cor_mat <- cor(
  vst_mat,
  method = "spearman",
  use = "pairwise.complete.obs"
)

diag(cor_mat) <- 1

if (anyNA(cor_mat)) {
  stop("NA detected in sample correlation matrix.")
}

write_matrix <- data.frame(
  sampleID = rownames(cor_mat),
  cor_mat,
  check.names = FALSE
)

write.table(
  write_matrix,
  file = file.path(
    outdir,
    "03_2026_sample_Spearman_correlation_matrix.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 13. Hierarchical clustering based on 1-Spearman
# ============================================================

# Convert Spearman correlation to a proper square distance matrix.
# Do not use pmax(0, 1 - cor_mat) here because pmax() may drop
# the matrix dimensions when the scalar is supplied first.
cor_distance_mat <- 1 - cor_mat

# Protect against tiny negative values caused by floating-point precision
cor_distance_mat[cor_distance_mat < 0] <- 0

# A sample must have zero distance to itself
diag(cor_distance_mat) <- 0

# Explicit validation before conversion to a dist object
if (!is.matrix(cor_distance_mat)) {
  stop("Correlation-derived distance object is not a matrix.")
}

if (nrow(cor_distance_mat) != ncol(cor_distance_mat)) {
  stop(
    "Correlation-derived distance matrix is not square: ",
    nrow(cor_distance_mat),
    " x ",
    ncol(cor_distance_mat)
  )
}

if (!identical(
  rownames(cor_distance_mat),
  colnames(cor_distance_mat)
)) {
  stop("Row and column sample orders differ in distance matrix.")
}

if (anyNA(cor_distance_mat) ||
    any(!is.finite(cor_distance_mat))) {
  stop("Invalid values detected in correlation-derived distance matrix.")
}

cor_distance <- as.dist(
  cor_distance_mat,
  diag = FALSE,
  upper = FALSE
)

expected_dist_length <- ncol(cor_mat) * (ncol(cor_mat) - 1) / 2

if (length(cor_distance) != expected_dist_length) {
  stop(
    "Invalid distance-vector length. Expected ",
    expected_dist_length,
    ", observed ",
    length(cor_distance),
    "."
  )
}

cat(
  "Correlation distance matrix validation: ",
  nrow(cor_distance_mat),
  " x ",
  ncol(cor_distance_mat),
  " square matrix: PASS\n",
  sep = ""
)

cat(
  "Distance-vector length: ",
  length(cor_distance),
  " (expected ",
  expected_dist_length,
  "): PASS\n\n",
  sep = ""
)

sample_hclust <- hclust(
  cor_distance,
  method = "average"
)


# ============================================================
# 14. Heatmap annotation
# ============================================================

annotation <- data.frame(
  Stage = meta$group_label,
  Condition = meta$condition,
  row.names = meta$sampleID,
  stringsAsFactors = FALSE
)

stage_colors <- setNames(
  group_palette,
  unname(group_labels)
)

condition_colors <- setNames(
  grDevices::hcl.colors(
    2,
    palette = "Dark 3"
  ),
  c("Control", "Injury")
)

annotation_colors <- list(
  Stage = stage_colors,
  Condition = condition_colors
)


# ============================================================
# 15. Spearman correlation heatmap
# ============================================================

cat("Generating Spearman correlation heatmap...\n")

pdf(
  file.path(
    outdir,
    "03_2026_sample_Spearman_correlation_heatmap.pdf"
  ),
  width = 14,
  height = 13
)

pheatmap(
  cor_mat,
  cluster_rows = sample_hclust,
  cluster_cols = sample_hclust,
  annotation_row = annotation,
  annotation_col = annotation,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 5,
  fontsize_col = 5,
  border_color = NA,
  main = "2026 sample-to-sample Spearman correlation"
)

dev.off()


tiff(
  file.path(
    outdir,
    "03_2026_sample_Spearman_correlation_heatmap_600dpi.tiff"
  ),
  width = 14,
  height = 13,
  units = "in",
  res = 600,
  type = "cairo",
  compression = "lzw"
)

pheatmap(
  cor_mat,
  cluster_rows = sample_hclust,
  cluster_cols = sample_hclust,
  annotation_row = annotation,
  annotation_col = annotation,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 5,
  fontsize_col = 5,
  border_color = NA,
  main = "2026 sample-to-sample Spearman correlation"
)

dev.off()


# ============================================================
# 16. Standalone hierarchical clustering dendrogram
# ============================================================

cat("Generating hierarchical clustering dendrogram...\n")

pdf(
  file.path(
    outdir,
    "03_2026_sample_hierarchical_clustering.pdf"
  ),
  width = 14,
  height = 7
)

par(
  mar = c(8, 4, 3, 1)
)

plot(
  sample_hclust,
  main = "Hierarchical clustering of 2026 samples",
  xlab = "",
  sub = "",
  ylab = "1 - Spearman correlation",
  cex = 0.55,
  hang = -1
)

dev.off()


# ============================================================
# 17. Quantitative outlier diagnostics
# ============================================================

cat("Calculating quantitative expression-level outlier metrics...\n")

sample_ids <- meta$sampleID

# ------------------------------------------------------------
# 17A. Global median correlation
# ------------------------------------------------------------

global_median_corr <- sapply(
  sample_ids,
  function(s) {

    x <- cor_mat[
      s,
      setdiff(sample_ids, s)
    ]

    median(
      x,
      na.rm = TRUE
    )
  }
)


# ------------------------------------------------------------
# 17B. Within-group median correlation
# ------------------------------------------------------------

within_group_median_corr <- numeric(
  length(sample_ids)
)

names(within_group_median_corr) <- sample_ids


# ------------------------------------------------------------
# 17C. Leave-one-out full-expression centroid RMS distance
# ------------------------------------------------------------

loo_centroid_rms_distance <- numeric(
  length(sample_ids)
)

names(loo_centroid_rms_distance) <- sample_ids


# ------------------------------------------------------------
# 17D. Leave-one-out PCA-space centroid RMS distance
# ------------------------------------------------------------

n_pc_use <- min(
  N_PCS_OUTLIER,
  ncol(pca$x)
)

pca_for_distance <- pca$x[
  ,
  seq_len(n_pc_use),
  drop = FALSE
]

loo_pca_rms_distance <- numeric(
  length(sample_ids)
)

names(loo_pca_rms_distance) <- sample_ids


# ------------------------------------------------------------
# Calculate within each biological group
# ------------------------------------------------------------

for (grp in group_levels) {

  grp_samples <- meta$sampleID[
    meta$group == grp
  ]

  if (length(grp_samples) < 2) {
    stop(
      "Group ",
      grp,
      " has fewer than two samples."
    )
  }

  for (s in grp_samples) {

    other_samples <- setdiff(
      grp_samples,
      s
    )

    # Within-group median Spearman correlation
    within_group_median_corr[s] <- median(
      cor_mat[
        s,
        other_samples
      ],
      na.rm = TRUE
    )

    # Leave-one-out centroid in full VST expression space
    loo_centroid <- rowMeans(
      vst_mat[
        ,
        other_samples,
        drop = FALSE
      ]
    )

    sample_vector <- vst_mat[
      ,
      s
    ]

    loo_centroid_rms_distance[s] <- sqrt(
      mean(
        (sample_vector - loo_centroid)^2
      )
    )

    # Leave-one-out centroid in first PCs
    loo_pca_centroid <- colMeans(
      pca_for_distance[
        other_samples,
        ,
        drop = FALSE
      ]
    )

    sample_pc_vector <- pca_for_distance[
      s,
      ,
      drop = TRUE
    ]

    loo_pca_rms_distance[s] <- sqrt(
      mean(
        (sample_pc_vector -
           loo_pca_centroid)^2
      )
    )
  }
}


# ============================================================
# 18. Robust within-group Z scores
# ============================================================

robust_z <- function(x) {

  med <- median(
    x,
    na.rm = TRUE
  )

  mad_value <- mad(
    x,
    center = med,
    constant = 1.4826,
    na.rm = TRUE
  )

  if (!is.finite(mad_value) ||
      mad_value <= .Machine$double.eps) {

    return(
      rep(
        0,
        length(x)
      )
    )
  }

  (x - med) / mad_value
}


diagnostics <- data.frame(
  sampleID = sample_ids,
  group = as.character(meta$group),
  group_label = as.character(meta$group_label),
  condition = as.character(meta$condition),
  time = meta$time,
  Global_median_Spearman = global_median_corr[
    sample_ids
  ],
  Within_group_median_Spearman =
    within_group_median_corr[
      sample_ids
    ],
  LOO_group_centroid_RMS_VST =
    loo_centroid_rms_distance[
      sample_ids
    ],
  LOO_group_centroid_RMS_PCA10 =
    loo_pca_rms_distance[
      sample_ids
    ],
  stringsAsFactors = FALSE
)

diagnostics$Within_group_correlation_robustZ <- NA_real_
diagnostics$VST_centroid_distance_robustZ <- NA_real_
diagnostics$PCA10_centroid_distance_robustZ <- NA_real_


for (grp in group_levels) {

  idx <- diagnostics$group == grp

  diagnostics$Within_group_correlation_robustZ[idx] <-
    robust_z(
      diagnostics$Within_group_median_Spearman[idx]
    )

  diagnostics$VST_centroid_distance_robustZ[idx] <-
    robust_z(
      diagnostics$LOO_group_centroid_RMS_VST[idx]
    )

  diagnostics$PCA10_centroid_distance_robustZ[idx] <-
    robust_z(
      diagnostics$LOO_group_centroid_RMS_PCA10[idx]
    )
}


# ============================================================
# 19. Predefined diagnostic flags
# ============================================================

diagnostics$Flag_low_within_group_correlation <-
  diagnostics$Within_group_correlation_robustZ <=
  -ROBUST_Z_THRESHOLD

diagnostics$Flag_high_VST_centroid_distance <-
  diagnostics$VST_centroid_distance_robustZ >=
  ROBUST_Z_THRESHOLD

diagnostics$Flag_high_PCA10_centroid_distance <-
  diagnostics$PCA10_centroid_distance_robustZ >=
  ROBUST_Z_THRESHOLD

diagnostics$Number_of_expression_flags <- rowSums(
  diagnostics[
    ,
    c(
      "Flag_low_within_group_correlation",
      "Flag_high_VST_centroid_distance",
      "Flag_high_PCA10_centroid_distance"
    )
  ]
)

diagnostics$Review_level <- ifelse(
  diagnostics$Number_of_expression_flags >= 2,
  "Multi_metric_review",
  ifelse(
    diagnostics$Number_of_expression_flags == 1,
    "Single_metric_review",
    "No_expression_flag"
  )
)


# ============================================================
# 20. Order diagnostic table
# ============================================================

group_order_number <- match(
  diagnostics$group,
  group_levels
)

sample_numeric <- suppressWarnings(
  as.numeric(
    sub(
      "^.*_",
      "",
      diagnostics$sampleID
    )
  )
)

sample_numeric[
  is.na(sample_numeric)
] <- seq_len(
  sum(is.na(sample_numeric))
)

diagnostics <- diagnostics[
  order(
    group_order_number,
    sample_numeric
  ),
  ,
  drop = FALSE
]

write.table(
  diagnostics,
  file = file.path(
    outdir,
    "03_2026_expression_outlier_diagnostics.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 21. Expression flag summary
# ============================================================

flag_summary <- as.data.frame(
  table(
    diagnostics$Review_level
  )
)

colnames(flag_summary) <- c(
  "Review_level",
  "Sample_number"
)

write.table(
  flag_summary,
  file = file.path(
    outdir,
    "03_2026_expression_outlier_flag_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 22. Group-level reproducibility summary
# ============================================================

group_reproducibility <- do.call(
  rbind,
  lapply(
    group_levels,
    function(grp) {

      x <- diagnostics[
        diagnostics$group == grp,
        ,
        drop = FALSE
      ]

      data.frame(
        Group = grp,
        Time_label = group_labels[grp],
        N = nrow(x),

        Median_within_group_Spearman = median(
          x$Within_group_median_Spearman,
          na.rm = TRUE
        ),

        Minimum_within_group_Spearman = min(
          x$Within_group_median_Spearman,
          na.rm = TRUE
        ),

        Maximum_within_group_Spearman = max(
          x$Within_group_median_Spearman,
          na.rm = TRUE
        ),

        Median_LOO_centroid_RMS_VST = median(
          x$LOO_group_centroid_RMS_VST,
          na.rm = TRUE
        ),

        Median_LOO_centroid_RMS_PCA10 = median(
          x$LOO_group_centroid_RMS_PCA10,
          na.rm = TRUE
        ),

        Flagged_samples = sum(
          x$Number_of_expression_flags > 0
        ),

        stringsAsFactors = FALSE
      )
    }
  )
)

write.table(
  group_reproducibility,
  file = file.path(
    outdir,
    "03_2026_group_reproducibility_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 23. Diagnostic scatter plot
# ============================================================

diagnostics$group_label_factor <- factor(
  group_labels[diagnostics$group],
  levels = unname(group_labels)
)

diagnostics$Review_level <- factor(
  diagnostics$Review_level,
  levels = c(
    "No_expression_flag",
    "Single_metric_review",
    "Multi_metric_review"
  )
)

p_diag <- ggplot(
  diagnostics,
  aes(
    x = Within_group_median_Spearman,
    y = LOO_group_centroid_RMS_VST,
    colour = group_label_factor,
    shape = Review_level
  )
) +
  geom_point(
    size = 3.1,
    alpha = 0.9
  ) +
  geom_text(
    data = diagnostics[
      diagnostics$Number_of_expression_flags > 0,
      ,
      drop = FALSE
    ],
    aes(label = sampleID),
    size = 2.6,
    vjust = -0.8,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = group_palette,
    drop = FALSE
  ) +
  labs(
    x = "Within-group median Spearman correlation",
    y = "Leave-one-out group-centroid RMS distance (VST)",
    colour = "Sampling stage",
    shape = "Expression QC"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    axis.title = element_text(
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    )
  )

ggsave(
  filename = file.path(
    outdir,
    "03_2026_expression_outlier_diagnostic_plot.pdf"
  ),
  plot = p_diag,
  width = 9,
  height = 6.8,
  units = "in"
)

ggsave(
  filename = file.path(
    outdir,
    "03_2026_expression_outlier_diagnostic_plot_600dpi.tiff"
  ),
  plot = p_diag,
  width = 9,
  height = 6.8,
  units = "in",
  dpi = 600,
  type = "cairo",
  compression = "lzw"
)


# ============================================================
# 24. PCA values added to outlier diagnostic file
# ============================================================

pc_add <- pca_coordinates[
  ,
  c(
    "sampleID",
    "PC1",
    "PC2"
  )
]

diagnostics_with_pca <- merge(
  diagnostics,
  pc_add,
  by = "sampleID",
  sort = FALSE
)

diagnostics_with_pca <- diagnostics_with_pca[
  match(
    diagnostics$sampleID,
    diagnostics_with_pca$sampleID
  ),
  ,
  drop = FALSE
]

write.table(
  diagnostics_with_pca,
  file = file.path(
    outdir,
    "03_2026_expression_outlier_diagnostics_with_PCA.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 25. Identify samples requiring review
# ============================================================

review_samples <- diagnostics[
  diagnostics$Number_of_expression_flags > 0,
  c(
    "sampleID",
    "group",
    "Within_group_median_Spearman",
    "LOO_group_centroid_RMS_VST",
    "LOO_group_centroid_RMS_PCA10",
    "Within_group_correlation_robustZ",
    "VST_centroid_distance_robustZ",
    "PCA10_centroid_distance_robustZ",
    "Number_of_expression_flags",
    "Review_level"
  ),
  drop = FALSE
]

write.table(
  review_samples,
  file = file.path(
    outdir,
    "03_2026_samples_requiring_expression_review.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 26. Save PCA object and clustering object
# ============================================================

saveRDS(
  pca,
  file = file.path(
    outdir,
    "03_2026_PCA_object.rds"
  )
)

saveRDS(
  sample_hclust,
  file = file.path(
    outdir,
    "03_2026_sample_hclust_object.rds"
  )
)


# ============================================================
# 27. Save R session information
# ============================================================

sink(
  file.path(
    outdir,
    "03_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 28. Final report
# ============================================================

n_single <- sum(
  diagnostics$Review_level ==
    "Single_metric_review"
)

n_multi <- sum(
  diagnostics$Review_level ==
    "Multi_metric_review"
)

n_no_flag <- sum(
  diagnostics$Review_level ==
    "No_expression_flag"
)

cat("\n")
cat("============================================================\n")
cat("03 2026 INTERNAL EXPRESSION STRUCTURE COMPLETED\n")
cat("============================================================\n")

cat(
  "Samples analyzed                  :",
  ncol(vst_mat),
  "\n"
)

cat(
  "Genes analyzed                    :",
  nrow(vst_mat),
  "\n"
)

cat(
  "PC1 variance                      :",
  sprintf("%.2f%%", variance_percent[1]),
  "\n"
)

cat(
  "PC2 variance                      :",
  sprintf("%.2f%%", variance_percent[2]),
  "\n"
)

cat(
  "Expression samples without flag   :",
  n_no_flag,
  "\n"
)

cat(
  "Single-metric review samples      :",
  n_single,
  "\n"
)

cat(
  "Multi-metric review samples       :",
  n_multi,
  "\n"
)

cat(
  "Robust-Z review threshold         : +/-",
  ROBUST_Z_THRESHOLD,
  "\n"
)

cat("------------------------------------------------------------\n")

if (nrow(review_samples) == 0) {

  cat("Samples requiring expression review: None\n")

} else {

  cat("Samples requiring expression review:\n")

  for (i in seq_len(nrow(review_samples))) {

    cat(
      "  ",
      review_samples$sampleID[i],
      "  [",
      review_samples$group[i],
      "]  ",
      review_samples$Review_level[i],
      "\n",
      sep = ""
    )
  }
}

cat("------------------------------------------------------------\n")
cat("IMPORTANT:\n")
cat("No sample was removed in this step.\n")
cat("Expression-level flags are diagnostic review flags only.\n")
cat("Sample exclusion requires integration with previous raw/QC evidence.\n")
cat("------------------------------------------------------------\n")
cat("PCA                              : PASS\n")
cat("Spearman correlation             : PASS\n")
cat("Hierarchical clustering          : PASS\n")
cat("Quantitative outlier diagnostics : PASS\n")
cat("============================================================\n\n")

cat("Output directory:\n")
cat(outdir, "\n\n")

cat("PASS\n")
cat(
  "Step 03 is complete. Review PCA, correlation structure, ",
  "and flagged samples before proceeding to Step 04 ",
  "(temporal structure: PERMANOVA, PERMDISP and centroid trajectory).\n"
)
