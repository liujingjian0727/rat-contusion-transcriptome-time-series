#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  bitmapType = "cairo"
)

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(vegan)
})

cat("\n")
cat("============================================================\n")
cat("12C COVERAGE-ADJUSTED EXPRESSION-SPACE SENSITIVITY\n")
cat("Final 77-sample 2026 primary cohort\n")
cat("============================================================\n\n")


# ============================================================
# 0. Inputs / outputs
# ============================================================

metadata_file <-
  "rattus_meta_2026_77samples.tsv"

vst_file <- paste0(
  "02_count_filtering_and_VST_result/",
  "2026_primary/",
  "02_2026_VST_matrix.tsv"
)

coverage_file <- paste0(
  "supplementary_figures/",
  "Supplementary_Figure_S1_gene_body_coverage/",
  "S1_per_sample_gene_body_coverage_metrics.tsv"
)

outdir <- paste0(
  "supplementary_figures/",
  "Supplementary_Figure_S1_gene_body_coverage/",
  "12C_coverage_adjusted_PCA_sensitivity"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. Fixed stage definitions
# ============================================================

groups <- c(
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

stage_labels <- c(
  Baseline = "Baseline",
  R0h  = "0 h",
  R1h  = "1 h",
  R2h  = "2 h",
  R6h  = "6 h",
  R10h = "10 h",
  R14h = "14 h",
  R18h = "18 h",
  R36h = "36 h",
  R60h = "60 h",
  R72h = "72 h"
)

stage_colors <- c(
  Baseline = "#7A7A7A",
  R0h  = "#4C78A8",
  R1h  = "#5A88B5",
  R2h  = "#6998C1",
  R6h  = "#72A5B5",
  R10h = "#79B39C",
  R14h = "#8DBB78",
  R18h = "#B4B85D",
  R36h = "#D1A354",
  R60h = "#D68155",
  R72h = "#C85A5A"
)


# ============================================================
# 2. Validate inputs
# ============================================================

required_files <- c(
  metadata_file,
  vst_file,
  coverage_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing required file(s):\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

cat("Required input files found: PASS\n")


# ============================================================
# 3. Read metadata
# ============================================================

meta <- read.delim(
  metadata_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

if (!all(
  c(
    "sampleID",
    "group"
  ) %in%
    colnames(meta)
)) {
  stop(
    "Metadata must contain sampleID and group."
  )
}

if (nrow(meta) != 77) {
  stop(
    "Expected 77 metadata rows; observed ",
    nrow(meta)
  )
}

if (anyDuplicated(meta$sampleID)) {
  stop(
    "Duplicated sampleID in metadata."
  )
}

meta$group <- factor(
  meta$group,
  levels = groups
)

final_samples <- meta$sampleID


# ============================================================
# 4. Read frozen VST matrix
# ============================================================

vst_df <- read.delim(
  vst_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

colnames(vst_df)[1] <- "GeneID"

if (anyDuplicated(vst_df$GeneID)) {
  stop(
    "Duplicated GeneID in VST matrix."
  )
}

vst_mat <- as.matrix(
  vst_df[
    ,
    -1,
    drop = FALSE
  ]
)

storage.mode(vst_mat) <- "numeric"

rownames(vst_mat) <- vst_df$GeneID

if (!setequal(
  colnames(vst_mat),
  final_samples
)) {
  stop(
    "VST sample set differs from final metadata."
  )
}

vst_mat <- vst_mat[
  ,
  final_samples,
  drop = FALSE
]

if (anyNA(vst_mat) ||
    any(!is.finite(vst_mat))) {
  stop(
    "Invalid values in VST matrix."
  )
}

cat(
  "VST matrix : ",
  nrow(vst_mat),
  " genes x ",
  ncol(vst_mat),
  " samples\n",
  sep = ""
)


# ============================================================
# 5. Read coverage metrics
# ============================================================

coverage <- read.delim(
  coverage_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

required_cov_cols <- c(
  "sampleID",
  "Ratio_3prime_to_5prime",
  "Spearman_to_overall_median"
)

if (!all(
  required_cov_cols %in%
    colnames(coverage)
)) {
  stop(
    "Coverage table lacks required columns."
  )
}

if (anyDuplicated(coverage$sampleID)) {
  stop(
    "Duplicated coverage sample IDs."
  )
}

if (!setequal(
  coverage$sampleID,
  final_samples
)) {
  stop(
    "Coverage sample set differs from final metadata."
  )
}

coverage <- coverage[
  match(
    final_samples,
    coverage$sampleID
  ),
  ,
  drop = FALSE
]


# ============================================================
# 6. Build model-data table
# ============================================================

model_data <- data.frame(
  sampleID = final_samples,
  group = meta$group,

  ratio =
    coverage$Ratio_3prime_to_5prime,

  profile_similarity =
    coverage$Spearman_to_overall_median,

  stringsAsFactors = FALSE
)

if (anyNA(
  model_data[
    ,
    c(
      "ratio",
      "profile_similarity"
    )
  ]
)) {
  stop(
    "NA detected in coverage metrics."
  )
}

model_data$z_ratio <- as.numeric(
  scale(
    model_data$ratio
  )
)

model_data$z_profile <- as.numeric(
  scale(
    model_data$profile_similarity
  )
)

rownames(model_data) <-
  model_data$sampleID


# ============================================================
# 7. Final cohort sanity checks
# ============================================================

cat(
  "Aligned samples : ",
  nrow(model_data),
  "\n",
  sep = ""
)

cat(
  "R2h_1 present   : ",
  "R2h_1" %in%
    model_data$sampleID,
  "\n",
  sep = ""
)

cat(
  "R2h_8 present   : ",
  "R2h_8" %in%
    model_data$sampleID,
  "\n\n",
  sep = ""
)

if ("R2h_1" %in%
    model_data$sampleID) {
  stop(
    "R2h_1 unexpectedly present."
  )
}

if (!("R2h_8" %in%
      model_data$sampleID)) {
  stop(
    "R2h_8 unexpectedly absent."
  )
}


# ============================================================
# 8. Matrix regression helper
#
# Expression input:
# genes x samples
#
# Regression is performed as:
# samples x genes
# ============================================================

remove_selected_effects <- function(
    expression_matrix,
    design_matrix,
    effect_columns) {

  Y <- t(
    expression_matrix
  )

  if (nrow(Y) !=
      nrow(design_matrix)) {
    stop(
      "Design/expression sample mismatch."
    )
  }

  design_rank <- qr(
    design_matrix
  )$rank

  if (design_rank !=
      ncol(design_matrix)) {
    stop(
      "Design matrix is not full rank."
    )
  }

  beta <- qr.solve(
    design_matrix,
    Y
  )

  cov_effect <- design_matrix[
    ,
    effect_columns,
    drop = FALSE
  ] %*%
    beta[
      effect_columns,
      ,
      drop = FALSE
    ]

  Y_adjusted <- Y -
    cov_effect

  adjusted <- t(
    Y_adjusted
  )

  rownames(adjusted) <-
    rownames(expression_matrix)

  colnames(adjusted) <-
    colnames(expression_matrix)

  list(
    adjusted = adjusted,
    beta = beta,
    design_rank = design_rank
  )
}


# ============================================================
# 9. Adjustment strategy 1
#
# CONDITIONAL / STAGE-PRESERVING
#
# expression ~ stage + z_ratio + z_profile
#
# Only the conditional coverage contribution is removed.
# ============================================================

design_conditional <- model.matrix(
  ~ group + z_ratio + z_profile,
  data = model_data
)

cov_columns_conditional <- which(
  colnames(design_conditional) %in%
    c(
      "z_ratio",
      "z_profile"
    )
)

if (length(
  cov_columns_conditional
) != 2) {
  stop(
    "Could not identify coverage columns ",
    "in conditional design."
  )
}

conditional_result <-
  remove_selected_effects(
    expression_matrix =
      vst_mat,

    design_matrix =
      design_conditional,

    effect_columns =
      cov_columns_conditional
  )

vst_conditional <-
  conditional_result$adjusted


# ============================================================
# 10. Adjustment strategy 2
#
# AGGRESSIVE / COVERAGE-ONLY
#
# expression ~ z_ratio + z_profile
#
# All linear expression signal associated with coverage
# metrics is removed, irrespective of stage.
#
# This may also remove true stage-related biology because
# coverage metrics and stage are correlated.
# ============================================================

design_aggressive <- model.matrix(
  ~ z_ratio + z_profile,
  data = model_data
)

cov_columns_aggressive <- which(
  colnames(design_aggressive) %in%
    c(
      "z_ratio",
      "z_profile"
    )
)

if (length(
  cov_columns_aggressive
) != 2) {
  stop(
    "Could not identify coverage columns ",
    "in aggressive design."
  )
}

aggressive_result <-
  remove_selected_effects(
    expression_matrix =
      vst_mat,

    design_matrix =
      design_aggressive,

    effect_columns =
      cov_columns_aggressive
  )

vst_aggressive <-
  aggressive_result$adjusted


# ============================================================
# 11. Sanity check dimensions
# ============================================================

stopifnot(
  identical(
    dim(vst_mat),
    dim(vst_conditional)
  )
)

stopifnot(
  identical(
    dim(vst_mat),
    dim(vst_aggressive)
  )
)

cat(
  "Conditional design rank : ",
  conditional_result$design_rank,
  "/",
  ncol(design_conditional),
  "\n",
  sep = ""
)

cat(
  "Aggressive design rank  : ",
  aggressive_result$design_rank,
  "/",
  ncol(design_aggressive),
  "\n\n",
  sep = ""
)


# ============================================================
# 12. PCA helper
# ============================================================

run_pca <- function(
    expression_matrix,
    version_name) {

  p <- prcomp(
    t(expression_matrix),
    center = TRUE,
    scale. = FALSE
  )

  variance <- 100 *
    p$sdev^2 /
    sum(
      p$sdev^2
    )

  coords <- data.frame(
    sampleID = rownames(
      p$x
    ),
    PC1 = p$x[, 1],
    PC2 = p$x[, 2],
    PC3 = p$x[, 3],
    PC4 = p$x[, 4],
    PC5 = p$x[, 5],
    Version = version_name,
    stringsAsFactors = FALSE
  )

  coords$group <-
    model_data$group[
      match(
        coords$sampleID,
        model_data$sampleID
      )
    ]

  list(
    object = p,
    variance = variance,
    coordinates = coords
  )
}


pca_original <- run_pca(
  vst_mat,
  "Original"
)

pca_conditional <- run_pca(
  vst_conditional,
  "Conditional"
)

pca_aggressive <- run_pca(
  vst_aggressive,
  "Aggressive"
)


cat(
  "Original PCA    : PC1 = ",
  sprintf(
    "%.2f%%",
    pca_original$variance[1]
  ),
  "; PC2 = ",
  sprintf(
    "%.2f%%",
    pca_original$variance[2]
  ),
  "\n",
  sep = ""
)

cat(
  "Conditional PCA : PC1 = ",
  sprintf(
    "%.2f%%",
    pca_conditional$variance[1]
  ),
  "; PC2 = ",
  sprintf(
    "%.2f%%",
    pca_conditional$variance[2]
  ),
  "\n",
  sep = ""
)

cat(
  "Aggressive PCA  : PC1 = ",
  sprintf(
    "%.2f%%",
    pca_aggressive$variance[1]
  ),
  "; PC2 = ",
  sprintf(
    "%.2f%%",
    pca_aggressive$variance[2]
  ),
  "\n\n",
  sep = ""
)


# ============================================================
# 13. Euclidean expression-space distances
# ============================================================

dist_original <- dist(
  t(vst_mat),
  method = "euclidean"
)

dist_conditional <- dist(
  t(vst_conditional),
  method = "euclidean"
)

dist_aggressive <- dist(
  t(vst_aggressive),
  method = "euclidean"
)


upper_original <-
  as.vector(
    dist_original
  )

upper_conditional <-
  as.vector(
    dist_conditional
  )

upper_aggressive <-
  as.vector(
    dist_aggressive
  )


distance_concordance <- data.frame(
  Comparison = c(
    "Original_vs_Conditional",
    "Original_vs_Aggressive"
  ),

  Pearson_r = c(
    cor(
      upper_original,
      upper_conditional,
      method = "pearson"
    ),
    cor(
      upper_original,
      upper_aggressive,
      method = "pearson"
    )
  ),

  Spearman_rho = c(
    cor(
      upper_original,
      upper_conditional,
      method = "spearman"
    ),
    cor(
      upper_original,
      upper_aggressive,
      method = "spearman"
    )
  ),

  Mean_absolute_distance_difference = c(
    mean(
      abs(
        upper_original -
          upper_conditional
      )
    ),
    mean(
      abs(
        upper_original -
          upper_aggressive
      )
    )
  ),

  stringsAsFactors = FALSE
)


# ============================================================
# 14. PERMANOVA / PERMDISP helper
# ============================================================

run_permanova <- function(
    distance_object,
    version_name) {

  # Match frozen Step 04 implementation:
  # 9999 permutations
  # fixed random seed = 20260910
  # spatial median
  # bias adjustment for unequal group sizes

  set.seed(20260910)

  perm <- vegan::adonis2(
    distance_object ~ group,
    data = model_data,
    permutations = 9999,
    by = "margin"
  )

  dispersion <- vegan::betadisper(
    distance_object,
    group = model_data$group,
    type = "median",
    bias.adjust = TRUE
  )

  set.seed(20260910)

  dispersion_perm <- vegan::permutest(
    dispersion,
    permutations = 9999
  )

  data.frame(
    Version = version_name,

    PERMANOVA_R2 =
      perm[
        "group",
        "R2"
      ],

    PERMANOVA_F =
      perm[
        "group",
        "F"
      ],

    PERMANOVA_P =
      perm[
        "group",
        "Pr(>F)"
      ],

    PERMDISP_F =
      dispersion_perm$tab[
        "Groups",
        "F"
      ],

    PERMDISP_P =
      dispersion_perm$tab[
        "Groups",
        "Pr(>F)"
      ],

    stringsAsFactors = FALSE
  )
}


permanova_summary <- rbind(
  run_permanova(
    dist_original,
    "Original"
  ),

  run_permanova(
    dist_conditional,
    "Conditional"
  ),

  run_permanova(
    dist_aggressive,
    "Aggressive"
  )
)


# ============================================================
# 15. Stage centroids in full expression space
# ============================================================

calculate_centroids <- function(
    expression_matrix,
    version_name) {

  result <- matrix(
    NA_real_,
    nrow = length(groups),
    ncol = nrow(expression_matrix),
    dimnames = list(
      groups,
      rownames(expression_matrix)
    )
  )

  for (g in groups) {

    samples_g <- model_data$sampleID[
      model_data$group == g
    ]

    result[
      g,
    ] <- rowMeans(
      expression_matrix[
        ,
        samples_g,
        drop = FALSE
      ]
    )
  }

  attr(
    result,
    "Version"
  ) <- version_name

  result
}


centroids_original <-
  calculate_centroids(
    vst_mat,
    "Original"
  )

centroids_conditional <-
  calculate_centroids(
    vst_conditional,
    "Conditional"
  )

centroids_aggressive <-
  calculate_centroids(
    vst_aggressive,
    "Aggressive"
  )


# ============================================================
# 16. Stage-centroid distance matrices
# ============================================================

centroid_dist_original <- dist(
  centroids_original
)

centroid_dist_conditional <- dist(
  centroids_conditional
)

centroid_dist_aggressive <- dist(
  centroids_aggressive
)


centroid_distance_concordance <- data.frame(
  Comparison = c(
    "Original_vs_Conditional",
    "Original_vs_Aggressive"
  ),

  Pearson_r = c(
    cor(
      as.vector(
        centroid_dist_original
      ),
      as.vector(
        centroid_dist_conditional
      ),
      method = "pearson"
    ),

    cor(
      as.vector(
        centroid_dist_original
      ),
      as.vector(
        centroid_dist_aggressive
      ),
      method = "pearson"
    )
  ),

  Spearman_rho = c(
    cor(
      as.vector(
        centroid_dist_original
      ),
      as.vector(
        centroid_dist_conditional
      ),
      method = "spearman"
    ),

    cor(
      as.vector(
        centroid_dist_original
      ),
      as.vector(
        centroid_dist_aggressive
      ),
      method = "spearman"
    )
  ),

  stringsAsFactors = FALSE
)


# ============================================================
# 17. Adjacent-stage RMS distances
#
# Same general concept used for Figure 3:
# sqrt(mean((centroid_2 - centroid_1)^2))
# ============================================================

calculate_adjacent_rms <- function(
    centroid_matrix,
    version_name) {

  out <- vector(
    "list",
    length(groups) - 1
  )

  for (i in seq_len(
    length(groups) - 1
  )) {

    g1 <- groups[i]
    g2 <- groups[i + 1]

    diff_vec <-
      centroid_matrix[
        g2,
      ] -
      centroid_matrix[
        g1,
      ]

    rms <- sqrt(
      mean(
        diff_vec^2
      )
    )

    out[[i]] <- data.frame(
      Version = version_name,
      From = g1,
      To = g2,
      Transition = paste0(
        unname(
          stage_labels[g1]
        ),
        " → ",
        unname(
          stage_labels[g2]
        )
      ),
      Order = i,
      RMS_distance = rms,
      stringsAsFactors = FALSE
    )
  }

  do.call(
    rbind,
    out
  )
}


adjacent_rms <- rbind(
  calculate_adjacent_rms(
    centroids_original,
    "Original"
  ),

  calculate_adjacent_rms(
    centroids_conditional,
    "Conditional"
  ),

  calculate_adjacent_rms(
    centroids_aggressive,
    "Aggressive"
  )
)


# ============================================================
# 18. Concordance of adjacent-stage RMS profiles
# ============================================================

rms_original <- adjacent_rms$
  RMS_distance[
    adjacent_rms$Version ==
      "Original"
  ]

rms_conditional <- adjacent_rms$
  RMS_distance[
    adjacent_rms$Version ==
      "Conditional"
  ]

rms_aggressive <- adjacent_rms$
  RMS_distance[
    adjacent_rms$Version ==
      "Aggressive"
  ]


adjacent_rms_concordance <- data.frame(
  Comparison = c(
    "Original_vs_Conditional",
    "Original_vs_Aggressive"
  ),

  Pearson_r = c(
    cor(
      rms_original,
      rms_conditional,
      method = "pearson"
    ),
    cor(
      rms_original,
      rms_aggressive,
      method = "pearson"
    )
  ),

  Spearman_rho = c(
    cor(
      rms_original,
      rms_conditional,
      method = "spearman"
    ),
    cor(
      rms_original,
      rms_aggressive,
      method = "spearman"
    )
  ),

  stringsAsFactors = FALSE
)


# ============================================================
# 19. Coverage associations with first five PCs
# ============================================================

coverage_pc_associations <- list()

idx <- 1

for (version_name in c(
  "Original",
  "Conditional",
  "Aggressive"
)) {

  coord <- switch(
    version_name,
    Original =
      pca_original$coordinates,
    Conditional =
      pca_conditional$coordinates,
    Aggressive =
      pca_aggressive$coordinates
  )

  coord <- coord[
    match(
      final_samples,
      coord$sampleID
    ),
    ,
    drop = FALSE
  ]

  for (pc_name in paste0(
    "PC",
    1:5
  )) {

    for (metric_name in c(
      "z_ratio",
      "z_profile"
    )) {

      rho <- suppressWarnings(
        cor(
          model_data[[metric_name]],
          coord[[pc_name]],
          method = "spearman"
        )
      )

      coverage_pc_associations[[idx]] <- data.frame(
        Version = version_name,
        PCA_axis = pc_name,
        Coverage_metric =
          metric_name,
        Spearman_rho =
          rho,
        stringsAsFactors = FALSE
      )

      idx <- idx + 1
    }
  }
}

coverage_pc_associations <- do.call(
  rbind,
  coverage_pc_associations
)

rownames(
  coverage_pc_associations
) <- NULL


# ============================================================
# 20. PCA centroid coordinates for plotting
# ============================================================

add_centroids_to_pca <- function(
    coords) {

  cent <- do.call(
    rbind,
    lapply(
      groups,
      function(g) {

        x <- coords[
          coords$group == g,
          ,
          drop = FALSE
        ]

        data.frame(
          group = g,
          PC1 = mean(
            x$PC1
          ),
          PC2 = mean(
            x$PC2
          ),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  cent$group <- factor(
    cent$group,
    levels = groups
  )

  cent
}


pca_cent_original <-
  add_centroids_to_pca(
    pca_original$coordinates
  )

pca_cent_conditional <-
  add_centroids_to_pca(
    pca_conditional$coordinates
  )

pca_cent_aggressive <-
  add_centroids_to_pca(
    pca_aggressive$coordinates
  )


# ============================================================
# 21. Plot theme
# ============================================================

theme_pub <- theme_classic(
  base_size = 10
) +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 16
    ),
    plot.title = element_text(
      face = "bold",
      size = 10.8
    ),
    axis.title = element_text(
      size = 9.5
    ),
    axis.text = element_text(
      size = 8.3
    ),
    legend.title = element_text(
      size = 8.5
    ),
    legend.text = element_text(
      size = 7.5
    ),
    legend.key.height = unit(
      0.33,
      "cm"
    ),
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )


# ============================================================
# 22. PCA plot helper
# ============================================================

make_pca_plot <- function(
    pca_result,
    centroid_df,
    title,
    tag) {

  coords <-
    pca_result$coordinates

  ggplot(
    coords,
    aes(
      x = PC1,
      y = PC2,
      colour = group
    )
  ) +

    geom_path(
      data = centroid_df,
      aes(
        x = PC1,
        y = PC2,
        group = 1
      ),
      inherit.aes = FALSE,
      colour = "grey35",
      linewidth = 0.75,
      arrow = arrow(
        length = unit(
          0.12,
          "cm"
        ),
        type = "closed"
      )
    ) +

    geom_point(
      size = 2.2,
      alpha = 0.72
    ) +

    geom_point(
      data = centroid_df,
      aes(
        x = PC1,
        y = PC2,
        fill = group
      ),
      inherit.aes = FALSE,
      shape = 21,
      size = 4.0,
      stroke = 0.7,
      colour = "black"
    ) +

    scale_colour_manual(
      values = stage_colors,
      breaks = groups,
      labels = unname(
        stage_labels[
          groups
        ]
      ),
      name = "Stage"
    ) +

    scale_fill_manual(
      values = stage_colors,
      guide = "none"
    ) +

    labs(
      tag = tag,
      title = title,
      x = paste0(
        "PC1 (",
        sprintf(
          "%.2f",
          pca_result$variance[1]
        ),
        "%)"
      ),
      y = paste0(
        "PC2 (",
        sprintf(
          "%.2f",
          pca_result$variance[2]
        ),
        "%)"
      )
    ) +

    theme_pub +

    theme(
      legend.position = "right"
    )
}


# ============================================================
# 23. Main diagnostic panels
# ============================================================

pA <- make_pca_plot(
  pca_original,
  pca_cent_original,
  "Original VST expression space",
  "A"
)

pB <- make_pca_plot(
  pca_conditional,
  pca_cent_conditional,
  "Stage-preserving coverage adjustment",
  "B"
)

pC <- make_pca_plot(
  pca_aggressive,
  pca_cent_aggressive,
  "Aggressive coverage-only adjustment",
  "C"
)


adjacent_rms$Version <- factor(
  adjacent_rms$Version,
  levels = c(
    "Original",
    "Conditional",
    "Aggressive"
  )
)

adjacent_rms$Transition <- factor(
  adjacent_rms$Transition,
  levels = unique(
    adjacent_rms$Transition[
      adjacent_rms$Version ==
        "Original"
    ]
  )
)


pD <- ggplot(
  adjacent_rms,
  aes(
    x = Transition,
    y = RMS_distance,
    group = Version,
    linetype = Version,
    shape = Version
  )
) +

  geom_line(
    linewidth = 0.75
  ) +

  geom_point(
    size = 2.3
  ) +

  scale_linetype_manual(
    values = c(
      Original = "solid",
      Conditional = "dashed",
      Aggressive = "dotdash"
    ),
    labels = c(
      Original = "Original",
      Conditional =
        "Stage-preserving",
      Aggressive =
        "Coverage-only"
    ),
    name = "Expression space"
  ) +

  scale_shape_manual(
    values = c(
      Original = 16,
      Conditional = 17,
      Aggressive = 15
    ),
    labels = c(
      Original = "Original",
      Conditional =
        "Stage-preserving",
      Aggressive =
        "Coverage-only"
    ),
    name = "Expression space"
  ) +

  labs(
    tag = "D",
    title =
      "Adjacent-stage centroid distances",
    x = NULL,
    y =
      "Full-VST centroid RMS distance"
  ) +

  theme_pub +

  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7
    ),
    legend.position = "top"
  )


# ============================================================
# 24. Distance-concordance secondary plots
# ============================================================

distance_df_conditional <- data.frame(
  Original_distance =
    upper_original,
  Adjusted_distance =
    upper_conditional
)

distance_df_aggressive <- data.frame(
  Original_distance =
    upper_original,
  Adjusted_distance =
    upper_aggressive
)


pE <- ggplot(
  distance_df_conditional,
  aes(
    x = Original_distance,
    y = Adjusted_distance
  )
) +

  geom_point(
    size = 0.9,
    alpha = 0.22
  ) +

  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    colour = "black",
    linewidth = 0.7
  ) +

  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = paste0(
      "Spearman ρ = ",
      sprintf(
        "%.3f",
        distance_concordance$
          Spearman_rho[
            distance_concordance$
              Comparison ==
              "Original_vs_Conditional"
          ]
      )
    ),
    hjust = 1.05,
    vjust = -0.7,
    size = 3
  ) +

  labs(
    title =
      "Original vs stage-preserving sample distances",
    x = "Original Euclidean distance",
    y =
      "Stage-preserving adjusted distance"
  ) +

  theme_pub


pF <- ggplot(
  distance_df_aggressive,
  aes(
    x = Original_distance,
    y = Adjusted_distance
  )
) +

  geom_point(
    size = 0.9,
    alpha = 0.22
  ) +

  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    colour = "black",
    linewidth = 0.7
  ) +

  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = paste0(
      "Spearman ρ = ",
      sprintf(
        "%.3f",
        distance_concordance$
          Spearman_rho[
            distance_concordance$
              Comparison ==
              "Original_vs_Aggressive"
          ]
      )
    ),
    hjust = 1.05,
    vjust = -0.7,
    size = 3
  ) +

  labs(
    title =
      "Original vs coverage-only sample distances",
    x = "Original Euclidean distance",
    y =
      "Coverage-only adjusted distance"
  ) +

  theme_pub


# ============================================================
# 25. Saving helpers
# ============================================================

open_tiff <- function(
    file,
    width,
    height) {

  if (requireNamespace(
    "ragg",
    quietly = TRUE
  )) {

    ragg::agg_tiff(
      filename = file,
      width = width,
      height = height,
      units = "in",
      res = 600,
      compression = "lzw"
    )

  } else {

    grDevices::tiff(
      filename = file,
      width = width,
      height = height,
      units = "in",
      res = 600,
      type = "cairo",
      compression = "lzw"
    )
  }
}


save_main_figure <- function(
    pA,
    pB,
    pC,
    pD,
    width = 13,
    height = 10) {

  draw <- function() {

    grid.newpage()

    pushViewport(
      viewport(
        layout = grid.layout(
          nrow = 2,
          ncol = 2,
          widths = unit(
            c(
              1,
              1
            ),
            "null"
          ),
          heights = unit(
            c(
              1,
              1
            ),
            "null"
          )
        )
      )
    )

    print(
      pA,
      vp = viewport(
        layout.pos.row = 1,
        layout.pos.col = 1
      )
    )

    print(
      pB,
      vp = viewport(
        layout.pos.row = 1,
        layout.pos.col = 2
      )
    )

    print(
      pC,
      vp = viewport(
        layout.pos.row = 2,
        layout.pos.col = 1
      )
    )

    print(
      pD,
      vp = viewport(
        layout.pos.row = 2,
        layout.pos.col = 2
      )
    )
  }


  pdf_file <- file.path(
    outdir,
    "12C_coverage_adjusted_PCA_sensitivity.pdf"
  )

  grDevices::cairo_pdf(
    filename = pdf_file,
    width = width,
    height = height
  )

  draw()
  dev.off()


  tif_file <- file.path(
    outdir,
    paste0(
      "12C_coverage_adjusted_PCA_",
      "sensitivity_600dpi.tiff"
    )
  )

  open_tiff(
    tif_file,
    width,
    height
  )

  draw()
  dev.off()
}


save_distance_figure <- function(
    pE,
    pF,
    width = 11,
    height = 5) {

  draw <- function() {

    grid.newpage()

    pushViewport(
      viewport(
        layout = grid.layout(
          1,
          2
        )
      )
    )

    print(
      pE,
      vp = viewport(
        layout.pos.row = 1,
        layout.pos.col = 1
      )
    )

    print(
      pF,
      vp = viewport(
        layout.pos.row = 1,
        layout.pos.col = 2
      )
    )
  }

  grDevices::cairo_pdf(
    filename = file.path(
      outdir,
      "12C_sample_distance_concordance.pdf"
    ),
    width = width,
    height = height
  )

  draw()
  dev.off()
}


# ============================================================
# 26. Save figures
# ============================================================

save_main_figure(
  pA,
  pB,
  pC,
  pD
)

save_distance_figure(
  pE,
  pF
)


# ============================================================
# 27. Combine PCA coordinates
# ============================================================

pca_coordinates_all <- rbind(
  pca_original$coordinates,
  pca_conditional$coordinates,
  pca_aggressive$coordinates
)


# ============================================================
# 28. Save analysis tables
# ============================================================

write.table(
  permanova_summary,
  file = file.path(
    outdir,
    "12C_PERMANOVA_PERMDISP_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  distance_concordance,
  file = file.path(
    outdir,
    "12C_sample_distance_concordance.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  centroid_distance_concordance,
  file = file.path(
    outdir,
    "12C_stage_centroid_distance_concordance.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  adjacent_rms,
  file = file.path(
    outdir,
    "12C_adjacent_stage_centroid_RMS.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  adjacent_rms_concordance,
  file = file.path(
    outdir,
    "12C_adjacent_stage_RMS_concordance.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  coverage_pc_associations,
  file = file.path(
    outdir,
    "12C_coverage_vs_PC1_to_PC5_associations.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  pca_coordinates_all,
  file = file.path(
    outdir,
    "12C_PCA_coordinates_original_and_adjusted.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 29. Save adjusted expression matrices
# gzipped because matrices are large
# ============================================================

write_matrix_gz <- function(
    mat,
    file) {

  con <- gzfile(
    file,
    open = "wt"
  )

  out <- data.frame(
    GeneID = rownames(mat),
    mat,
    check.names = FALSE
  )

  write.table(
    out,
    file = con,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  close(con)
}


write_matrix_gz(
  vst_conditional,
  file.path(
    outdir,
    "12C_VST_stage_preserving_coverage_adjusted.tsv.gz"
  )
)


write_matrix_gz(
  vst_aggressive,
  file.path(
    outdir,
    "12C_VST_coverage_only_adjusted.tsv.gz"
  )
)


# ============================================================
# 30. PCA variance summary
# ============================================================

pca_variance_summary <- data.frame(
  Version = c(
    "Original",
    "Conditional",
    "Aggressive"
  ),

  PC1_percent = c(
    pca_original$variance[1],
    pca_conditional$variance[1],
    pca_aggressive$variance[1]
  ),

  PC2_percent = c(
    pca_original$variance[2],
    pca_conditional$variance[2],
    pca_aggressive$variance[2]
  ),

  PC3_percent = c(
    pca_original$variance[3],
    pca_conditional$variance[3],
    pca_aggressive$variance[3]
  ),

  stringsAsFactors = FALSE
)


write.table(
  pca_variance_summary,
  file = file.path(
    outdir,
    "12C_PCA_variance_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 31. Save plot objects / session
# ============================================================

saveRDS(
  list(
    A_original_PCA = pA,
    B_conditional_PCA = pB,
    C_aggressive_PCA = pC,
    D_adjacent_RMS = pD,
    E_distance_conditional = pE,
    F_distance_aggressive = pF
  ),
  file = file.path(
    outdir,
    "12C_plot_objects.rds"
  )
)


sink(
  file.path(
    outdir,
    "12C_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 32. Console output
# ============================================================

cat("\n")
cat("============================================================\n")
cat("PCA VARIANCE SUMMARY\n")
cat("============================================================\n\n")

print(
  pca_variance_summary,
  row.names = FALSE,
  digits = 4
)


cat("\n")
cat("============================================================\n")
cat("PERMANOVA / PERMDISP\n")
cat("============================================================\n\n")

print(
  permanova_summary,
  row.names = FALSE,
  digits = 5
)


cat("\n")
cat("============================================================\n")
cat("SAMPLE-DISTANCE CONCORDANCE\n")
cat("============================================================\n\n")

print(
  distance_concordance,
  row.names = FALSE,
  digits = 5
)


cat("\n")
cat("============================================================\n")
cat("STAGE-CENTROID DISTANCE CONCORDANCE\n")
cat("============================================================\n\n")

print(
  centroid_distance_concordance,
  row.names = FALSE,
  digits = 5
)


cat("\n")
cat("============================================================\n")
cat("ADJACENT-STAGE RMS CONCORDANCE\n")
cat("============================================================\n\n")

print(
  adjacent_rms_concordance,
  row.names = FALSE,
  digits = 5
)


cat("\n")
cat("============================================================\n")
cat("COVERAGE ASSOCIATIONS WITH PC1-PC5\n")
cat("============================================================\n\n")

print(
  coverage_pc_associations,
  row.names = FALSE,
  digits = 4
)


cat("\n")
cat("============================================================\n")
cat("INTERPRETATION FRAMEWORK\n")
cat("============================================================\n")

cat(
  "Original      : frozen VST expression space.\n"
)

cat(
  "Conditional   : coverage effects removed conditional on stage; ",
  "stage effects are deliberately preserved.\n"
)

cat(
  "Aggressive    : all linear signal associated with the two ",
  "coverage metrics is removed without protecting stage.\n"
)

cat(
  "\nThe aggressive adjustment is intentionally conservative ",
  "and may remove genuine biology because coverage and stage ",
  "are correlated.\n"
)

cat(
  "Therefore, it should be interpreted as a sensitivity bound, ",
  "not as the preferred expression matrix.\n"
)


cat("\n")
cat("============================================================\n")
cat("12C COMPLETED\n")
cat("============================================================\n")

cat(
  "Samples analysed : ",
  ncol(vst_mat),
  "\n",
  sep = ""
)

cat(
  "Genes analysed   : ",
  nrow(vst_mat),
  "\n",
  sep = ""
)

cat(
  "R2h_1 excluded   : ",
  !("R2h_1" %in%
      colnames(vst_mat)),
  "\n",
  sep = ""
)

cat(
  "R2h_8 retained   : ",
  "R2h_8" %in%
    colnames(vst_mat),
  "\n",
  sep = ""
)

cat("\nOutput directory:\n")
cat(outdir, "\n\n")

cat("PASS\n")
