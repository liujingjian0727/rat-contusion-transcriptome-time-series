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
cat("13D FASTP-QC-ADJUSTED EXPRESSION-SPACE SENSITIVITY\n")
cat("Final 77-sample 2026 primary cohort\n")
cat("============================================================\n\n")


# ============================================================
# 0. Fixed paths
# ============================================================

metadata_file <-
  "rattus_meta_2026_77samples.tsv"

vst_file <- paste0(
  "02_count_filtering_and_VST_result/",
  "2026_primary/",
  "02_2026_VST_matrix.tsv"
)

qc_file <- paste0(
  "supplementary_figures/",
  "Supplementary_Figure_S3_sequencing_mapping_assignment_QC/",
  "13B_fastp_QC/",
  "Supplementary_Table_S2_2026_primary_77sample_",
  "full_sequencing_QC.tsv"
)

outdir <- paste0(
  "supplementary_figures/",
  "Supplementary_Figure_S3_sequencing_mapping_assignment_QC/",
  "13D_fastp_QC_adjusted_expression_space_sensitivity"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

RANDOM_SEED <- 20260910


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
  qc_file
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
# 3. Read final metadata
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
  ) %in% colnames(meta)
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

if (anyDuplicated(
  meta$sampleID
)) {
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

if (anyDuplicated(
  vst_df$GeneID
)) {
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

storage.mode(
  vst_mat
) <- "numeric"

rownames(vst_mat) <-
  vst_df$GeneID

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
    "Invalid values detected in VST matrix."
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
# 5. Read integrated sequencing QC
# ============================================================

qc <- read.delim(
  qc_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

required_qc_cols <- c(
  "sampleID",
  "Q30_clean_percent",
  "GC_clean_percent",
  "Duplication_rate_percent"
)

if (!all(
  required_qc_cols %in%
    colnames(qc)
)) {
  stop(
    "Integrated QC table lacks required columns."
  )
}

if (anyDuplicated(
  qc$sampleID
)) {
  stop(
    "Duplicated sampleID in QC table."
  )
}

if (!setequal(
  qc$sampleID,
  final_samples
)) {
  stop(
    "QC sample set differs from final metadata."
  )
}

qc <- qc[
  match(
    final_samples,
    qc$sampleID
  ),
  ,
  drop = FALSE
]


# ============================================================
# 6. Build model-data table
# ============================================================

model_data <- data.frame(
  sampleID =
    final_samples,

  group =
    meta$group,

  Q30 =
    qc$Q30_clean_percent,

  GC =
    qc$GC_clean_percent,

  Duplication =
    qc$Duplication_rate_percent,

  stringsAsFactors = FALSE
)

model_data$group <- factor(
  model_data$group,
  levels = groups
)

model_data$z_Q30 <- as.numeric(
  scale(
    model_data$Q30
  )
)

model_data$z_GC <- as.numeric(
  scale(
    model_data$GC
  )
)

model_data$z_Duplication <- as.numeric(
  scale(
    model_data$Duplication
  )
)

if (anyNA(
  model_data[
    ,
    c(
      "Q30",
      "GC",
      "Duplication",
      "z_Q30",
      "z_GC",
      "z_Duplication"
    )
  ]
)) {
  stop(
    "NA detected in QC covariates."
  )
}

rownames(model_data) <-
  model_data$sampleID


# ============================================================
# 7. Cohort sanity checks
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
# 8. QC collinearity diagnostics
# ============================================================

qc_cor_spearman <- cor(
  model_data[
    ,
    c(
      "Q30",
      "GC",
      "Duplication"
    )
  ],
  method = "spearman"
)

qc_cor_pearson <- cor(
  model_data[
    ,
    c(
      "Q30",
      "GC",
      "Duplication"
    )
  ],
  method = "pearson"
)


# ============================================================
# 9. Matrix regression helper
#
# Expression:
# genes x samples
#
# We subtract only selected design-column effects.
# Intercept and all nonselected effects are retained.
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

  effect_matrix <-
    design_matrix[
      ,
      effect_columns,
      drop = FALSE
    ] %*%
    beta[
      effect_columns,
      ,
      drop = FALSE
    ]

  adjusted_Y <- Y -
    effect_matrix

  adjusted <- t(
    adjusted_Y
  )

  rownames(adjusted) <-
    rownames(expression_matrix)

  colnames(adjusted) <-
    colnames(expression_matrix)

  list(
    adjusted = adjusted,
    beta = beta,
    design_rank = design_rank,
    design_columns =
      colnames(design_matrix)
  )
}


# ============================================================
# 10. CONDITIONAL / STAGE-PRESERVING adjustment
#
# expression ~ stage + Q30 + GC + duplication
#
# Only QC-associated effects conditional on stage are removed.
# ============================================================

design_conditional <- model.matrix(
  ~ group +
    z_Q30 +
    z_GC +
    z_Duplication,
  data = model_data
)

conditional_effect_columns <- which(
  colnames(
    design_conditional
  ) %in%
    c(
      "z_Q30",
      "z_GC",
      "z_Duplication"
    )
)

if (length(
  conditional_effect_columns
) != 3) {
  stop(
    "Could not identify three QC covariates ",
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
      conditional_effect_columns
  )

vst_conditional <-
  conditional_result$adjusted


# ============================================================
# 11. AGGRESSIVE / QC-ONLY adjustment
#
# expression ~ Q30 + GC + duplication
#
# All linear signal associated with the three QC variables
# is removed without explicitly protecting sampling stage.
#
# Because these QC metrics are stage-associated, this
# adjustment can also remove genuine biological variation.
# ============================================================

design_aggressive <- model.matrix(
  ~ z_Q30 +
    z_GC +
    z_Duplication,
  data = model_data
)

aggressive_effect_columns <- which(
  colnames(
    design_aggressive
  ) %in%
    c(
      "z_Q30",
      "z_GC",
      "z_Duplication"
    )
)

if (length(
  aggressive_effect_columns
) != 3) {
  stop(
    "Could not identify three QC covariates ",
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
      aggressive_effect_columns
  )

vst_aggressive <-
  aggressive_result$adjusted


# ============================================================
# 12. Dimension checks
# ============================================================

if (!identical(
  dim(vst_mat),
  dim(vst_conditional)
)) {
  stop(
    "Conditional matrix dimension mismatch."
  )
}

if (!identical(
  dim(vst_mat),
  dim(vst_aggressive)
)) {
  stop(
    "Aggressive matrix dimension mismatch."
  )
}

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
# 13. PCA helper
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
    sampleID =
      rownames(
        p$x
      ),

    PC1 =
      p$x[
        ,
        1
      ],

    PC2 =
      p$x[
        ,
        2
      ],

    PC3 =
      p$x[
        ,
        3
      ],

    PC4 =
      p$x[
        ,
        4
      ],

    PC5 =
      p$x[
        ,
        5
      ],

    Version =
      version_name,

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
  "Stage_preserving"
)

pca_aggressive <- run_pca(
  vst_aggressive,
  "QC_only"
)


cat(
  "Original PCA         : PC1 = ",
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
  "Stage-preserving PCA : PC1 = ",
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
  "QC-only PCA          : PC1 = ",
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
# 14. Euclidean sample distances
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

vec_original <- as.vector(
  dist_original
)

vec_conditional <- as.vector(
  dist_conditional
)

vec_aggressive <- as.vector(
  dist_aggressive
)


distance_concordance <- data.frame(
  Comparison = c(
    "Original_vs_Stage_preserving",
    "Original_vs_QC_only"
  ),

  Pearson_r = c(
    cor(
      vec_original,
      vec_conditional,
      method = "pearson"
    ),

    cor(
      vec_original,
      vec_aggressive,
      method = "pearson"
    )
  ),

  Spearman_rho = c(
    cor(
      vec_original,
      vec_conditional,
      method = "spearman"
    ),

    cor(
      vec_original,
      vec_aggressive,
      method = "spearman"
    )
  ),

  Mean_absolute_distance_difference = c(
    mean(
      abs(
        vec_original -
          vec_conditional
      )
    ),

    mean(
      abs(
        vec_original -
          vec_aggressive
      )
    )
  ),

  stringsAsFactors = FALSE
)


# ============================================================
# 15. PERMANOVA / PERMDISP helper
#
# Exactly matches frozen Step 04:
# - Euclidean VST distance
# - 9999 permutations
# - seed 20260910
# - betadisper median
# - bias.adjust TRUE
# ============================================================

run_permanova <- function(
    distance_object,
    version_name) {

  set.seed(
    RANDOM_SEED
  )

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

  set.seed(
    RANDOM_SEED
  )

  dispersion_perm <- vegan::permutest(
    dispersion,
    permutations = 9999
  )

  data.frame(
    Version =
      version_name,

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
    "Stage_preserving"
  ),

  run_permanova(
    dist_aggressive,
    "QC_only"
  )
)


# ============================================================
# 16. Stage centroids in full expression space
# ============================================================

calculate_centroids <- function(
    expression_matrix) {

  result <- matrix(
    NA_real_,
    nrow = length(groups),
    ncol = nrow(
      expression_matrix
    ),
    dimnames = list(
      groups,
      rownames(
        expression_matrix
      )
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

  result
}


centroids_original <-
  calculate_centroids(
    vst_mat
  )

centroids_conditional <-
  calculate_centroids(
    vst_conditional
  )

centroids_aggressive <-
  calculate_centroids(
    vst_aggressive
  )


# ============================================================
# 17. Stage-centroid distance concordance
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
    "Original_vs_Stage_preserving",
    "Original_vs_QC_only"
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
# 18. Adjacent-stage centroid RMS
# ============================================================

calculate_adjacent_rms <- function(
    centroid_matrix,
    version_name) {

  result_list <- vector(
    "list",
    length(groups) - 1
  )

  for (i in seq_len(
    length(groups) - 1
  )) {

    g1 <- groups[i]
    g2 <- groups[i + 1]

    delta <-
      centroid_matrix[
        g2,
      ] -
      centroid_matrix[
        g1,
      ]

    rms <- sqrt(
      mean(
        delta^2
      )
    )

    result_list[i] <- list(
      data.frame(
        Version =
          version_name,

        From =
          g1,

        To =
          g2,

        Transition =
          paste0(
            unname(
              stage_labels[g1]
            ),
            " -> ",
            unname(
              stage_labels[g2]
            )
          ),

        Order =
          i,

        RMS_distance =
          rms,

        stringsAsFactors = FALSE
      )
    )
  }

  do.call(
    rbind,
    result_list
  )
}


adjacent_rms <- rbind(
  calculate_adjacent_rms(
    centroids_original,
    "Original"
  ),

  calculate_adjacent_rms(
    centroids_conditional,
    "Stage_preserving"
  ),

  calculate_adjacent_rms(
    centroids_aggressive,
    "QC_only"
  )
)


rms_original <- adjacent_rms$RMS_distance[
  adjacent_rms$Version ==
    "Original"
]

rms_conditional <- adjacent_rms$RMS_distance[
  adjacent_rms$Version ==
    "Stage_preserving"
]

rms_aggressive <- adjacent_rms$RMS_distance[
  adjacent_rms$Version ==
    "QC_only"
]


adjacent_rms_concordance <- data.frame(
  Comparison = c(
    "Original_vs_Stage_preserving",
    "Original_vs_QC_only"
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
# 19. QC association with PC1-PC5 after adjustment
# ============================================================

qc_pc_results <- list()
counter <- 1

versions <- c(
  "Original",
  "Stage_preserving",
  "QC_only"
)

for (version_name in versions) {

  if (version_name == "Original") {

    coords <-
      pca_original$coordinates

  } else if (
    version_name ==
      "Stage_preserving"
  ) {

    coords <-
      pca_conditional$coordinates

  } else {

    coords <-
      pca_aggressive$coordinates
  }

  coords <- coords[
    match(
      final_samples,
      coords$sampleID
    ),
    ,
    drop = FALSE
  ]


  for (pc in paste0(
    "PC",
    1:5
  )) {

    for (metric in c(
      "z_Q30",
      "z_GC",
      "z_Duplication"
    )) {

      x <- model_data[
        ,
        metric
      ]

      y <- coords[
        ,
        pc
      ]

      sp <- suppressWarnings(
        cor.test(
          x,
          y,
          method = "spearman",
          exact = FALSE
        )
      )

      qc_pc_results[counter] <- list(
        data.frame(
          Version =
            version_name,

          PCA_axis =
            pc,

          QC_metric =
            metric,

          Spearman_rho =
            unname(
              sp$estimate
            ),

          P_value =
            sp$p.value,

          stringsAsFactors = FALSE
        )
      )

      counter <- counter + 1
    }
  }
}

qc_pc_associations <- do.call(
  rbind,
  qc_pc_results
)

qc_pc_associations$BH_FDR <- ave(
  qc_pc_associations$P_value,
  qc_pc_associations$Version,
  FUN = function(x) {
    p.adjust(
      x,
      method = "BH"
    )
  }
)


# ============================================================
# 20. PCA variance summary
# ============================================================

pca_variance_summary <- data.frame(
  Version = c(
    "Original",
    "Stage_preserving",
    "QC_only"
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


# ============================================================
# 21. PCA centroids for plotting
# ============================================================

make_pca_centroids <- function(
    coords) {

  result_list <- vector(
    "list",
    length(groups)
  )

  for (i in seq_along(
    groups
  )) {

    g <- groups[i]

    tmp <- coords[
      coords$group == g,
      ,
      drop = FALSE
    ]

    result_list[i] <- list(
      data.frame(
        group =
          g,

        PC1 =
          mean(
            tmp$PC1
          ),

        PC2 =
          mean(
            tmp$PC2
          ),

        stringsAsFactors = FALSE
      )
    )
  }

  result <- do.call(
    rbind,
    result_list
  )

  result$group <- factor(
    result$group,
    levels = groups
  )

  result
}


pca_centroids_original <-
  make_pca_centroids(
    pca_original$coordinates
  )

pca_centroids_conditional <-
  make_pca_centroids(
    pca_conditional$coordinates
  )

pca_centroids_aggressive <-
  make_pca_centroids(
    pca_aggressive$coordinates
  )


# ============================================================
# 22. Plot theme
# ============================================================

theme_pub <- theme_classic(
  base_size = 10
) +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 17
    ),

    plot.title = element_text(
      face = "bold",
      size = 10.8
    ),

    plot.subtitle = element_text(
      size = 8.2,
      colour = "grey35"
    ),

    axis.title = element_text(
      size = 9.5
    ),

    axis.text = element_text(
      size = 8
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
# 23. PCA plotting helper
# ============================================================

make_pca_plot <- function(
    pca_result,
    centroid_df,
    tag,
    title,
    subtitle) {

  coords <- pca_result$coordinates

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
      linewidth = 0.70,
      arrow = arrow(
        length = unit(
          0.11,
          "cm"
        ),
        type = "closed"
      )
    ) +

    geom_point(
      size = 2.15,
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
      size = 3.7,
      stroke = 0.65,
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
      name = "Sampling stage"
    ) +

    scale_fill_manual(
      values = stage_colors,
      guide = "none"
    ) +

    labs(
      tag = tag,
      title = title,
      subtitle = subtitle,

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
# 24. Panels A-C
# ============================================================

pA <- make_pca_plot(
  pca_original,
  pca_centroids_original,
  "A",
  "Original VST expression space",
  "No fastp QC adjustment"
)


pB <- make_pca_plot(
  pca_conditional,
  pca_centroids_conditional,
  "B",
  "Stage-preserving fastp-QC adjustment",
  "QC-associated effects removed conditional on sampling stage"
)


pC <- make_pca_plot(
  pca_aggressive,
  pca_centroids_aggressive,
  "C",
  "QC-only adjustment",
  "All linear signal associated with Q30, GC and duplication removed"
)


# ============================================================
# 25. Panel D: stage PERMANOVA sensitivity
# ============================================================

perm_plot <- permanova_summary

perm_plot$Display_version <- factor(
  perm_plot$Version,
  levels = c(
    "Original",
    "Stage_preserving",
    "QC_only"
  ),
  labels = c(
    "Original",
    "Stage-preserving",
    "QC-only"
  )
)

perm_plot$Top_label <- paste0(
  "R² = ",
  sprintf(
    "%.3f",
    perm_plot$PERMANOVA_R2
  ),
  "\nPERMANOVA P ",
  ifelse(
    perm_plot$PERMANOVA_P < 0.001,
    "<0.001",
    paste0(
      "= ",
      sprintf(
        "%.3f",
        perm_plot$PERMANOVA_P
      )
    )
  )
)

x_labels <- setNames(
  paste0(
    as.character(
      perm_plot$Display_version
    ),
    "\nPERMDISP P = ",
    sprintf(
      "%.3f",
      perm_plot$PERMDISP_P
    )
  ),
  as.character(
    perm_plot$Display_version
  )
)


pD <- ggplot(
  perm_plot,
  aes(
    x = Display_version,
    y = PERMANOVA_R2
  )
) +

  geom_col(
    width = 0.62,
    fill = "grey67",
    colour = "grey20",
    linewidth = 0.45
  ) +

  geom_text(
    aes(
      label = Top_label
    ),
    vjust = -0.35,
    size = 2.9,
    lineheight = 0.95
  ) +

  scale_x_discrete(
    labels = x_labels
  ) +

  scale_y_continuous(
    limits = c(
      0,
      max(
        perm_plot$PERMANOVA_R2
      ) * 1.22
    ),
    expand = c(
      0,
      0
    )
  ) +

  labs(
    tag = "D",
    title =
      "Sampling-stage structure after fastp-QC adjustment",
    subtitle =
      "PERMANOVA based on full VST Euclidean expression space",
    x = NULL,
    y =
      "PERMANOVA R² for sampling stage"
  ) +

  theme_pub +

  theme(
    axis.text.x = element_text(
      size = 7.8,
      face = "bold",
      lineheight = 1.05
    )
  )


# ============================================================
# 26. Save combined figure
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


draw_combined <- function() {

  grid.newpage()

  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 2,
        ncol = 2
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
  "13D_fastp_QC_adjusted_expression_space_sensitivity.pdf"
)

grDevices::cairo_pdf(
  filename = pdf_file,
  width = 13,
  height = 10
)

draw_combined()
dev.off()


tiff_file <- file.path(
  outdir,
  paste0(
    "13D_fastp_QC_adjusted_expression_space_",
    "sensitivity_600dpi.tiff"
  )
)

open_tiff(
  tiff_file,
  13,
  10
)

draw_combined()
dev.off()


# ============================================================
# 27. Save adjusted matrices
# ============================================================

write_matrix_gz <- function(
    matrix_object,
    output_file) {

  out <- data.frame(
    GeneID =
      rownames(
        matrix_object
      ),
    matrix_object,
    check.names = FALSE
  )

  con <- gzfile(
    output_file,
    open = "wt"
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
    "13D_VST_stage_preserving_fastp_adjusted.tsv.gz"
  )
)

write_matrix_gz(
  vst_aggressive,
  file.path(
    outdir,
    "13D_VST_QC_only_adjusted.tsv.gz"
  )
)


# ============================================================
# 28. Save PCA coordinates
# ============================================================

pca_coordinates_all <- rbind(
  pca_original$coordinates,
  pca_conditional$coordinates,
  pca_aggressive$coordinates
)

write.table(
  pca_coordinates_all,
  file = file.path(
    outdir,
    "13D_PCA_coordinates_original_and_adjusted.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 29. Save summary tables
# ============================================================

write.table(
  pca_variance_summary,
  file = file.path(
    outdir,
    "13D_PCA_variance_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  permanova_summary,
  file = file.path(
    outdir,
    "13D_PERMANOVA_PERMDISP_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  distance_concordance,
  file = file.path(
    outdir,
    "13D_sample_distance_concordance.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  centroid_distance_concordance,
  file = file.path(
    outdir,
    "13D_stage_centroid_distance_concordance.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  adjacent_rms,
  file = file.path(
    outdir,
    "13D_adjacent_stage_centroid_RMS.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  adjacent_rms_concordance,
  file = file.path(
    outdir,
    "13D_adjacent_stage_RMS_concordance.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  qc_pc_associations,
  file = file.path(
    outdir,
    "13D_QC_vs_PC1_to_PC5_after_adjustment.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  qc_cor_spearman,
  file = file.path(
    outdir,
    "13D_QC_metric_Spearman_matrix.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


write.table(
  qc_cor_pearson,
  file = file.path(
    outdir,
    "13D_QC_metric_Pearson_matrix.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 30. Save session info and plot objects
# ============================================================

saveRDS(
  list(
    Panel_A = pA,
    Panel_B = pB,
    Panel_C = pC,
    Panel_D = pD
  ),
  file = file.path(
    outdir,
    "13D_plot_objects.rds"
  )
)

sink(
  file.path(
    outdir,
    "13D_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 31. Console report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("QC METRIC SPEARMAN CORRELATION\n")
cat("============================================================\n\n")

print(
  qc_cor_spearman,
  digits = 4
)


cat("\n")
cat("============================================================\n")
cat("PCA VARIANCE SUMMARY\n")
cat("============================================================\n\n")

print(
  pca_variance_summary,
  row.names = FALSE,
  digits = 5
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
cat("QC ASSOCIATIONS WITH PC1-PC5 AFTER ADJUSTMENT\n")
cat("============================================================\n\n")

print(
  qc_pc_associations,
  row.names = FALSE,
  digits = 4
)


cat("\n")
cat("============================================================\n")
cat("INTERPRETATION FRAMEWORK\n")
cat("============================================================\n")

cat(
  "Original:\n",
  "  Frozen VST expression space.\n\n",
  sep = ""
)

cat(
  "Stage_preserving:\n",
  "  Q30-, GC- and duplication-associated components are removed\n",
  "  conditional on sampling stage. Stage effects are deliberately\n",
  "  retained.\n\n",
  sep = ""
)

cat(
  "QC_only:\n",
  "  All linear expression components associated with Q30, GC and\n",
  "  duplication are removed without explicitly preserving stage.\n\n",
  sep = ""
)

cat(
  "Because the fastp QC variables are themselves associated with\n",
  "sampling stage, QC-only adjustment may also remove genuine\n",
  "stage-related biology. It is therefore a conservative sensitivity\n",
  "bound rather than a preferred expression matrix.\n"
)


cat("\n")
cat("============================================================\n")
cat("13D COMPLETED\n")
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

cat("Primary figure:\n")
cat(pdf_file, "\n")
cat(tiff_file, "\n\n")

cat("PASS\n")
