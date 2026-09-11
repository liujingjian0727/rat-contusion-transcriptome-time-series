#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  bitmapType = "cairo"
)

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

cat("\n")
cat("============================================================\n")
cat("12B GENE-BODY COVERAGE VS PCA DIAGNOSTIC\n")
cat("Final 77-sample 2026 primary cohort\n")
cat("============================================================\n\n")


# ============================================================
# 0. Input / output
# ============================================================

metadata_file <-
  "rattus_meta_2026_77samples.tsv"

coverage_metric_file <- paste0(
  "supplementary_figures/",
  "Supplementary_Figure_S1_gene_body_coverage/",
  "S1_per_sample_gene_body_coverage_metrics.tsv"
)

pca_file <- paste0(
  "final_submission_figures/",
  "panels/",
  "Figure_2B_PCA_coordinates.tsv"
)

vst_file <- paste0(
  "02_count_filtering_and_VST_result/",
  "2026_primary/",
  "02_2026_VST_matrix.tsv"
)

outdir <- paste0(
  "supplementary_figures/",
  "Supplementary_Figure_S1_gene_body_coverage/",
  "12B_coverage_vs_PCA_diagnostic"
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
# 2. Check required files
# ============================================================

required_files <- c(
  metadata_file,
  coverage_metric_file
)

missing_required <- required_files[
  !file.exists(required_files)
]

if (length(missing_required) > 0) {

  stop(
    "Missing required file(s):\n",
    paste(
      missing_required,
      collapse = "\n"
    )
  )
}

cat("Required input files found: PASS\n")


# ============================================================
# 3. Read canonical final metadata
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
# 4. Read gene-body coverage metrics
# ============================================================

coverage <- read.delim(
  coverage_metric_file,
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
    "Coverage metric table is missing required columns."
  )
}

if (anyDuplicated(
  coverage$sampleID
)) {

  stop(
    "Duplicated sampleID in coverage metric table."
  )
}

if (!setequal(
  coverage$sampleID,
  final_samples
)) {

  stop(
    "Coverage metric sample set differs ",
    "from final 77-sample metadata."
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
# 5. Obtain frozen PCA coordinates
#
# Prefer Figure 2 PCA source table.
# If absent, recompute exactly from frozen VST matrix.
# ============================================================

if (file.exists(pca_file)) {

  cat(
    "Using frozen Figure 2 PCA coordinates:\n",
    pca_file,
    "\n\n",
    sep = ""
  )

  pca_df <- read.delim(
    pca_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    stringsAsFactors = FALSE
  )

  if (!all(
    c(
      "sampleID",
      "PC1",
      "PC2"
    ) %in%
      colnames(pca_df)
  )) {

    stop(
      "Frozen PCA file lacks sampleID/PC1/PC2."
    )
  }

} else {

  cat(
    "Frozen PCA coordinate table not found.\n",
    "Recomputing PCA from frozen VST matrix.\n\n",
    sep = ""
  )

  if (!file.exists(vst_file)) {
    stop(
      "Neither frozen PCA coordinates nor VST matrix found."
    )
  }

  vst_df <- read.delim(
    vst_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    stringsAsFactors = FALSE
  )

  colnames(vst_df)[1] <- "GeneID"

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

  pca <- prcomp(
    t(vst_mat),
    center = TRUE,
    scale. = FALSE
  )

  pca_var <- 100 *
    pca$sdev^2 /
    sum(
      pca$sdev^2
    )

  pca_df <- data.frame(
    sampleID = rownames(
      pca$x
    ),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )

  cat(
    "Recomputed PCA variance: PC1 = ",
    sprintf(
      "%.2f%%",
      pca_var[1]
    ),
    "; PC2 = ",
    sprintf(
      "%.2f%%",
      pca_var[2]
    ),
    "\n",
    sep = ""
  )
}


# ============================================================
# 6. Align all 77 samples
# ============================================================

if (!setequal(
  pca_df$sampleID,
  final_samples
)) {

  stop(
    "PCA sample set differs from final metadata."
  )
}

pca_df <- pca_df[
  match(
    final_samples,
    pca_df$sampleID
  ),
  ,
  drop = FALSE
]


diagnostic <- data.frame(
  sampleID = final_samples,

  group = meta$group,

  Stage = unname(
    stage_labels[
      as.character(
        meta$group
      )
    ]
  ),

  Ratio_3prime_to_5prime =
    coverage$
      Ratio_3prime_to_5prime,

  Profile_Spearman_to_median =
    coverage$
      Spearman_to_overall_median,

  PC1 = pca_df$PC1,

  PC2 = pca_df$PC2,

  stringsAsFactors = FALSE
)


if (anyNA(
  diagnostic[
    ,
    c(
      "Ratio_3prime_to_5prime",
      "Profile_Spearman_to_median",
      "PC1",
      "PC2"
    )
  ]
)) {

  stop(
    "NA detected in diagnostic variables."
  )
}


cat(
  "Aligned diagnostic samples : ",
  nrow(diagnostic),
  "\n",
  sep = ""
)

cat(
  "R2h_1 present             : ",
  "R2h_1" %in%
    diagnostic$sampleID,
  "\n",
  sep = ""
)

cat(
  "R2h_8 present             : ",
  "R2h_8" %in%
    diagnostic$sampleID,
  "\n\n",
  sep = ""
)


# ============================================================
# 7. Helper functions
# ============================================================

safe_cor_test <- function(
    x,
    y,
    method) {

  suppressWarnings(
    cor.test(
      x,
      y,
      method = method,
      exact = FALSE
    )
  )
}


format_p <- function(p) {

  if (is.na(p)) {
    return("NA")
  }

  if (p < 0.001) {
    return("<0.001")
  }

  sprintf(
    "%.3f",
    p
  )
}


extract_lm_stats <- function(
    x,
    y) {

  model <- lm(
    y ~ x
  )

  sm <- summary(
    model
  )

  coef_tab <- sm$coefficients

  data.frame(
    Linear_model_R2 =
      unname(
        sm$r.squared
      ),

    Linear_model_adjusted_R2 =
      unname(
        sm$adj.r.squared
      ),

    Linear_model_slope =
      unname(
        coef_tab[
          "x",
          "Estimate"
        ]
      ),

    Linear_model_slope_P =
      unname(
        coef_tab[
          "x",
          "Pr(>|t|)"
        ]
      ),

    stringsAsFactors = FALSE
  )
}


# ============================================================
# 8. Pairwise correlation / regression diagnostics
# ============================================================

pair_specs <- data.frame(
  Coverage_metric = c(
    "Ratio_3prime_to_5prime",
    "Ratio_3prime_to_5prime",
    "Profile_Spearman_to_median",
    "Profile_Spearman_to_median"
  ),

  PCA_axis = c(
    "PC1",
    "PC2",
    "PC1",
    "PC2"
  ),

  stringsAsFactors = FALSE
)


pairwise_results <- list()

for (i in seq_len(
  nrow(pair_specs)
)) {

  x_name <-
    pair_specs$
      Coverage_metric[i]

  y_name <-
    pair_specs$
      PCA_axis[i]

  x <- diagnostic[[x_name]]

  y <- diagnostic[[y_name]]

  pearson <- safe_cor_test(
    x,
    y,
    "pearson"
  )

  spearman <- safe_cor_test(
    x,
    y,
    "spearman"
  )

  lm_stats <- extract_lm_stats(
    x,
    y
  )

  pairwise_results[[i]] <- data.frame(
    Coverage_metric = x_name,
    PCA_axis = y_name,
    N = length(x),

    Pearson_r =
      unname(
        pearson$estimate
      ),

    Pearson_P =
      pearson$p.value,

    Spearman_rho =
      unname(
        spearman$estimate
      ),

    Spearman_P =
      spearman$p.value,

    lm_stats,

    stringsAsFactors = FALSE
  )
}

pairwise_results <- do.call(
  rbind,
  pairwise_results
)

rownames(
  pairwise_results
) <- NULL


# ============================================================
# 9. Relationship between the two coverage metrics
# ============================================================

coverage_metric_correlation <-
  safe_cor_test(
    diagnostic$
      Ratio_3prime_to_5prime,
    diagnostic$
      Profile_Spearman_to_median,
    "spearman"
  )

coverage_metric_relationship <-
  data.frame(
    Metric1 =
      "Ratio_3prime_to_5prime",

    Metric2 =
      "Profile_Spearman_to_median",

    N = nrow(
      diagnostic
    ),

    Spearman_rho =
      unname(
        coverage_metric_correlation$
          estimate
      ),

    Spearman_P =
      coverage_metric_correlation$
        p.value,

    stringsAsFactors = FALSE
  )


# ============================================================
# 10. Standardized metrics for multivariable models
# ============================================================

diagnostic$z_ratio <-
  as.numeric(
    scale(
      diagnostic$
        Ratio_3prime_to_5prime
    )
  )

diagnostic$z_profile_similarity <-
  as.numeric(
    scale(
      diagnostic$
        Profile_Spearman_to_median
    )
  )


# ============================================================
# 11. Coverage-only multivariable models
#
# PC ~ 3'/5' ratio + profile similarity
# ============================================================

coverage_only_models <- list()

for (axis_name in c(
  "PC1",
  "PC2"
)) {

  form <- as.formula(
    paste0(
      axis_name,
      " ~ z_ratio + z_profile_similarity"
    )
  )

  fit <- lm(
    form,
    data = diagnostic
  )

  sm <- summary(
    fit
  )

  fstat <- sm$fstatistic

  overall_p <- pf(
    fstat[1],
    fstat[2],
    fstat[3],
    lower.tail = FALSE
  )

  coverage_only_models[[axis_name]] <- data.frame(
    PCA_axis = axis_name,

    R2 =
      sm$r.squared,

    Adjusted_R2 =
      sm$adj.r.squared,

    Overall_model_P =
      overall_p,

    Ratio_standardized_beta =
      coef(
        fit
      )[
        "z_ratio"
      ],

    Ratio_P =
      coef(
        sm
      )[
        "z_ratio",
        "Pr(>|t|)"
      ],

    Profile_similarity_standardized_beta =
      coef(
        fit
      )[
        "z_profile_similarity"
      ],

    Profile_similarity_P =
      coef(
        sm
      )[
        "z_profile_similarity",
        "Pr(>|t|)"
      ],

    stringsAsFactors = FALSE
  )
}

coverage_only_models <- do.call(
  rbind,
  coverage_only_models
)

rownames(
  coverage_only_models
) <- NULL


# ============================================================
# 12. Stage-adjusted incremental tests
#
# Critical question:
# After accounting for sampling stage,
# do coverage metrics explain additional PC variation?
#
# Reduced:
#   PC ~ group
#
# Full:
#   PC ~ group + z_ratio + z_profile_similarity
#
# No causal interpretation is implied.
# ============================================================

stage_adjusted_results <- list()

for (axis_name in c(
  "PC1",
  "PC2"
)) {

  reduced_formula <-
    as.formula(
      paste0(
        axis_name,
        " ~ group"
      )
    )

  full_formula <-
    as.formula(
      paste0(
        axis_name,
        " ~ group + ",
        "z_ratio + ",
        "z_profile_similarity"
      )
    )

  reduced_fit <- lm(
    reduced_formula,
    data = diagnostic
  )

  full_fit <- lm(
    full_formula,
    data = diagnostic
  )

  comparison <- anova(
    reduced_fit,
    full_fit
  )

  reduced_summary <- summary(
    reduced_fit
  )

  full_summary <- summary(
    full_fit
  )

  stage_adjusted_results[[axis_name]] <- data.frame(
    PCA_axis = axis_name,

    Stage_only_R2 =
      reduced_summary$
        r.squared,

    Stage_only_adjusted_R2 =
      reduced_summary$
        adj.r.squared,

    Stage_plus_coverage_R2 =
      full_summary$
        r.squared,

    Stage_plus_coverage_adjusted_R2 =
      full_summary$
        adj.r.squared,

    Delta_R2 =
      full_summary$
        r.squared -
      reduced_summary$
        r.squared,

    Delta_adjusted_R2 =
      full_summary$
        adj.r.squared -
      reduced_summary$
        adj.r.squared,

    Incremental_F =
      comparison$
        F[2],

    Incremental_P =
      comparison$
        `Pr(>F)`[2],

    Ratio_P_in_stage_adjusted_model =
      coef(
        full_summary
      )[
        "z_ratio",
        "Pr(>|t|)"
      ],

    Profile_similarity_P_in_stage_adjusted_model =
      coef(
        full_summary
      )[
        "z_profile_similarity",
        "Pr(>|t|)"
      ],

    stringsAsFactors = FALSE
  )
}

stage_adjusted_results <- do.call(
  rbind,
  stage_adjusted_results
)

rownames(
  stage_adjusted_results
) <- NULL


# ============================================================
# 13. Stage-level summaries
# ============================================================

stage_summary <- do.call(
  rbind,
  lapply(
    groups,
    function(g) {

      x <- diagnostic[
        diagnostic$group == g,
        ,
        drop = FALSE
      ]

      data.frame(
        Group = g,
        Stage = unname(
          stage_labels[g]
        ),
        N = nrow(x),

        Median_3prime_to_5prime_ratio =
          median(
            x$
              Ratio_3prime_to_5prime
          ),

        IQR_3prime_to_5prime_ratio =
          IQR(
            x$
              Ratio_3prime_to_5prime
          ),

        Median_profile_similarity =
          median(
            x$
              Profile_Spearman_to_median
          ),

        IQR_profile_similarity =
          IQR(
            x$
              Profile_Spearman_to_median
          ),

        Median_PC1 =
          median(
            x$PC1
          ),

        Median_PC2 =
          median(
            x$PC2
          ),

        stringsAsFactors = FALSE
      )
    }
  )
)


# ============================================================
# 14. Leave-one-out sensitivity of pairwise Spearman
#
# Used only to determine whether one sample dominates
# the observed metric-PCA association.
# ============================================================

loo_results <- list()

counter <- 1

for (i in seq_len(
  nrow(pair_specs)
)) {

  x_name <-
    pair_specs$
      Coverage_metric[i]

  y_name <-
    pair_specs$
      PCA_axis[i]

  full_rho <- suppressWarnings(
    cor(
      diagnostic[[x_name]],
      diagnostic[[y_name]],
      method = "spearman"
    )
  )

  loo_rho <- numeric(
    nrow(diagnostic)
  )

  for (j in seq_len(
    nrow(diagnostic)
  )) {

    keep <- seq_len(
      nrow(diagnostic)
    ) != j

    loo_rho[j] <-
      suppressWarnings(
        cor(
          diagnostic[
            keep,
            x_name
          ],
          diagnostic[
            keep,
            y_name
          ],
          method = "spearman"
        )
      )
  }

  tmp <- data.frame(
    Coverage_metric = x_name,
    PCA_axis = y_name,
    Full_Spearman_rho = full_rho,
    LOO_min_rho = min(
      loo_rho,
      na.rm = TRUE
    ),
    LOO_max_rho = max(
      loo_rho,
      na.rm = TRUE
    ),
    Maximum_absolute_LOO_change =
      max(
        abs(
          loo_rho -
            full_rho
        ),
        na.rm = TRUE
      ),
    Sample_causing_maximum_change =
      diagnostic$sampleID[
        which.max(
          abs(
            loo_rho -
              full_rho
          )
        )
      ],
    stringsAsFactors = FALSE
  )

  loo_results[[counter]] <- tmp

  counter <- counter + 1
}

loo_results <- do.call(
  rbind,
  loo_results
)

rownames(
  loo_results
) <- NULL


# ============================================================
# 15. Plot helpers
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
      size = 11
    ),
    axis.title = element_text(
      size = 10
    ),
    axis.text = element_text(
      size = 8.5
    ),
    legend.title = element_text(
      size = 8.5
    ),
    legend.text = element_text(
      size = 7.7
    ),
    legend.key.height = unit(
      0.35,
      "cm"
    ),
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )


get_pair_result <- function(
    metric,
    axis) {

  x <- pairwise_results[
    pairwise_results$
      Coverage_metric == metric &
      pairwise_results$
        PCA_axis == axis,
    ,
    drop = FALSE
  ]

  if (nrow(x) != 1) {
    stop(
      "Pairwise result lookup failed."
    )
  }

  x
}


make_scatter <- function(
    x_var,
    y_var,
    x_label,
    y_label,
    title,
    tag) {

  stat <- get_pair_result(
    x_var,
    y_var
  )

  annotation_text <- paste0(
    "Spearman ρ = ",
    sprintf(
      "%.3f",
      stat$Spearman_rho
    ),
    ", P ",
    ifelse(
      stat$Spearman_P < 0.001,
      "< 0.001",
      paste0(
        "= ",
        sprintf(
          "%.3f",
          stat$Spearman_P
        )
      )
    ),
    "\n",
    "Pearson r = ",
    sprintf(
      "%.3f",
      stat$Pearson_r
    ),
    ", P ",
    ifelse(
      stat$Pearson_P < 0.001,
      "< 0.001",
      paste0(
        "= ",
        sprintf(
          "%.3f",
          stat$Pearson_P
        )
      )
    ),
    "\n",
    "Linear R² = ",
    sprintf(
      "%.3f",
      stat$Linear_model_R2
    )
  )

  ggplot(
    diagnostic,
    aes(
      x = .data[[x_var]],
      y = .data[[y_var]],
      colour = group
    )
  ) +

    geom_point(
      size = 2.5,
      alpha = 0.88
    ) +

    geom_smooth(
      aes(
        group = 1
      ),
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      colour = "black",
      fill = "grey75",
      linewidth = 0.7,
      alpha = 0.25
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

    annotate(
      "label",
      x = -Inf,
      y = Inf,
      label = annotation_text,
      hjust = -0.05,
      vjust = 1.08,
      size = 2.7,
      fill = "white",
      colour = "grey20",
      linewidth = 0.25
    ) +

    labs(
      tag = tag,
      title = title,
      x = x_label,
      y = y_label
    ) +

    theme_pub +

    theme(
      legend.position = "right"
    )
}


# ============================================================
# 16. Four diagnostic panels
# ============================================================

pA <- make_scatter(
  x_var =
    "Ratio_3prime_to_5prime",

  y_var =
    "PC1",

  x_label =
    "3'/5' raw-coverage ratio",

  y_label =
    "PC1 score",

  title =
    "Gene-body 3'/5' coverage ratio versus PC1",

  tag = "A"
)


pB <- make_scatter(
  x_var =
    "Ratio_3prime_to_5prime",

  y_var =
    "PC2",

  x_label =
    "3'/5' raw-coverage ratio",

  y_label =
    "PC2 score",

  title =
    "Gene-body 3'/5' coverage ratio versus PC2",

  tag = "B"
)


pC <- make_scatter(
  x_var =
    "Profile_Spearman_to_median",

  y_var =
    "PC1",

  x_label =
    "Gene-body profile Spearman correlation to cohort median",

  y_label =
    "PC1 score",

  title =
    "Coverage-profile similarity versus PC1",

  tag = "C"
)


pD <- make_scatter(
  x_var =
    "Profile_Spearman_to_median",

  y_var =
    "PC2",

  x_label =
    "Gene-body profile Spearman correlation to cohort median",

  y_label =
    "PC2 score",

  title =
    "Coverage-profile similarity versus PC2",

  tag = "D"
)


# ============================================================
# 17. Save functions
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


save_combined <- function(
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
    "12B_gene_body_coverage_vs_PCA_diagnostic.pdf"
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
      "12B_gene_body_coverage_vs_PCA_",
      "diagnostic_600dpi.tiff"
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


# ============================================================
# 18. Save diagnostic figure
# ============================================================

save_combined(
  pA,
  pB,
  pC,
  pD,
  width = 13,
  height = 10
)


# ============================================================
# 19. Save tables
# ============================================================

write.table(
  diagnostic,
  file = file.path(
    outdir,
    "12B_sample_level_coverage_PCA_data.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  pairwise_results,
  file = file.path(
    outdir,
    "12B_pairwise_coverage_PCA_associations.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  coverage_metric_relationship,
  file = file.path(
    outdir,
    "12B_relationship_between_coverage_metrics.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  coverage_only_models,
  file = file.path(
    outdir,
    "12B_coverage_only_multivariable_models.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  stage_adjusted_results,
  file = file.path(
    outdir,
    "12B_stage_adjusted_incremental_tests.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  stage_summary,
  file = file.path(
    outdir,
    "12B_stage_level_coverage_PCA_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  loo_results,
  file = file.path(
    outdir,
    "12B_leave_one_out_correlation_sensitivity.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 20. Save plot objects / session info
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
    "12B_plot_objects.rds"
  )
)


sink(
  file.path(
    outdir,
    "12B_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 21. Console summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("PAIRWISE ASSOCIATIONS\n")
cat("============================================================\n\n")

print(
  pairwise_results,
  row.names = FALSE,
  digits = 4
)


cat("\n")
cat("============================================================\n")
cat("COVERAGE METRIC RELATIONSHIP\n")
cat("============================================================\n\n")

print(
  coverage_metric_relationship,
  row.names = FALSE,
  digits = 4
)


cat("\n")
cat("============================================================\n")
cat("COVERAGE-ONLY MULTIVARIABLE MODELS\n")
cat("============================================================\n\n")

print(
  coverage_only_models,
  row.names = FALSE,
  digits = 4
)


cat("\n")
cat("============================================================\n")
cat("STAGE-ADJUSTED INCREMENTAL TESTS\n")
cat("============================================================\n\n")

print(
  stage_adjusted_results,
  row.names = FALSE,
  digits = 4
)


cat("\n")
cat("============================================================\n")
cat("LEAVE-ONE-OUT CORRELATION SENSITIVITY\n")
cat("============================================================\n\n")

print(
  loo_results,
  row.names = FALSE,
  digits = 4
)


cat("\n")
cat("============================================================\n")
cat("KEY INTERPRETATION GUIDE\n")
cat("============================================================\n")

cat(
  "1. Pairwise correlations quantify association only; ",
  "they do not establish technical causation.\n"
)

cat(
  "2. Coverage-only model R2 estimates how much PC variation ",
  "is associated with the two coverage metrics together.\n"
)

cat(
  "3. The stage-adjusted incremental test is the critical result.\n"
)

cat(
  "   Reduced model: PC ~ sampling stage\n"
)

cat(
  "   Full model   : PC ~ sampling stage + coverage metrics\n"
)

cat(
  "4. Delta R2 and incremental P indicate whether coverage metrics ",
  "add explanatory information beyond sampling stage.\n"
)

cat(
  "5. No arbitrary pass/fail threshold is applied.\n"
)

cat(
  "6. Leave-one-out results indicate whether a single sample ",
  "dominates the observed correlation.\n"
)


cat("\n")
cat("============================================================\n")
cat("12B DIAGNOSTIC COMPLETED\n")
cat("============================================================\n")

cat(
  "Samples analysed : ",
  nrow(diagnostic),
  "\n",
  sep = ""
)

cat(
  "R2h_1 excluded   : ",
  !("R2h_1" %in%
      diagnostic$sampleID),
  "\n",
  sep = ""
)

cat(
  "R2h_8 retained   : ",
  "R2h_8" %in%
    diagnostic$sampleID,
  "\n",
  sep = ""
)

cat("\nOutput directory:\n")
cat(outdir, "\n\n")

cat("PASS\n")
