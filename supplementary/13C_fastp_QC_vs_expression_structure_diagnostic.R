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
cat("13C FASTP QC VS EXPRESSION-STRUCTURE DIAGNOSTIC\n")
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

full_qc_file <- paste0(
  "supplementary_figures/",
  "Supplementary_Figure_S3_sequencing_mapping_assignment_QC/",
  "13B_fastp_QC/",
  "Supplementary_Table_S2_2026_primary_77sample_",
  "full_sequencing_QC.tsv"
)

outdir <- paste0(
  "supplementary_figures/",
  "Supplementary_Figure_S3_sequencing_mapping_assignment_QC/",
  "13C_fastp_vs_expression_structure"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

RANDOM_SEED <- 20260910


# ============================================================
# 1. Stage definitions
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


# ============================================================
# 2. Validate files
# ============================================================

required_files <- c(
  metadata_file,
  vst_file,
  full_qc_file
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

cat("Required inputs found: PASS\n")


# ============================================================
# 3. Metadata
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
    "Expected 77 metadata samples; observed ",
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
# 4. Read integrated sequencing QC
# ============================================================

qc <- read.delim(
  full_qc_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

required_qc <- c(
  "sampleID",
  "Q30_clean_percent",
  "GC_clean_percent",
  "Duplication_rate_percent"
)

if (!all(
  required_qc %in%
    colnames(qc)
)) {

  stop(
    "Integrated QC table lacks required columns."
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
# 5. Read frozen VST matrix
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
    "Invalid numeric values in VST matrix."
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
# 6. Recompute PCA from frozen VST matrix
# ============================================================

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

cat(
  "PCA variance : PC1 = ",
  sprintf(
    "%.2f%%",
    pca_var[1]
  ),
  "; PC2 = ",
  sprintf(
    "%.2f%%",
    pca_var[2]
  ),
  "\n\n",
  sep = ""
)

if (abs(
  pca_var[1] -
    28.28
) > 0.05 ||
    abs(
      pca_var[2] -
        12.85
    ) > 0.05) {

  warning(
    "PCA variance differs from frozen Figure 2 values."
  )
}


# ============================================================
# 7. Build aligned diagnostic table
# ============================================================

dat <- data.frame(
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

  PC1 =
    pca$x[
      final_samples,
      1
    ],

  PC2 =
    pca$x[
      final_samples,
      2
    ],

  PC3 =
    pca$x[
      final_samples,
      3
    ],

  PC4 =
    pca$x[
      final_samples,
      4
    ],

  PC5 =
    pca$x[
      final_samples,
      5
    ],

  stringsAsFactors = FALSE
)

dat$group <- factor(
  dat$group,
  levels = groups
)

dat$z_Q30 <- as.numeric(
  scale(
    dat$Q30
  )
)

dat$z_GC <- as.numeric(
  scale(
    dat$GC
  )
)

dat$z_Duplication <- as.numeric(
  scale(
    dat$Duplication
  )
)


if (anyNA(
  dat[
    ,
    c(
      "Q30",
      "GC",
      "Duplication",
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
  "Aligned samples : ",
  nrow(dat),
  "\n",
  sep = ""
)

cat(
  "R2h_1 present   : ",
  "R2h_1" %in%
    dat$sampleID,
  "\n",
  sep = ""
)

cat(
  "R2h_8 present   : ",
  "R2h_8" %in%
    dat$sampleID,
  "\n\n",
  sep = ""
)


# ============================================================
# 8. Technical-metric correlations
# ============================================================

tech_matrix <- as.matrix(
  dat[
    ,
    c(
      "Q30",
      "GC",
      "Duplication"
    )
  ]
)

tech_spearman <- cor(
  tech_matrix,
  method = "spearman"
)

tech_pearson <- cor(
  tech_matrix,
  method = "pearson"
)


# ============================================================
# 9. Collinearity diagnostics
#
# VIF-like calculation:
# each technical metric is regressed on sampling stage
# plus the other two technical metrics.
# ============================================================

vif_like <- function(
    response,
    others) {

  formula_text <- paste(
    response,
    "~ group +",
    paste(
      others,
      collapse = " + "
    )
  )

  fit <- lm(
    as.formula(
      formula_text
    ),
    data = dat
  )

  r2 <- summary(
    fit
  )$r.squared

  data.frame(
    Variable = response,
    Auxiliary_R2 = r2,
    VIF_like = 1 / (
      1 -
      r2
    ),
    stringsAsFactors = FALSE
  )
}


vif_table <- rbind(
  vif_like(
    "z_Q30",
    c(
      "z_GC",
      "z_Duplication"
    )
  ),

  vif_like(
    "z_GC",
    c(
      "z_Q30",
      "z_Duplication"
    )
  ),

  vif_like(
    "z_Duplication",
    c(
      "z_Q30",
      "z_GC"
    )
  )
)


design_full <- model.matrix(
  ~ group +
    z_Q30 +
    z_GC +
    z_Duplication,
  data = dat
)

design_rank <- qr(
  design_full
)$rank

design_condition_number <- kappa(
  design_full
)


# ============================================================
# 10. Raw QC-vs-PC associations
# ============================================================

tech_metrics <- c(
  "Q30",
  "GC",
  "Duplication"
)

pc_names <- paste0(
  "PC",
  1:5
)

raw_assoc_list <- list()
counter <- 1

for (metric in tech_metrics) {

  for (pc in pc_names) {

    x <- dat[
      ,
      metric
    ]

    y <- dat[
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

    pe <- suppressWarnings(
      cor.test(
        x,
        y,
        method = "pearson"
      )
    )

    raw_assoc_list[[counter]] <- data.frame(
      Technical_metric = metric,
      PCA_axis = pc,
      N = length(x),

      Spearman_rho =
        unname(
          sp$estimate
        ),

      Spearman_P =
        sp$p.value,

      Pearson_r =
        unname(
          pe$estimate
        ),

      Pearson_P =
        pe$p.value,

      stringsAsFactors = FALSE
    )

    counter <- counter + 1
  }
}

raw_associations <- do.call(
  rbind,
  raw_assoc_list
)

raw_associations$Spearman_BH_FDR <-
  p.adjust(
    raw_associations$Spearman_P,
    method = "BH"
  )


# ============================================================
# 11. Stage-adjusted residual associations
#
# Both QC metric and PC are residualized against stage.
# These correlations therefore quantify within-stage
# association rather than simple between-stage co-variation.
# ============================================================

residual_assoc_list <- list()
counter <- 1

for (metric in tech_metrics) {

  metric_fit <- lm(
    as.formula(
      paste(
        metric,
        "~ group"
      )
    ),
    data = dat
  )

  metric_resid <- residuals(
    metric_fit
  )


  for (pc in pc_names) {

    pc_fit <- lm(
      as.formula(
        paste(
          pc,
          "~ group"
        )
      ),
      data = dat
    )

    pc_resid <- residuals(
      pc_fit
    )

    sp <- suppressWarnings(
      cor.test(
        metric_resid,
        pc_resid,
        method = "spearman",
        exact = FALSE
      )
    )

    pe <- suppressWarnings(
      cor.test(
        metric_resid,
        pc_resid,
        method = "pearson"
      )
    )

    residual_assoc_list[[counter]] <- data.frame(
      Technical_metric = metric,
      PCA_axis = pc,
      N = length(metric_resid),

      Stage_adjusted_Spearman_rho =
        unname(
          sp$estimate
        ),

      Stage_adjusted_Spearman_P =
        sp$p.value,

      Stage_adjusted_Pearson_r =
        unname(
          pe$estimate
        ),

      Stage_adjusted_Pearson_P =
        pe$p.value,

      stringsAsFactors = FALSE
    )

    counter <- counter + 1
  }
}

residual_associations <- do.call(
  rbind,
  residual_assoc_list
)

residual_associations$
  Stage_adjusted_Spearman_BH_FDR <-
  p.adjust(
    residual_associations$
      Stage_adjusted_Spearman_P,
    method = "BH"
  )


# ============================================================
# 12. PC1 / PC2 stage-adjusted incremental models
#
# Reduced:
#   PC ~ stage
#
# Full:
#   PC ~ stage + Q30 + GC + duplication
# ============================================================

incremental_list <- list()

for (axis_name in c(
  "PC1",
  "PC2"
)) {

  reduced_fit <- lm(
    as.formula(
      paste(
        axis_name,
        "~ group"
      )
    ),
    data = dat
  )

  full_fit <- lm(
    as.formula(
      paste(
        axis_name,
        "~ group + z_Q30 + z_GC + z_Duplication"
      )
    ),
    data = dat
  )

  reduced_summary <- summary(
    reduced_fit
  )

  full_summary <- summary(
    full_fit
  )

  cmp <- anova(
    reduced_fit,
    full_fit
  )

  coef_table <- coef(
    full_summary
  )

  incremental_list[[axis_name]] <- data.frame(
    PCA_axis = axis_name,

    Stage_only_R2 =
      reduced_summary$
        r.squared,

    Stage_only_adjusted_R2 =
      reduced_summary$
        adj.r.squared,

    Stage_plus_fastp_R2 =
      full_summary$
        r.squared,

    Stage_plus_fastp_adjusted_R2 =
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
      cmp$F[2],

    Incremental_P =
      cmp$
        `Pr(>F)`[2],

    Q30_standardized_beta =
      coef_table[
        "z_Q30",
        "Estimate"
      ],

    Q30_P =
      coef_table[
        "z_Q30",
        "Pr(>|t|)"
      ],

    GC_standardized_beta =
      coef_table[
        "z_GC",
        "Estimate"
      ],

    GC_P =
      coef_table[
        "z_GC",
        "Pr(>|t|)"
      ],

    Duplication_standardized_beta =
      coef_table[
        "z_Duplication",
        "Estimate"
      ],

    Duplication_P =
      coef_table[
        "z_Duplication",
        "Pr(>|t|)"
      ],

    stringsAsFactors = FALSE
  )
}

incremental_models <- do.call(
  rbind,
  incremental_list
)

rownames(
  incremental_models
) <- NULL


# ============================================================
# 13. Full-expression-space PERMANOVA
# ============================================================

sample_distance <- dist(
  t(vst_mat),
  method = "euclidean"
)


# Stage-only model: reproduce frozen Step 04
set.seed(
  RANDOM_SEED
)

perm_stage <- vegan::adonis2(
  sample_distance ~ group,
  data = dat,
  permutations = 9999,
  by = "margin"
)


# Full model: stage + fastp technical metrics
set.seed(
  RANDOM_SEED
)

perm_full <- vegan::adonis2(
  sample_distance ~
    group +
    z_Q30 +
    z_GC +
    z_Duplication,
  data = dat,
  permutations = 9999,
  by = "margin"
)


stage_only_summary <- data.frame(
  Model =
    "Stage_only",

  Term =
    "group",

  R2 =
    perm_stage[
      "group",
      "R2"
    ],

  F =
    perm_stage[
      "group",
      "F"
    ],

  P =
    perm_stage[
      "group",
      "Pr(>F)"
    ],

  stringsAsFactors = FALSE
)


full_perm_terms <- c(
  "group",
  "z_Q30",
  "z_GC",
  "z_Duplication"
)

full_permanova_summary <- data.frame(
  Model =
    "Stage_plus_fastp",

  Term =
    full_perm_terms,

  R2 =
    perm_full[
      full_perm_terms,
      "R2"
    ],

  F =
    perm_full[
      full_perm_terms,
      "F"
    ],

  P =
    perm_full[
      full_perm_terms,
      "Pr(>F)"
    ],

  stringsAsFactors = FALSE
)

full_permanova_summary$BH_FDR <-
  p.adjust(
    full_permanova_summary$P,
    method = "BH"
  )


# ============================================================
# 14. Stage PERMDISP check
#
# Match frozen Step 04 definition exactly.
# ============================================================

dispersion <- vegan::betadisper(
  sample_distance,
  group = dat$group,
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

permdisp_summary <- data.frame(
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


# ============================================================
# 15. Prepare heatmap data
# ============================================================

raw_heatmap <- raw_associations[
  ,
  c(
    "Technical_metric",
    "PCA_axis",
    "Spearman_rho"
  )
]

colnames(
  raw_heatmap
)[3] <- "rho"

raw_heatmap$Type <-
  "Raw association"


adjusted_heatmap <- residual_associations[
  ,
  c(
    "Technical_metric",
    "PCA_axis",
    "Stage_adjusted_Spearman_rho"
  )
]

colnames(
  adjusted_heatmap
)[3] <- "rho"

adjusted_heatmap$Type <-
  "Stage-adjusted"


raw_heatmap$Technical_metric <- factor(
  raw_heatmap$Technical_metric,
  levels = tech_metrics,
  labels = c(
    "Q30",
    "GC",
    "Duplication"
  )
)

adjusted_heatmap$Technical_metric <- factor(
  adjusted_heatmap$Technical_metric,
  levels = tech_metrics,
  labels = c(
    "Q30",
    "GC",
    "Duplication"
  )
)

raw_heatmap$PCA_axis <- factor(
  raw_heatmap$PCA_axis,
  levels = pc_names
)

adjusted_heatmap$PCA_axis <- factor(
  adjusted_heatmap$PCA_axis,
  levels = pc_names
)


# ============================================================
# 16. Publication-style plotting theme
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
      size = 8.3,
      colour = "grey35"
    ),
    axis.title = element_text(
      size = 9.5
    ),
    axis.text = element_text(
      size = 8
    ),
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )


# ============================================================
# 17. Panel A: raw QC-PC associations
# ============================================================

pA <- ggplot(
  raw_heatmap,
  aes(
    x = PCA_axis,
    y = Technical_metric,
    fill = rho
  )
) +

  geom_tile(
    colour = "white",
    linewidth = 0.7
  ) +

  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        rho
      )
    ),
    size = 3.2
  ) +

  scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "white",
    high = "#C9504D",
    midpoint = 0,
    limits = c(
      -1,
      1
    ),
    name = "Spearman\nrho"
  ) +

  labs(
    tag = "A",
    title =
      "Raw associations between fastp metrics and principal components",
    subtitle =
      "Signs of PCA axes are arbitrary; association magnitude is the primary diagnostic",
    x = NULL,
    y = NULL
  ) +

  theme_pub


# ============================================================
# 18. Panel B: stage-adjusted associations
# ============================================================

pB <- ggplot(
  adjusted_heatmap,
  aes(
    x = PCA_axis,
    y = Technical_metric,
    fill = rho
  )
) +

  geom_tile(
    colour = "white",
    linewidth = 0.7
  ) +

  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        rho
      )
    ),
    size = 3.2
  ) +

  scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "white",
    high = "#C9504D",
    midpoint = 0,
    limits = c(
      -1,
      1
    ),
    name = "Residual\nSpearman rho"
  ) +

  labs(
    tag = "B",
    title =
      "Associations after adjustment for sampling stage",
    subtitle =
      "Both QC metrics and PCs were residualized against sampling stage",
    x = NULL,
    y = NULL
  ) +

  theme_pub


# ============================================================
# 19. Panel C: incremental R2 for PC1 / PC2
# ============================================================

r2_plot <- rbind(
  data.frame(
    PCA_axis =
      incremental_models$
        PCA_axis,
    Model =
      "Stage only",
    R2 =
      incremental_models$
        Stage_only_R2,
    stringsAsFactors = FALSE
  ),

  data.frame(
    PCA_axis =
      incremental_models$
        PCA_axis,
    Model =
      "Stage + fastp QC",
    R2 =
      incremental_models$
        Stage_plus_fastp_R2,
    stringsAsFactors = FALSE
  )
)

r2_plot$Model <- factor(
  r2_plot$Model,
  levels = c(
    "Stage only",
    "Stage + fastp QC"
  )
)


pC <- ggplot(
  r2_plot,
  aes(
    x = PCA_axis,
    y = R2,
    fill = Model
  )
) +

  geom_col(
    position = position_dodge(
      width = 0.72
    ),
    width = 0.62,
    colour = "grey25",
    linewidth = 0.4
  ) +

  geom_text(
    aes(
      label = sprintf(
        "%.3f",
        R2
      )
    ),
    position = position_dodge(
      width = 0.72
    ),
    vjust = -0.4,
    size = 3
  ) +

  scale_fill_manual(
    values = c(
      "grey75",
      "grey40"
    ),
    name = NULL
  ) +

  scale_y_continuous(
    limits = c(
      0,
      max(
        r2_plot$R2
      ) * 1.16
    ),
    expand = c(
      0,
      0
    )
  ) +

  labs(
    tag = "C",
    title =
      "Incremental explanatory value of fastp QC metrics",
    subtitle =
      "PC ~ stage versus PC ~ stage + Q30 + GC + duplication",
    x = NULL,
    y = "Model R²"
  ) +

  theme_pub +

  theme(
    legend.position = "top"
  )


# ============================================================
# 20. Panel D: PERMANOVA marginal R2
# ============================================================

perm_plot <- full_permanova_summary

perm_plot$Term_label <- factor(
  perm_plot$Term,
  levels = c(
    "group",
    "z_Q30",
    "z_GC",
    "z_Duplication"
  ),
  labels = c(
    "Sampling stage",
    "Q30",
    "GC",
    "Duplication"
  )
)


pD <- ggplot(
  perm_plot,
  aes(
    x = Term_label,
    y = R2
  )
) +

  geom_col(
    width = 0.62,
    fill = "grey65",
    colour = "grey20",
    linewidth = 0.45
  ) +

  geom_text(
    aes(
      label = paste0(
        "R²=",
        sprintf(
          "%.3f",
          R2
        ),
        "\nP",
        ifelse(
          P < 0.001,
          "<0.001",
          paste0(
            "=",
            sprintf(
              "%.3f",
              P
            )
          )
        )
      )
    ),
    vjust = -0.35,
    size = 2.8,
    lineheight = 0.95
  ) +

  scale_y_continuous(
    limits = c(
      0,
      max(
        perm_plot$R2
      ) * 1.25
    ),
    expand = c(
      0,
      0
    )
  ) +

  labs(
    tag = "D",
    title =
      "Unique contributions in full VST expression space",
    subtitle =
      "Marginal PERMANOVA: stage + Q30 + GC + duplication",
    x = NULL,
    y = "Marginal PERMANOVA R²"
  ) +

  theme_pub +

  theme(
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    )
  )


# ============================================================
# 21. Save combined figure
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
  "13C_fastp_QC_vs_expression_structure.pdf"
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
  "13C_fastp_QC_vs_expression_structure_600dpi.tiff"
)

open_tiff(
  tiff_file,
  13,
  10
)

draw_combined()
dev.off()


# ============================================================
# 22. Save all tables
# ============================================================

write.table(
  dat,
  file = file.path(
    outdir,
    "13C_sample_level_fastp_PCA_data.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  tech_spearman,
  file = file.path(
    outdir,
    "13C_fastp_metric_Spearman_matrix.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


write.table(
  tech_pearson,
  file = file.path(
    outdir,
    "13C_fastp_metric_Pearson_matrix.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


write.table(
  vif_table,
  file = file.path(
    outdir,
    "13C_fastp_collinearity_VIF_like.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  raw_associations,
  file = file.path(
    outdir,
    "13C_raw_fastp_PC_associations.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  residual_associations,
  file = file.path(
    outdir,
    "13C_stage_adjusted_fastp_PC_associations.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  incremental_models,
  file = file.path(
    outdir,
    "13C_stage_adjusted_incremental_PC_models.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  stage_only_summary,
  file = file.path(
    outdir,
    "13C_stage_only_PERMANOVA.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  full_permanova_summary,
  file = file.path(
    outdir,
    "13C_stage_plus_fastp_marginal_PERMANOVA.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  permdisp_summary,
  file = file.path(
    outdir,
    "13C_PERMDISP_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 23. Session info
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
    "13C_plot_objects.rds"
  )
)

sink(
  file.path(
    outdir,
    "13C_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 24. Console report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("TECHNICAL-METRIC SPEARMAN CORRELATION\n")
cat("============================================================\n\n")

print(
  tech_spearman,
  digits = 4
)


cat("\n")
cat("============================================================\n")
cat("COLLINEARITY DIAGNOSTICS\n")
cat("============================================================\n\n")

print(
  vif_table,
  row.names = FALSE,
  digits = 4
)

cat(
  "\nFull design rank : ",
  design_rank,
  "/",
  ncol(design_full),
  "\n",
  sep = ""
)

cat(
  "Design condition number : ",
  sprintf(
    "%.4f",
    design_condition_number
  ),
  "\n",
  sep = ""
)


cat("\n")
cat("============================================================\n")
cat("RAW FASTP-PC ASSOCIATIONS\n")
cat("============================================================\n\n")

print(
  raw_associations,
  row.names = FALSE,
  digits = 4
)


cat("\n")
cat("============================================================\n")
cat("STAGE-ADJUSTED FASTP-PC ASSOCIATIONS\n")
cat("============================================================\n\n")

print(
  residual_associations,
  row.names = FALSE,
  digits = 4
)


cat("\n")
cat("============================================================\n")
cat("STAGE-ADJUSTED INCREMENTAL PC MODELS\n")
cat("============================================================\n\n")

print(
  incremental_models,
  row.names = FALSE,
  digits = 5
)


cat("\n")
cat("============================================================\n")
cat("STAGE-ONLY PERMANOVA\n")
cat("============================================================\n\n")

print(
  stage_only_summary,
  row.names = FALSE,
  digits = 5
)


cat("\n")
cat("============================================================\n")
cat("STAGE + FASTP MARGINAL PERMANOVA\n")
cat("============================================================\n\n")

print(
  full_permanova_summary,
  row.names = FALSE,
  digits = 5
)


cat("\n")
cat("============================================================\n")
cat("PERMDISP\n")
cat("============================================================\n\n")

print(
  permdisp_summary,
  row.names = FALSE,
  digits = 5
)


cat("\n")
cat("============================================================\n")
cat("13C COMPLETED\n")
cat("============================================================\n")

cat(
  "Samples analysed : ",
  nrow(dat),
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
      dat$sampleID),
  "\n",
  sep = ""
)

cat(
  "R2h_8 retained   : ",
  "R2h_8" %in%
    dat$sampleID,
  "\n",
  sep = ""
)

cat(
  "Original PC1     : ",
  sprintf(
    "%.2f%%",
    pca_var[1]
  ),
  "\n",
  sep = ""
)

cat(
  "Original PC2     : ",
  sprintf(
    "%.2f%%",
    pca_var[2]
  ),
  "\n",
  sep = ""
)

cat("\nOutput directory:\n")
cat(outdir, "\n\n")

cat("PASS\n")
